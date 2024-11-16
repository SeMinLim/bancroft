import FIFO::*;
import FIFOF::*;
import Vector::*;

import Serializer::*;


typedef 64 KmerLength;
typedef struct {
	Bit#(64) addr;
	Bit#(32) bytes;
} MemPortReq deriving (Eq,Bits);


interface DecompressorIfc;
	method Action start(Bit#(32) param);
	method ActionValue#(Bool) done;
endinterface
module mkDecompressor( DecompressorIfc );
	Reg#(Bool) started <- mkReg(False);
	FIFO#(Bool) doneQ <- mkFIFO;

	SerializerIfc#(512, 16) serializer512b32b <- mkSerializer;
	//------------------------------------------------------------------------------------
	// Cycle Counter
	//------------------------------------------------------------------------------------
	Reg#(Bit#(32)) cycleCounter <- mkReg(0);
	rule incCycle;
		cycleCounter <= cycleCounter + 1;
	endrule
	//------------------------------------------------------------------------------------
	// Memory Read & Write
	//------------------------------------------------------------------------------------
	Vector#(MemPortCnt, FIFO#(MemPortReq)) readReqQs <- replicateM(mkSizedFIFO(1024));
	Vector#(MemPortCnt, FIFO#(Bit#(512))) readWordQs <- replicateM(mkSizedFIFO(1024));
	FIFOF#(Bit#(32)) packet32bQ <- mkSizedFIFOF(1024);
	// Read compressed sequence
	Reg#(Bit#(64)) pcktAddrOff <- mkReg(0);
	rule sendReadReqPckt;
		readReqQs[0].enq(MemPortReq{addr:pcktAddrOff, bytes:64});
		pcktAddrOff <= pcktAddrOff + 64;
	endrule
	rule readWordPckt;
		readWordQs[0].deq;
		let r = readWordQs[0].first;
		serializer512b32b.put(r);
	endrule
	rule getWordPckt;
		let r <- serializer512b32b.get;
		packet32bQ.enq(r);
	endrule
	// Read reference sequence
	FIFO#(MemPortReq) readReqRefQ <- mkSizedFIFO(1024);
	FIFO#(Bit#(512)) readWordRefQ <- mkSizedFIFO(1024);
	Reg#(Bit#(64)) addrBuf <- mkReg(0);
	Reg#(Bit#(32)) bytesBuf <- mkReg(0);
	Reg#(Bool) readReqRefDEQUEUE <- mkReg(True);
	rule sendReadReqRef;
		Bit#(64) addr = 0;
		Bit#(32) bytes = 0;
		if ( readReqRefDEQUEUE ) begin
			readReqRefQ.deq;
			let r = readReqRefQ.first;
			
			addr = r.addr;
			bytes = r.bytes;

			if ( bytes > 64 ) begin
				addrBuf <= addr + 64;
				bytesBuf <= bytes - 64;
				readReqRefDEQUEUE <= False;
			end
		end else begin
			addr = addrBuf;
			bytes = bytesBuf;

			if ( bytes == 64 ) begin
				addrBuf <= 0;
				bytesBuf <= 0;
				readReqRefDEQUEUE <= True;
			end else begin
				addrBuf <= addr + 64;
				bytesBuf <= bytes - 64;
			end
		end

		readReqQs[1].enq(MemPortReq{addr:addr, bytes:64});
	endrule
	rule readWordRef;
		readWordQs[1].deq;
		let r = readWordQs[1].first;
		readWordRefQ.enq(r);
	endrule
	//------------------------------------------------------------------------------------
	// Decompressor
	//------------------------------------------------------------------------------------
	Reg#(Bool) headNeeded <- mkReg(True);
	Vector#(4, FIFO#(Bit#(8))) relayHeaderQs <- replicateM(mkFIFO);
	rule relayHeader( headNeeded );
		packet32bQ.deq;
		for ( Integer i = 0; i < 4; i = i + 1 ) begin
			relayHeaderQs[i].enq(truncate(packet32bQ.first >> (i * 8)));
		end
		headNeeded <= False;
	endrule
	
	Vector#(16, Bit#(8)) sendPcktReqQs <- replicateM(mkFIFO);
	Vector#(4, Bit#(32)) pcktCntQs <- replicateM(mkFIFO);
	for ( Integer i = 0; i < 4; i = i + 1 ) begin
		rule interpretHeader;
			relayHeaderQs[i].deq;
			let heads = relayHeaderQs[i].first;
			
			let pcktCnt = 0;
			for ( Integer j = 0; j < 4; j = j + 1 ) begin
				Bit#(2) head = truncate(heads >> (j * 2));
				if ( head == 2'b00 ) begin 	    	// [UNMATCH]
					sendPcktReqQs[(i*4) + j].enq(0);
					pcktCnt = pcktCnt + 1;
				end else if ( head == 2'b01 ) begin 	// [FORWARD NORMAL MATCH]
					sendPcktReqQs[(i*4) + j].enq(1);
					pcktCnt = pcktCnt + 1;
				end else if ( head == 2'b10 ) begin 	// [CONTINUOUS MATCH]
					sendPcktReqQs[(i*4) + j].enq(2);
				end else if ( head == 2'b11 ) begin 	// [REVERSE NORMAL MATCH]
					sendPcktReqQs[(i*4) + j].enq(3);
					pcktCnt = pcktCnt + 1;
				end
			end

			pcktCntQs[i].enq(pcktCnt);
		endrule
	end
	
	Reg#(Bit#(32)) pcktCntBuf <- mkReg(0);
	rule calPcktCnt;
		Vector#(4, Bit#(32)) pcktCntEach = replicate(0);
		for ( Integer i = 0; i < 4; i = i + 1 ) begin
			pcktCntQs[i].deq;
			pcktCntEach[i] = pcktCntQs.first;
		end
		Bit#(32) pcktCntTotal = pcktCntEach[0] + pcktCntEach[1] + pcktCntEach[2] + pcktCntEach[3];
		pcktCntBuf <= pcktCntTotal;
	endrule
	
	FIFO#(Bit#(8)) orderQ <- mkFIFO;
	FIFO#(Bit#(32)) verbatimQ <- mkFIFO;
	FIFO#(Bit#(32)) bytesNeededQ <- mkFIFO;
	Reg#(Bit#(32)) idxStartBuf <- mkReg(0);
	Reg#(Bit#(32)) idxEndBuf <- mkReg(0);
	Reg#(Bit#(8)) idxDirection <- mkReg(0);
	Reg#(Bit#(32)) readPcktCnt <- mkReg(0);
	Reg#(Bit#(32)) pcktCntTmpBuf <- mkReg(0);
	rule readPckt;
		sendPcktReqQs[readPcktCnt].deq;
		let r = sendPcktReqQs[readPcktCnt].first;

		Bit#(32) pcktCnt = pcktCntTmpBuf;
		if ( r == 0 ) begin				// [UNMATCH]
			if ( idxDirection != 0 ) begin
				Bit#(64) addr 	 = (zeroExtend(idxStartBuf) * 2) / 8;
				Bit#(64) idxPntr = (zeroExtend(idxStartBuf) * 2) % 8;

				Bit#(32) totalLength 	= (idxEndBuf - idxStartBuf) * 2;
				Bit#(32) bytesNeededTmp = (totalLength / 512) + 1;
				Bit#(32) bytesAddtnl 	= 0;
				Bit#(32) bytesNeeded 	= 0;
				if ( idxEndBuf > idxStartBuf ) begin
					bytesAddtnl 	= ((totalLength % 512) + idxPntr) / 512;
					bytesNeeded 	= (bytesNeededTmp + bytesAddtnl) * 64;
				end else begin
					bytesNeeded 	= 64;
				end
				
				readReqRefQ.enq(MemPortReq{addr:addr, bytes:bytesNeeded});
				bytesNeededQ.enq(bytesNeeded);
				
				idxDirection <= 0;

				orderQ.enq(1);
			end else begin
				orderQ.enq(0);
			end

			packet32bQ.deq;
			verbatimQ.enq(packet32bQ.first);
			
			pcktCnt = pcktCnt + 1;
		end else if ( r == 1 ) begin			// [FORWARD NORMAL MATCH]
			if ( idxDirection != 0 ) begin
				Bit#(64) addr 	 = (zeroExtend(idxStartBuf) * 2) / 8;
				Bit#(64) idxPntr = (zeroExtend(idxStartBuf) * 2) % 8;

				Bit#(32) totalLength 	= (idxEndBuf - idxStartBuf) * 2;
				Bit#(32) bytesNeededTmp = (totalLength / 512) + 1;
				Bit#(32) bytesAddtnl 	= 0;
				Bit#(32) bytesNeeded 	= 0;
				if ( idxEndBuf > idxStartBuf ) begin
					bytesAddtnl 	= ((totalLength % 512) + idxPntr) / 512;
					bytesNeeded 	= (bytesNeededTmp + bytesAddtnl) * 64;
				end else begin
					bytesNeeded 	= 64;
				end
				
				readReqRefQ.enq(MemPortReq{addr:addr, bytes:bytesNeeded});
				bytesNeededQ.enq(bytesNeeded);

				orderQ.enq(2);
			end

			packet32bQ.deq;
			
			idxStartBuf <= packet32bQ.first;
			idxEndBuf <= packet32bQ.first;
			idxDirection <= 1;
			
			pcktCnt = pcktCnt + 1;
		end else if ( r == 2 ) begin 			// [CONTINUOUS MATCH]
			if ( idxDirection == 1 ) begin		// [FORWARD]
				idxEndBuf <= idxEndBuf + fromInteger(valueOf(KmerLength));
			end else if ( idxDirection == 2 ) begin	// [REVERSE]
				idxStartBuf <= idxStartBuf - fromInteger(valueOf(KmerLength));
			end
		end else if ( r == 3 ) begin			// [REVERSE NORMAL MATCH]
			if ( idxDirection != 0 ) begin
				Bit#(64) addr 	 = (zeroExtend(idxStartBuf) * 2) / 8;
				Bit#(64) idxPntr = (zeroExtend(idxStartBuf) * 2) % 8;

				Bit#(32) totalLength 	= (idxEndBuf - idxStartBuf) * 2;
				Bit#(32) bytesNeededTmp = (totalLength / 512) + 1;
				Bit#(32) bytesAddtnl 	= 0;
				Bit#(32) bytesNeeded 	= 0;
				if ( idxEndBuf > idxStartBuf ) begin
					bytesAddtnl 	= ((totalLength % 512) + idxPntr) / 512;
					bytesNeeded 	= (bytesNeededTmp + bytesAddtnl) * 64;
				end else begin
					bytesNeeded 	= 64;
				end
				
				readReqRefQ.enq(MemPortReq{addr:addr, bytes:bytesNeeded});
				bytesNeededQ.enq(bytesNeeded);

				orderQ.enq(2);
			end

			packet32bQ.deq;
			
			idxStartBuf <= packet32bQ.first;
			idxEndBuf <= packet32bQ.first;

			idxDirection <= 2;
			
			pcktCnt = pcktCnt + 1;
		end

		if ( readPcktCnt + 1 == 16 ) begin
			if ( pcktCnt != pcktCntBuf ) begin
				$write( "\033[0;33m[%1d]\033[0m -> \033[1;32m[Decompressor]\033[0m SYSTEM ERROR [TYPE 1]\n", cycleCount );
			end
			readPcktCnt <= 0;
			headerNeeded <= True;
		end else begin
			readPcktCnt <= readPcktCnt + 1;
			pcktCntTmpBuf <= pcktCnt;
		end
	endrule

	Reg#(Bit#(32)) decodeCnt <- mkReg(0);
	Reg#(Bit#(32)) decodeResult <- mkReg(0);
	Reg#(Bit#(32)) bytesNeededBuf <- mkReg(0);
	Reg#(Bit#(8)) orderBuf <- mkReg(0);
	Reg#(Bool) orderDEQUEUE <- mkReg(True);
	rule decode;
		Bit#(8) order = 0;
		if ( orderDEQUEUE ) begin
			orderQ.deq;
			order = orderQ.first;
		end else begin
			order = orderBuf;
		end

		if ( order == 0 ) begin			// [UNMATCH FIRST]
			verbatimQ.deq;
			decodeResult <= decodeResult ^ verbatimQ.first;
			decodeCnt <= decodeCnt + 1;
		end else begin				// [MATCH FIRST]
			Bit#(32) bytesNeeded = 0;
			if ( orderDEQUEUE ) begin
				bytesNeededQ.deq;
				bytesNeeded = bytesNeededQ.first;
			end else begin
				bytesNeeded = bytesNeededBuf;
			end

			readWordRefQ.deq;
			Bit#(32) decodeResultTmp = decodeResult ^ truncate(readWordRefQ.first);
			
			if ( bytesNeeded > 64 ) begin
				decodeResult <= decodeResultTmp;
				orderDEQUEUE <= False;
				bytesNeddedBuf <= byesNeeded - 64;
			end else begin
				if ( order == 1 ) begin // [MATCH THEN UNMATCH]
					verbatimQ.deq;
					decodeResultTmp = decodeResultTmp ^ verbatimQ.first;
					decodeResult <= decodeResultTmp;
					decodeCnt <= decodeCnt + 2;
				end else begin		// [MATCH THEN ANOTHER MATCH]
					decodeResult <= decodeResultTmp;
					decodeCnt <= decodeCnt + 1;
				end
				orderDEQUEUE <= True;
			end
		end

		if ( decodeCnt + 1 == pcktCntBuf ) doneQ.enq(1);
	endrule
	//------------------------------------------------------------------------------------
	// Interface
	//------------------------------------------------------------------------------------
	// Read the compressed genomic data
	method Action readPckt(Bit#(512) pckt);
		pcktQ.enq(pckt);
	endmethod
	// Read the 2-bit encoded reference genome
	method ActionValue#(MemPortReq) rqstRdRefr;
		rqstRdRefrQ.deq;
		return rqstRdRefrQ.first;
	endmethod
	method Action readRefr(Bit#(512) refr);
		refrQ.enq(refr);
	endmethod
	// Write the decompressed 2-bit encoded sequence
	method ActionValue#(MemPortReq) rqstWrRslt;
		rqstWrRsltQ.deq;
		return rqstWrRsltQ.first;
	endmethod
	method ActionValue#(Bit#(512)) rslt;
		rsltQ.deq;
		return rsltQ.first;
	endmethod
endmodule

