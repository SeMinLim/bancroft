import FIFO::*;
import FIFOF::*;
import Vector::*;

import BRAM::*;
import BRAMFIFO::*;

import Serializer::*;


typedef 16	   Head;
typedef 64 	   Kmer;
typedef 3716988	   PcktCntHead32b;
typedef 25394544   PcktCntBody32b;
typedef 29111532   PcktCntTotal32b;
typedef 1819471    PcktCntTotal512b;
typedef 17990661   RsltCntTotal32b;
typedef 41481147   RsltCntTotal128b;
typedef 11494704   RsltCntTotal512b;
typedef 5885287968 RsltCntTotalLngt;
typedef TMul#(PcktCntHead32b, 16) TotalTrial;
typedef struct {
	Bit#(64) addr;
	Bit#(32) bytes;
} MemPortReq deriving (Eq,Bits);


function Tuple3#(Bit#(64), Bit#(32), Bit#(64)) calcPrmtRqst( Bit#(32) idxStrt, Bit#(32) idxCntn );
	Bit#(64) address = (zeroExtend(idxStrt) * 2) / 8;
	Bit#(64) idxPntr = (zeroExtend(idxStrt) * 2) % 8;

	Bit#(32) readLngtByte 	  = 0;
	Bit#(32) readLngtByteTmpl = ((idxCntn + 1) / 4);
	if ( (idxCntn + 1) % 4  == 0 ) begin
		if ( idxPntr > 0 ) readLngtByte = (readLngtByteTmpl + 1) * 64;
		else 		   readLngtByte = readLngtByteTmpl * 64;
	end else begin
		readLngtByte = (readLngtByteTmpl + 1) * 64;
	end
	
	return tuple3(address, readLngtByte, idxPntr);
endfunction


interface DecompressorIfc;
	method Action readPckt(Bit#(512) pckt);
	method ActionValue#(MemPortReq) rqstRdRefr;
	method Action readRefr(Bit#(512) refr);
	method ActionValue#(MemPortReq) rqstWrRslt;
	method ActionValue#(Bit#(512)) wrteRslt;
endinterface
(* synthesize *)
module mkDecompressor( DecompressorIfc );
	SerializerIfc#(512, 16) serializer512b32b <- mkSerializer;

	FIFO#(MemPortReq) rqstRdRefrQ <- mkSizedBRAMFIFO(1024);
	FIFO#(MemPortReq) rqstWrRsltQ <- mkSizedBRAMFIFO(1024);
	FIFO#(Bit#(32))  pcktQ <- mkSizedBRAMFIFO(1024);
	FIFO#(Bit#(32))  vbtmQ <- mkSizedBRAMFIFO(1024);
	FIFO#(Bit#(512)) refrQ <- mkSizedBRAMFIFO(1024);
	FIFO#(Bit#(512)) rsltQ <- mkSizedBRAMFIFO(1024);

	Reg#(Bit#(32)) pcktHeadTrck <- mkReg(0);
	Reg#(Bit#(32)) pcktBodyTrck <- mkReg(0);
	
	Reg#(Bool) readHeadOn <- mkReg(True);
	Reg#(Bool) intpHeadOn <- mkReg(False);
	Reg#(Bool) cntlOrdrOn <- mkReg(False);
	Reg#(Bool) dcomprssOn <- mkReg(False);
	//------------------------------------------------------------------------------------
	// Cycle Counter
	//------------------------------------------------------------------------------------
	Reg#(Bit#(32)) cycleCounter <- mkReg(0);
	rule incCycle;
		cycleCounter <= cycleCounter + 1;
	endrule
	//------------------------------------------------------------------------------------
	// Serializer
	//------------------------------------------------------------------------------------
	rule serialize;
		let pckt <- serializer512b32b.get;
		pcktQ.enq(pckt);
	endrule
	//------------------------------------------------------------------------------------
	// [STAGE 1]
	// Interpret the head section first
	//------------------------------------------------------------------------------------
	FIFO#(Bit#(2)) headQ <- mkFIFO;
	Reg#(Bit#(32)) headBuf 	   <- mkReg(0);
	Reg#(Bit#(32)) readHeadCnt <- mkReg(0);
	rule readHead( readHeadOn );
		if ( readHeadCnt == 0 ) begin
			pcktQ.deq;
			headQ.enq(truncate(pcktQ.first));

			headBuf     <= pcktQ.first >> 2;
			readHeadCnt <= readHeadCnt + 1;
			intpHeadOn  <= True;
		end else begin
			headQ.enq(truncate(headBuf));

			if ( readHeadCnt + 1 == fromInteger(valueOf(Head)) ) begin
				if ( pcktHeadTrck + 1 == fromInteger(valueOf(PcktCntHead32b)) ) begin
					headBuf      <= 1;
				end else begin
					headBuf      <= 0;
					readHeadCnt  <= 0;
					readHeadOn   <= False;
				end
				pcktHeadTrck <= pcktHeadTrck + 1;
			end else begin
				headBuf      <= headBuf >> 2;
				readHeadCnt  <= readHeadCnt + 1;
			end
		end
	endrule

	FIFO#(Tuple5#(Bit#(32), Bit#(64), Bit#(32), Bool, Bit#(8))) needPrmtDecdQ <- mkSizedBRAMFIFO(1024);
	FIFO#(Bit#(8))	ordrQ	      	<- mkSizedBRAMFIFO(1024);
	Reg#(Bit#(8))   idxDrctBuf    	<- mkReg(0);
	Reg#(Bit#(32))  idxStrtBuf    	<- mkReg(0);
	Reg#(Bit#(32))  idxCntnBuf    	<- mkReg(0);
	Reg#(Bit#(32))  intpHeadCurrCnt <- mkReg(0);
	Reg#(Bit#(32))  intpHeadTotlCnt <- mkReg(0);
	rule intpHead( intpHeadOn );
		headQ.deq;
		let head = headQ.first;
		
		if ( head == 2'b00 ) begin 	    // [UNMATCH]
			if ( idxDrctBuf != 0 ) begin
				let p = calcPrmtRqst(idxStrtBuf, idxCntnBuf);	

				rqstRdRefrQ.enq(MemPortReq{addr:tpl_1(p), bytes:tpl_2(p)});
				needPrmtDecdQ.enq(tuple5(tpl_2(p), tpl_3(p), idxCntnBuf, True, idxDrctBuf));
				ordrQ.enq(1);
			
				idxDrctBuf <= 0;
				idxStrtBuf <= 0;
				idxCntnBuf <= 0;
			end else begin
				ordrQ.enq(0);
			end

			pcktQ.deq;
			vbtmQ.enq(pcktQ.first);
		end else if ( head == 2'b01 ) begin // [FORWARD NORMAL MATCH]
			if ( idxDrctBuf != 0 ) begin
				let p = calcPrmtRqst(idxStrtBuf, idxCntnBuf);

				rqstRdRefrQ.enq(MemPortReq{addr:tpl_1(p), bytes:tpl_2(p)});
				needPrmtDecdQ.enq(tuple5(tpl_2(p), tpl_3(p), idxCntnBuf, True, idxDrctBuf));
				ordrQ.enq(2);
			end

			pcktQ.deq;
			idxStrtBuf <= pcktQ.first;
			idxDrctBuf <= 1;
			idxCntnBuf <= 1;
		end else if ( head == 2'b10 ) begin // [CONTINUOUS MATCH]
			if ( idxDrctBuf == 1 ) begin		// [FORWARD]
				idxCntnBuf <= idxCntnBuf + 1;
			end else if ( idxDrctBuf == 2 ) begin	// [REVERSE]
				idxStrtBuf <= idxStrtBuf - fromInteger(valueOf(Kmer));
				idxCntnBuf <= idxCntnBuf + 1;
			end
		end else if ( head == 2'b11 ) begin // [REVERSE NORMAL MATCH]
			if ( idxDrctBuf != 0 ) begin
				let p = calcPrmtRqst(idxStrtBuf, idxCntnBuf);

				rqstRdRefrQ.enq(MemPortReq{addr:tpl_1(p), bytes:tpl_2(p)});
				needPrmtDecdQ.enq(tuple5(tpl_2(p), tpl_3(p), idxCntnBuf, True, idxDrctBuf));
				ordrQ.enq(2);
			end

			pcktQ.deq;
			idxStrtBuf <= pcktQ.first;
			idxDrctBuf <= 2;
			idxCntnBuf <= 1;
		end

		if ( intpHeadCurrCnt + 1 == fromInteger(valueOf(Head)) ) begin
			if ( intpHeadTotlCnt + 1 == fromInteger(valueOf(TotalTrial)) ) begin
				intpHeadCurrCnt <= intpHeadCurrCnt;
			end else begin
				readHeadOn 	<= True;
				intpHeadOn 	<= False;
				intpHeadCurrCnt <= 0;
			end
			intpHeadTotlCnt <= intpHeadTotlCnt + 1;
		end else begin
			cntlOrdrOn 	<= True;
			intpHeadCurrCnt <= intpHeadCurrCnt + 1;
			intpHeadTotlCnt <= intpHeadTotlCnt + 1;
		end
	endrule
	//------------------------------------------------------------------------------------
	// [STAGE 2]
	// Decompress the sequence & Send the decompressed data out as 512-bit
	//------------------------------------------------------------------------------------
	FIFO#(Tuple5#(Bit#(32), Bit#(64), Bit#(32), Bool, Bit#(8))) dcomprssPrmtQ <- mkFIFO;
	FIFO#(Tuple4#(Bit#(512), Bit#(32), Bit#(32), Bit#(8)))  dcomprssPipeQ <- mkSizedBRAMFIFO(1024);
	FIFO#(Bit#(8)) dcomprssOrdrQ <- mkFIFO;
	rule cntlOrdr( cntlOrdrOn );
		ordrQ.deq;
		let ordr = ordrQ.first;
		
		if ( ordr != 0 ) begin
			needPrmtDecdQ.deq;
			let prmt = needPrmtDecdQ.first;

			dcomprssPrmtQ.enq(prmt);

			cntlOrdrOn <= False;
		end

		dcomprssOrdrQ.enq(ordr);
		
		dcomprssOn <= True;
	endrule

	FIFO#(Tuple3#(Bit#(512), Bit#(32), Bit#(8))) make512bQ <- mkSizedBRAMFIFO(1024);
	FIFO#(Tuple4#(Bit#(512), Bit#(512), Bit#(32), Bit#(8))) make512bPipeQ <- mkSizedBRAMFIFO(1024);
	rule dcomprssPipeStage_1( dcomprssOn );
		dcomprssOrdrQ.deq;
		let ordr = dcomprssOrdrQ.first;
	
		if ( ordr == 0 ) begin 	// [CONTINUOUS UNMATCH]
			vbtmQ.deq;
			let v = vbtmQ.first;

			make512bQ.enq(tuple3(zeroExtend(v), 32, 1));
			pcktBodyTrck <= pcktBodyTrck + 1;
		end else begin 	       	// [MATCH THEN UNMATCH OR ANOTHER MATCH]
			refrQ.deq;
			dcomprssPrmtQ.deq;
			let r = refrQ.first;
			let p = dcomprssPrmtQ.first;
			let bytes 	= tpl_1(p);
			let pointer 	= tpl_2(p);
			let continuous 	= tpl_3(p);
			let starter	= tpl_4(p);
			let direction	= tpl_5(p);

			if ( bytes > 64 ) begin
				if ( starter && (pointer > 0) ) begin
					r = r >> pointer;
					continuous = continuous - 3;
				end

				make512bQ.enq(tuple3(r, 512, direction));
				dcomprssOrdrQ.enq(ordr);
				dcomprssPrmtQ.enq(tuple5(bytes - 64, pointer, continuous, False, direction));
				
				pcktBodyTrck <= pcktBodyTrck + 1;
			end else begin
				Bit#(32) length = 0;
				if ( starter ) begin
					length = continuous * 128;
					r = r >> pointer;
				end else begin
					length = truncate(pointer);
				end

				Bit#(32) remove = 512 - length;

				dcomprssPipeQ.enq(tuple4(r, remove, length, direction));
				
				if ( ordr == 1 ) begin
					dcomprssOrdrQ.enq(0);
					dcomprssOn <= False;
				end
				cntlOrdrOn <= True;
			end
		end
	endrule

	rule dcomprssPipeStage_2;
		dcomprssPipeQ.deq;
		let p = dcomprssPipeQ.first;
		let r		= tpl_1(p);
		let remove 	= tpl_2(p);
		let length 	= tpl_3(p);
		let direction 	= tpl_4(p);

		Bit#(512) rFinal = (r << remove) >> remove;

		dcomprssOn <= True;

		make512bQ.enq(tuple3(rFinal, length, direction));
	endrule

	Reg#(Bit#(512)) rsltBuf     <- mkReg(0);
	Reg#(Bit#(32))  rsltLngtBuf <- mkReg(0);
	Reg#(Bit#(64))  rsltAddrBuf <- mkReg(268435456);
	rule make512bPipeStage_1;
		make512bQ.deq;
		let p = make512bQ.first;
		let r 	      = tpl_1(p);
		let length    = tpl_2(p);
		let direction = tpl_3(p);

		Bit#(512) rsltTmpl = 0;
		if ( direction == 2 ) begin
			rsltTmpl = ~(reverseBits(r));
		end else begin
			rsltTmpl = r;
		end

		Bit#(512) rslt = rsltBuf | (rsltTmpl << rsltLngtBuf);
		
		make512bPipeQ.enq(tuple4(rsltTmpl, rslt, length, direction));
	endrule
	
	rule make512bPipeStage_2;
		make512bPipeQ.deq;
		let p = make512bPipeQ.first;
		let r 	      = tpl_1(p);
		let rslt      = tpl_2(p);
		let length    = tpl_3(p);
		let direction = tpl_4(p);

		if ( rsltLngtBuf + length >= 512 ) begin
			Bit#(32) icld = 512 - rsltLngtBuf;
			rqstWrRsltQ.enq(MemPortReq{addr:rsltAddrBuf, bytes:64});
			rsltQ.enq(rslt);

			rsltBuf     <= r >> icld;
			rsltLngtBuf <= length - icld;
			rsltAddrBuf <= rsltAddrBuf + 64;	
		end else begin
			rsltBuf     <= rslt;
			rsltLngtBuf <= rsltLngtBuf + length;
		end
	endrule
	//------------------------------------------------------------------------------------
	// Interface
	//------------------------------------------------------------------------------------
	// Read the compressed genomic data
	method Action readPckt(Bit#(512) pckt);
		serializer512b32b.put(pckt);
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
	method ActionValue#(Bit#(512)) wrteRslt;
		rsltQ.deq;
		return rsltQ.first;
	endmethod
endmodule

