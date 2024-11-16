import FIFO::*;
import FIFOF::*;
import Vector::*;

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
	merhod ActionValue#(Bit#(512)) rslt;
endinterface
(* synthesize *)
module mkDecompressor( DecompressorIfc );
	SerializerIfc#(512, 16) serializer512b32b <- mkSerializer;

	FIFO#(Bit#(32)) pcktQ <- mkSizedFIFO(1024);
	FIFO#(Bit#(32)) vbtmQ <- mkSizedFIFO(1024);

	Reg#(Bit#(32)) pcktHeadTrck <- mkReg(0);
	Reg#(Bit#(32)) pcktBodyTrck <- mkReg(0);
	
	Reg#(Bool) readHeadOn <- mkReg(True);
	Reg#(Bool) intpHeadOn <- mkReg(False);
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

	FIFO#(Tuple3#(Bit#(32), Bit#(64), Bit#(32))) needPrmtDecdQ <- mkFIFO;
	FIFO#(Bit#(8))	ordrQ	      	<- mkFIFO;
	Reg#(Bit#(8))   idxDrctBuf    	<- mkReg(0);
	Reg#(Bit#(32))  idxStrtBuf    	<- mkReg(0);
	Reg#(Bit#(32))  idxCntnBuf    	<- mkReg(0);
	Reg#(Bit#(32))  intpHeadCurrCnt <- mkReg(0);
	Reg#(Bit#(32))  intpHeadTotlCnt <- mkReg(0);
	rule intpHead( intpHeadOn );
		headQ.deq;
		let head = headQ.first;
		
		case ( head )
			2'b00: // [UNMATCH]
			if ( idxDrctBuf != 0 ) begin
				let p = caclPrmtRqst(idxStrtBuf, idxCntnBuf);	

				rqstRdRefrQ.enq(MemPortReq{addr:tpl_1(p), bytes:tpl_2(p)});
				needPrmtDecdQ.enq(tuple3(tpl_2(p), tpl_3(p), idxCntnBuf));
				ordrQ.enq(1);
			
				idxDrctBuf <= 0;
				idxStrtBuf <= 0;
				idxCntnBuf <= 0;
			end else begin
				ordrQ.enq(0);
			end

			pcktQ.deq;
			vbtmQ.enq(pcktQ.first);
			2'b01: // [FORWARD NORMAL MATCH]
			if ( idxDrctBuf != 0 ) begin
				let p = caclPrmtRqst(idxStrtBuf, idxCntnBuf);

				rqstRdRefrQ.enq(MemPortReq{addr:tpl_1(p), bytes:tpl_2(p)});
				needPrmtDecdQ.enq(tuple3(tpl_2(p), tpl_3(p), idxCntnBuf));
				ordrQ.enq(2);

				idxCntnBuf <= 0;
			end

			pcktQ.deq;
			idxStrtBuf <= pcktQ.first;
			idxDrctBuf <= 1;
			2'b10: // [CONTINUOUS MATCH]
			if ( idxDrctBuf == 1 ) begin		// [FORWARD]
				idxCntnBuf <= idxCntnBuf + 1;
			end else if ( idxDirection == 2 ) begin	// [REVERSE]
				idxStrtBuf <= idxStrtBuf - fromInteger(valueOf(Kmer));
				idxCntnBuf <= idxCntnBuf + 1;
			end
			2'b11: // [REVERSE NORMAL MATCH]
			if ( idxDrctBuf != 0 ) begin
				let p = caclPrmtRqst(idxStrtBuf, idxCntnBuf);

				rqstRdRefrQ.enq(MemPortReq{addr:tpl_1(p), bytes:tpl_2(p)});
				needPrmtDecdQ.enq(tuple3(tpl_2(p), tpl_3(p), idxCntnBuf));
				orderQ.enq(2);

				idxCntnBuf <= 0;
			end

			pcktQ.deq;
			idxStrtBuf <= pcktQ.first;
			idxDrctBuf <= 2;
		endcase

		if ( intpHeadCurrCnt + 1 == fromInteger(valueOf(Head)) ) begin
			if ( intpHeadTotlCnt + 1 == fromInteger(TotalTrial) ) begin
				intpHeadCurrCnt <= intpHeadCurrCnt;
			end else begin
				readHeadOn 	<= True;
				intpHeadOn 	<= false;
				intpHeadCurrCnt <= 0;
			end
			intpHeadTotlCnt <= intpHeadTotlCnt + 1;
		end else begin
			dcomprssOn 	<= True;
			intpHeadCurrCnt <= intpHeadCurrCnt + 1;
			intpHeadTotlCnt <= intpHeadTotlCnt + 1;
		end
	endrule
	//------------------------------------------------------------------------------------
	// [STAGE 2]
	// Decompress the sequence & Send the decompressed data out as 512-bit
	//------------------------------------------------------------------------------------
	Reg#(Bit#(512)) rsltBuf <- mkReg(0);
	Reg#(Bit#(32)) rsltLngtBuf <- mkReg(0);
	Reg#(Bit#(32)) ordrBuf <- mkReg(0);
	Reg#(Bit#(32)) dcomprssCurrCnt <- mkReg(0);
	rule dcomprss( dcomprssOn );
		Bit#(8) r = 0;
		if ( dcomprssCurrCnt == 0 ) begin
			ordrQ.deq;
			r = ordrQ.first;
			ordrBuf <= r;
		end else begin
			r = ordrBuf;
		end

		case ( r )
			0: // [CONTINUOUS UNMATCH]
			vbtmQ.deq;
			let v = vbtmQ.first;

			if ( rsltLngtBuf + 32 >= 512 ) begin
				Bit#(512) rslt = rsltBuf | zeroExtend(v) << rsltLngtBuf;
				rsltQ.enq(rslt);

				rsltBuf     <= zeroExtend(v >> (512 - rsltLngtBuf));
				rsltLngtBuf <= 32 - (512 - rsltLngtBuf);
			end else begin
				rsltBuf     <= rsltBuf | zeroExtend(v) << rsltLngtBuf;
				rsltLngtBuf <= rsltLngtBuf + 32;
			end

			dcomprssCurrCnt <= 0;
			1: // [MATCH THEN UNMATCH]
			needPrmtDecdQ.deq;

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
	method ActionValue#(Bit#(512)) rslt;
		rsltQ.deq;
		return rsltQ.first;
	endmethod
endmodule

