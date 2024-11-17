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
	method ActionValue#(Bit#(512)) wrteRslt;
endinterface
(* synthesize *)
module mkDecompressor( DecompressorIfc );
	SerializerIfc#(512, 16) serializer512b32b <- mkSerializer;

	FIFO#(MemPortReq) rqstRdRefrQ <- mkFIFO;
	FIFO#(MemPortReq) rqstWrRsltQ <- mkFIFO;
	FIFO#(Bit#(32))  pcktQ <- mkSizedFIFO(1024);
	FIFO#(Bit#(32))  vbtmQ <- mkSizedFIFO(1024);
	FIFO#(Bit#(512)) refrQ <- mkFIFO;
	FIFO#(Bit#(512)) rsltQ <- mkFIFO;

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
		
		if ( head == 2'b00 ) begin 	    // [UNMATCH]
			if ( idxDrctBuf != 0 ) begin
				let p = calcPrmtRqst(idxStrtBuf, idxCntnBuf);	

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
		end else if ( head == 2'b01 ) begin // [FORWARD NORMAL MATCH]
			if ( idxDrctBuf != 0 ) begin
				let p = calcPrmtRqst(idxStrtBuf, idxCntnBuf);

				rqstRdRefrQ.enq(MemPortReq{addr:tpl_1(p), bytes:tpl_2(p)});
				needPrmtDecdQ.enq(tuple3(tpl_2(p), tpl_3(p), idxCntnBuf));
				ordrQ.enq(2);

				idxCntnBuf <= 0;
			end

			pcktQ.deq;
			idxStrtBuf <= pcktQ.first;
			idxDrctBuf <= 1;
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
				needPrmtDecdQ.enq(tuple3(tpl_2(p), tpl_3(p), idxCntnBuf));
				ordrQ.enq(2);

				idxCntnBuf <= 0;
			end

			pcktQ.deq;
			idxStrtBuf <= pcktQ.first;
			idxDrctBuf <= 2;
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
			dcomprssOn 	<= True;
			intpHeadCurrCnt <= intpHeadCurrCnt + 1;
			intpHeadTotlCnt <= intpHeadTotlCnt + 1;
		end
	endrule
	//------------------------------------------------------------------------------------
	// [STAGE 2]
	// Decompress the sequence & Send the decompressed data out as 512-bit
	//------------------------------------------------------------------------------------
	Reg#(Tuple3#(Bit#(32), Bit#(64), Bit#(32))) needPrmtDecdBuf <- mkReg(?);

	Reg#(Bit#(512)) rsltBuf <- mkReg(0);
	Reg#(Bit#(32))  rsltLngtBuf <- mkReg(0);
	Reg#(Bit#(64))  rsltAddrBuf <- mkReg(268435456);

	Reg#(Bit#(8))  ordrBuf <- mkReg(0);
	Reg#(Bool)     ordrDqueOn <- mkReg(True);

	Reg#(Bit#(32)) dcomprssCurrCnt <- mkReg(0);
	Reg#(Bit#(32)) dcomprssGoalCnt <- mkReg(0);
	rule dcomprss( dcomprssOn );
		Bit#(8) ordr = 0;
		if ( ordrDqueOn ) begin
			ordrQ.deq;
			ordr = ordrQ.first;
		end else begin
			ordr = ordrBuf;
		end

		if ( ordr == 0 ) begin // [CONTINUOUS UNMATCH]
			vbtmQ.deq;
			let v = vbtmQ.first;

			if ( rsltLngtBuf + 32 >= 512 ) begin
				Bit#(512) rslt = rsltBuf | zeroExtend(v) << rsltLngtBuf;
				rqstWrRsltQ.enq(MemPortReq{addr:rsltAddrBuf, bytes:64});
				rsltQ.enq(rslt);

				rsltBuf     <= zeroExtend(v >> (512 - rsltLngtBuf));
				rsltLngtBuf <= 32 - (512 - rsltLngtBuf);
				rsltAddrBuf <= rsltAddrBuf + 64;
			end else begin
				rsltBuf     <= rsltBuf | zeroExtend(v) << rsltLngtBuf;
				rsltLngtBuf <= rsltLngtBuf + 32;
			end

			ordrDqueOn   <= True;
		end else begin 		// [MATCH THEN UNMATCH OR ANOTHER MATCH]
			refrQ.deq;
			let r = refrQ.first;

			Tuple3#(Bit#(32), Bit#(64), Bit#(32)) p = tuple3(0, 0, 0);
			if ( dcomprssCurrCnt == 0 ) begin
				needPrmtDecdQ.deq;
				p = needPrmtDecdQ.first;
				r = r >> tpl_2(p);
			end else begin
				p = needPrmtDecdBuf;
			end
				
			if ( tpl_1(p) > 64 ) begin
				Bit#(512) rslt = rsltBuf | r << rsltLngtBuf;
				rqstWrRsltQ.enq(MemPortReq{addr:rsltAddrBuf, bytes:64});	
				rsltQ.enq(rslt);

				rsltBuf     <= r >> (512 - rsltLngtBuf);
				rsltAddrBuf <= rsltAddrBuf + 64;

				Bit#(32) idxCntn = 0;
				if ( tpl_2(p) > 0 ) idxCntn = tpl_3(p) - 3;
				else idxCntn = tpl_3(p) - 4;

				needPrmtDecdBuf <= tuple3(tpl_1(p) - 64, tpl_2(p), idxCntn);
				dcomprssCurrCnt <= dcomprssCurrCnt + 1;
				
				ordrBuf    <= ordr;
				ordrDqueOn <= False;
			end else if ( tpl_1(p) == 64 )begin
				Bit#(32) currLngt = 0;
				if ( dcomprssCurrCnt == 0 ) begin
					currLngt = tpl_3(p) * 128;
					r = r << (512 - currLngt);
					r = r >> (512 - currLngt);
				end else begin
					currLngt = truncate(tpl_2(p));
					r = r << (512 - currLngt);
					r = r >> (512 - currLngt);
				end

				if ( rsltLngtBuf + currLngt > 512 ) begin
					Bit#(512) rslt = rsltBuf | r << rsltLngtBuf;
					rqstWrRsltQ.enq(MemPortReq{addr:rsltAddrBuf, bytes:64});
					rsltQ.enq(rslt);

					rsltBuf     <= r >> (512 - rsltLngtBuf);
					rsltLngtBuf <= currLngt - (512 - rsltLngtBuf);
					rsltAddrBuf <= rsltAddrBuf + 64;
				end else begin
					rsltBuf     <= rsltBuf | r << rsltLngtBuf;
					rsltLngtBuf <= rsltLngtBuf + currLngt;
				end

				needPrmtDecdBuf <= tuple3(0, 0, 0);
				dcomprssCurrCnt <= 0;
				
				if ( ordr == 1 )      ordrBuf 	 <= 0;
				else if ( ordr == 2 ) ordrDqueOn <= True;
			end
		end

		if ( pcktBodyTrck + 1 == fromInteger(valueOf(PcktCntBody32b)) ) begin
			dcomprssOn   <= False;
			pcktBodyTrck <= 0;
		end else begin
			pcktBodyTrck <= pcktBodyTrck  + 1;
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

