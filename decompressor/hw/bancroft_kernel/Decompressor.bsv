import FIFO::*;
import FIFOF::*;
import Vector::*;

import BRAM::*;
import BRAMFIFO::*;

import Serializer::*;
import BLShifter::*;


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


function Tuple3#(Bit#(64), Bit#(32), Bit#(32)) calcPrmtRqst( Bit#(32) idxStrt, Bit#(32) idxCntn );
	Bit#(64) address = (zeroExtend(idxStrt) * 2) / 8;
	Bit#(32) idxPntr = (idxStrt * 2) % 8;

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
	BLShiftIfc#(Bit#(512), 32, 2) shifterLeftDecompressPipe_1  <- mkPipelinedShift(False); // Shift Left
	BLShiftIfc#(Bit#(512), 32, 2) shifterRightDecompressPipe_1 <- mkPipelinedShift(True);  // Shirt Right
	BLShiftIfc#(Bit#(512), 32, 2) shifterLeftDecompressPipe_2  <- mkPipelinedShift(False); // Shift Left
	BLShiftIfc#(Bit#(512), 32, 2) shifterRightDecompressPipe_2 <- mkPipelinedShift(True);  // Shirt Right
	BLShiftIfc#(Bit#(512), 32, 2) shifterRightDecompressPipe_3 <- mkPipelinedShift(True);  // Shirt Right
	BLShiftIfc#(Bit#(512), 32, 2) shifterLeft512bComposerPipe_1  <- mkPipelinedShift(False); // Shift Left
	BLShiftIfc#(Bit#(512), 32, 2) shifterRight512bComposerPipe_2 <- mkPipelinedShift(True);  // Shirt Right

	FIFO#(MemPortReq) rqstRdRefrQ <- mkSizedBRAMFIFO(1024);
	FIFO#(MemPortReq) rqstWrRsltQ <- mkSizedBRAMFIFO(1024);
	FIFO#(Bit#(32))  pcktQ <- mkSizedBRAMFIFO(1024);
	FIFO#(Bit#(32))  vbtmQ <- mkSizedBRAMFIFO(1024);
	FIFO#(Bit#(512)) refrQ <- mkSizedBRAMFIFO(1024);
	FIFO#(Bit#(512)) rsltQ <- mkSizedBRAMFIFO(1024);

	Reg#(Bit#(32)) pcktHeadTrck <- mkReg(0);
	
	Reg#(Bool) readHeadOn	<- mkReg(True);
	Reg#(Bool) intpHeadOn	<- mkReg(False);
	Reg#(Bool) cntlOrdrOn	<- mkReg(False);
	Reg#(Bool) dcomprssOn	<- mkReg(False);
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

	// Bytes, Pointer, Continuous, Direction, First, ShiftOn
	FIFO#(Tuple5#(Bit#(32), Bit#(32), Bit#(32), Bit#(8), Bool)) needPrmtDecdQ <- mkSizedBRAMFIFO(1024);
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
				let address = tpl_1(p);	
				let bytes   = tpl_2(p);
				let pointer = tpl_3(p);

				rqstRdRefrQ.enq(MemPortReq{addr:address, bytes:bytes});
				needPrmtDecdQ.enq(tuple5(bytes, pointer, idxCntnBuf, idxDrctBuf, True));
			
				idxDrctBuf <= 0;
				idxStrtBuf <= 0;
				idxCntnBuf <= 0;
			end

			pcktQ.deq;
			vbtmQ.enq(pcktQ.first);
			ordrQ.enq(0);
		end else if ( head == 2'b01 ) begin // [FORWARD NORMAL MATCH]
			if ( idxDrctBuf != 0 ) begin
				let p = calcPrmtRqst(idxStrtBuf, idxCntnBuf);
				let address = tpl_1(p);
				let bytes   = tpl_2(p);
				let pointer = tpl_3(p);

				rqstRdRefrQ.enq(MemPortReq{addr:address, bytes:bytes});
				needPrmtDecdQ.enq(tuple5(bytes, pointer, idxCntnBuf, idxDrctBuf, True));
			end

			pcktQ.deq;
			idxStrtBuf <= pcktQ.first;
			idxDrctBuf <= 1;
			idxCntnBuf <= 1;
			ordrQ.enq(1);
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
				let address = tpl_1(p);
				let bytes   = tpl_2(p);
				let pointer = tpl_3(p);

				rqstRdRefrQ.enq(MemPortReq{addr:address, bytes:bytes});
				needPrmtDecdQ.enq(tuple5(bytes, pointer, idxCntnBuf, idxDrctBuf, True));
			end

			pcktQ.deq;
			idxStrtBuf <= pcktQ.first;
			idxDrctBuf <= 2;
			idxCntnBuf <= 1;
			ordrQ.enq(2);
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
			cntlOrdrOn 	<= True;
			intpHeadCurrCnt <= intpHeadCurrCnt + 1;
			intpHeadTotlCnt <= intpHeadTotlCnt + 1;
		end
	endrule
	//------------------------------------------------------------------------------------
	// [STAGE 2]
	// Decompress the sequence
	//------------------------------------------------------------------------------------
	// Pipeline Stage 1
	// Type 1 => Right Shifter
	// Type 2 => Left  Shifter
	//------------------------------------------------------------------------------------
//	FIFO#(Tuple4#(Bit#(8), Bit#(32), Bit#(32), Bit#(8))) dcomprssPipe_1_Type_1_Q <- mkSizedBRAMFIFO(512);
//	FIFO#(Tuple4#(Bit#(8), Bit#(32), Bit#(32), Bit#(8))) dcomprssPipe_1_Type_2_Q <- mkSizedBRAMFIFO(512);
//	FIFO#(Tuple4#(Bit#(8), Bit#(32), Bit#(32), Bit#(8))) dcomprssPipe_2_Type_1_Q <- mkSizedBRAMFIFO(512);
//	FIFO#(Tuple4#(Bit#(8), Bit#(32), Bit#(32), Bit#(8))) dcomprssPipe_2_Type_2_Q <- mkSizedBRAMFIFO(512);
//	FIFO#(Tuple4#(Bit#(8), Bit#(32), Bit#(32), Bit#(8))) dcomprssPipe_3_Type_1_Q <- mkSizedBRAMFIFO(512);
	FIFO#(Tuple5#(Bit#(32), Bit#(32), Bit#(32), Bit#(8), Bool)) dcomprssPrmtQ    <- mkSizedBRAMFIFO(512);
	FIFO#(Bit#(8)) dcomprssOrdrQ <- mkSizedBRAMFIFO(512);
	Reg#(Bit#(64))  rsltAddrBuf <- mkReg(268435456);

//	FIFO#(Tuple3#(Bit#(512), Bit#(32), Bit#(8))) make512bQ <- mkSizedBRAMFIFO(512);
//	FIFO#(Tuple3#(Bit#(512), Bit#(32), Bit#(8))) make512b_Pipe_1_Q <- mkSizedBRAMFIFO(512);
//	FIFO#(Tuple3#(Bit#(512), Bit#(32), Bit#(8))) make512b_Pipe_2_Q <- mkSizedBRAMFIFO(512);
	rule dcomprssMain( dcomprssOn );
		Bit#(8) ordr = 0;
		if ( cntlOrdrOn ) begin
			ordrQ.deq;
			ordr = ordrQ.first;
		end else begin
			dcomprssOrdrQ.deq;
			ordr = dcomprssOrdrQ.first;
		end
	
		if ( ordr == 0 ) begin 	// [UNMATCH]
			vbtmQ.deq;
			let v = vbtmQ.first;

			//make512bQ.enq(tuple3(zeroExtend(v), 32, 1));
			rqstWrRsltQ.enq(MemPortReq{addr:rsltAddrBuf, bytes:64});
			rsltQ.enq(zeroExtend(v));
		end else begin 	       	// [MATCH]
			Tuple5#(Bit#(32), Bit#(32), Bit#(32), Bit#(8), Bool) param;
			if ( cntlOrdrOn ) begin
				needPrmtDecdQ.deq;
				param = needPrmtDecdQ.first;
			end else begin
				dcomprssPrmtQ.deq;
				param = dcomprssPrmtQ.first;
			end
			let bytes 	= tpl_1(param);
			let pointer 	= tpl_2(param);
			let continuous 	= tpl_3(param);
			let direction 	= tpl_4(param);
			let first	= tpl_5(param);
		
			refrQ.deq;
			let value = refrQ.first;
			if ( bytes > 64 ) begin
				if ( first && (pointer > 0) ) begin
					shifterRightDecompressPipe_1.enq(value, pointer);
					dcomprssOrdrQ.enq(ordr);
					dcomprssPrmtQ.enq(tuple5(bytes - 64, pointer, continuous, direction, False));
					//dcomprssPipe_1_Type_1_Q.enq(tuple4(0, pointer, continuous, direction));
				end else begin
					//make512bQ.enq(tuple3(value, 512, direction));
					dcomprssOrdrQ.enq(ordr);
					dcomprssPrmtQ.enq(tuple5(bytes - 64, pointer, continuous, direction, False));
					rqstWrRsltQ.enq(MemPortReq{addr:rsltAddrBuf, bytes:64});
					rsltQ.enq(value);
				end
				cntlOrdrOn <= False;
			end else begin
				if ( first ) begin
					shifterRightDecompressPipe_1.enq(value, pointer);
					//dcomprssPipe_1_Type_1_Q.enq(tuple4(1, pointer, continuous, direction));
				end else begin
					shifterLeftDecompressPipe_1.enq(value, (512 - pointer));
					//dcomprssPipe_1_Type_2_Q.enq(tuple4(0, pointer, continuous, direction));
				end
				cntlOrdrOn <= True;
			end
		end
	endrule
	rule shifterRight;
		shifterRightDecompressPipe_1.deq;
		let value = shifterRightDecompressPipe_1.first;
		rqstWrRsltQ.enq(MemPortReq{addr:rsltAddrBuf, bytes:64});
		rsltQ.enq(value);
	endrule
	rule shifterLeft;
		shifterLeftDecompressPipe_1.deq;
		let value = shifterLeftDecompressPipe_1.first;
		rqstWrRsltQ.enq(MemPortReq{addr:rsltAddrBuf, bytes:64});
		rsltQ.enq(value);
	endrule
	//------------------------------------------------------------------------------------
	// Pipeline Stage 2
	//------------------------------------------------------------------------------------
/*	rule dcomprssPipe_1_Type_1;
		shifterRightDecompressPipe_1.deq;
		dcomprssPipe_1_Type_1_Q.deq;
		let value = shifterRightDecompressPipe_1.first;
		let param = dcomprssPipe_1_Type_1_Q.first;
		let cases 	= tpl_1(param);
		let pointer 	= tpl_2(param);
		let continuous 	= tpl_3(param);
		let direction 	= tpl_4(param);

		if ( cases == 0 ) begin
			make512bQ.enq(tuple3(value, 512, direction));
		end else begin
			shifterLeftDecompressPipe_2.enq(value, (512 - (continuous * 128)));
			dcomprssPipe_2_Type_2_Q.enq(tuple4(cases, pointer, continuous, direction));
		end
	endrule
	rule dcomprssPipe_1_Type_2;
		shifterLeftDecompressPipe_1.deq;
		dcomprssPipe_1_Type_2_Q.deq;
		let value = shifterLeftDecompressPipe_1.first;
		let param = dcomprssPipe_1_Type_2_Q.first;
		let cases 	= tpl_1(param);
		let pointer	= tpl_2(param);
		let continuous	= tpl_3(param);
		let direction	= tpl_4(param);

		shifterRightDecompressPipe_2.enq(value, (512 - pointer));
		dcomprssPipe_2_Type_1_Q.enq(tuple4(cases, pointer, continuous, direction));
	endrule
	//------------------------------------------------------------------------------------
	// Pipeline Stage 3
	//------------------------------------------------------------------------------------
	rule dcomprssPipe_2_Type_1;
		shifterRightDecompressPipe_2.deq;
		dcomprssPipe_2_Type_1_Q.deq;
		let value = shifterRightDecompressPipe_2.first;
		let param = dcomprssPipe_2_Type_1_Q.first;
		let cases 	= tpl_1(param);
		let pointer 	= tpl_2(param);
		let continuous 	= tpl_3(param);
		let direction 	= tpl_4(param);
		
		make512bQ.enq(tuple3(value, pointer, direction));
	endrule
	rule dcomprssPipe_2_Type_2;
		shifterLeftDecompressPipe_2.deq;
		dcomprssPipe_2_Type_2_Q.deq;
		let value = shifterLeftDecompressPipe_2.first;
		let param = dcomprssPipe_2_Type_2_Q.first;
		let cases 	= tpl_1(param);
		let pointer	= tpl_2(param);
		let continuous	= tpl_3(param);
		let direction	= tpl_4(param);

		shifterRightDecompressPipe_3.enq(value, (512 - (continuous * 128)));
		dcomprssPipe_3_Type_1_Q.enq(tuple4(cases, pointer, continuous, direction));
	endrule
	//------------------------------------------------------------------------------------
	// Pipeline Stage 4
	//------------------------------------------------------------------------------------
	rule dcomprssPipe_3_Type_1;
		shifterRightDecompressPipe_3.deq;
		dcomprssPipe_3_Type_1_Q.deq;
		let value = shifterRightDecompressPipe_3.first;
		let param = dcomprssPipe_3_Type_1_Q.first;
		let cases 	= tpl_1(param);
		let pointer 	= tpl_2(param);
		let continuous 	= tpl_3(param);
		let direction 	= tpl_4(param);
		
		make512bQ.enq(tuple3(value, continuous * 128, direction));
	endrule
	//------------------------------------------------------------------------------------
	// [STAGE 3]
	// Compose & Send the decompressed data out as 512-bit
	//------------------------------------------------------------------------------------
	Reg#(Bit#(512)) rsltBuf     <- mkReg(0);
	Reg#(Bit#(32))  rsltLngtBuf <- mkReg(0);
	Reg#(Bit#(64))  rsltAddrBuf <- mkReg(268435456);
	rule make512bPipeStage_1;
		make512bQ.deq;
		let param	= make512bQ.first;
		let value	= tpl_1(param);
		let length    	= tpl_2(param);
		let direction 	= tpl_3(param);

		Bit#(512) resultTmp = 0;
		if ( direction == 2 ) begin
			resultTmp = ~(reverseBits(value));
		end else begin
			resultTmp = value;
		end
		
		shifterLeft512bComposerPipe_1.enq(resultTmp, rsltLngtBuf);
		make512b_Pipe_1_Q.enq(tuple3(value, length, direction));
	endrule

	rule make512bPipeStage_2;
		shifterLeft512bComposerPipe_1.deq;
		make512b_Pipe_1_Q.deq;
		let resultTmp 	= shifterLeft512bComposerPipe_1.first;
		let param	= make512b_Pipe_1_Q.first;
		let value	= tpl_1(param);
		let length    	= tpl_2(param);
		let direction 	= tpl_3(param);
		
		Bit#(512) results = resultTmp | rsltBuf;

		if ( rsltLngtBuf + length >= 512 ) begin
			rqstWrRsltQ.enq(MemPortReq{addr:rsltAddrBuf, bytes:64});
			rsltQ.enq(results);
			shifterRight512bComposerPipe_2.enq(value, (512 - rsltLngtBuf));
			make512b_Pipe_2_Q.enq(tuple3(value, length, direction));
		end else begin
			rsltBuf     <= results;
			rsltLngtBuf <= rsltLngtBuf + length;
		end
	endrule
	
	rule make512bPipeStage_3;
		shifterRight512bComposerPipe_2.deq;
		make512b_Pipe_2_Q.deq;
		let remain	= shifterRight512bComposerPipe_2.first;
		let param	= make512b_Pipe_2_Q.first;
		let value	= tpl_1(param);
		let length    	= tpl_2(param);
		let direction 	= tpl_3(param);

		rsltBuf     <= remain;
		rsltLngtBuf <= length - (512 - rsltLngtBuf);
		rsltAddrBuf <= rsltAddrBuf + 64;	
	endrule*/
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

