import FIFO::*;
import FIFOF::*;
import Vector::*;

import Decompressor::*;


typedef 64 Kmer;

typedef 343525 	DataCntHead32b;
typedef 2422857 DataCntBody32b;
typedef TAdd#(DataCntHead32b, DataCntBody32b) DataCntTotal32b;
typedef TDiv#(DataCntTotal32b, 16) DataCntTotal512b;

typedef 1747088 ResultCntTotalStride;
typedef 3749312 ResultCntTotalKmer;
typedef TMul#(ResultCntTotalStride, 32) ResultCntTotalStrideLength;
typedef TMul#(ResultCntTotalKmer, 128) ResultCntTotalKmerLength;
typedef TAdd#(ResultCntTotalStrideLength, ResultCntTotalKmerLength) ResultCntTotalLength;
typedef TDiv#(ResultCntTotalLength, 512) ResultCntTotal;

typedef 2 MemPortCnt;
typedef struct {
	Bit#(64) addr;
	Bit#(32) bytes;
} MemPortReq deriving (Eq,Bits);


interface MemPortIfc;
	method ActionValue#(MemPortReq) readReq;
	method ActionValue#(MemPortReq) writeReq;
	method ActionValue#(Bit#(512)) writeWord;
	method Action readWord(Bit#(512) word);
endinterface


interface KernelMainIfc;
	method Action start(Bit#(32) param);
	method ActionValue#(Bool) done;
	interface Vector#(MemPortCnt, MemPortIfc) mem;
endinterface
module mkKernelMain(KernelMainIfc);
	FIFO#(Bool) startQ <- mkFIFO;
	FIFO#(Bool) doneQ  <- mkFIFO;

	Reg#(Bool) started 		<- mkReg(False);
	Reg#(Bool) reqReadDataOn 	<- mkReg(False);
	Reg#(Bool) readDataOn 		<- mkReg(False);
	Reg#(Bool) reqReadRefOn 	<- mkReg(False);
	Reg#(Bool) readRefOn 		<- mkReg(False);
	Reg#(Bool) reqWriteResultOn 	<- mkReg(False);
	Reg#(Bool) writeResultOn 	<- mkReg(False);

	DecompressorIfc decompressor <- mkDecompressor;
	//------------------------------------------------------------------------------------
	// [Cycle Counter]
	//------------------------------------------------------------------------------------
	Reg#(Bit#(32)) cycleCounter <- mkReg(0);
	rule incCycle;
		cycleCounter <= cycleCounter + 1;
	endrule
	//------------------------------------------------------------------------------------
	// [System Start]
	//------------------------------------------------------------------------------------
	rule systemStart( !started );
		startQ.deq;
		started      		<= True;
		reqReadDataOn	 	<= True;
		reqReadRefOn 		<= True;
		reqWriteResultOn 	<= True;
	endrule
	//------------------------------------------------------------------------------------
	// [Memory Read]
	//------------------------------------------------------------------------------------
	Vector#(MemPortCnt, FIFO#(MemPortReq)) readReqQs <- replicateM(mkFIFO);
	Vector#(MemPortCnt, FIFO#(MemPortReq)) writeReqQs <- replicateM(mkFIFO);
	Vector#(MemPortCnt, FIFO#(Bit#(512))) writeWordQs <- replicateM(mkFIFO);
	Vector#(MemPortCnt, FIFO#(Bit#(512))) readWordQs <- replicateM(mkFIFO);

	// Read the compressed genomic data 	   [MEMPORT 0]
	Reg#(Bit#(32)) reqReadDataCnt <- mkReg(0);
	Reg#(Bit#(64)) addr <- mkReg(0);
	rule reqReadData( reqReadDataOn );
		readReqQs[0].enq(MemPortReq{addr:addr, bytes:64});
		if ( reqReadDataCnt + 1 == fromInteger(valueOf(DataCntTotal512b)) ) begin
			addr 		<= 0;
			reqReadDataCnt 	<= 0;
			reqReadDataOn 	<= False;
		end else begin
			addr 		<= addressR + 64;
			reqReadDataCnt 	<= requestReadDataCnt + 1;
		end
		readDataOn <= True;
	endrule
	Reg#(Bit#(32)) readDataCnt <- mkReg(0);
	rule readData( readDataOn );
		readWordQs[0].deq;
		let data = readWordQs[0].first;
	
		decompressor.readData(data);
		
		if ( readDataCnt + 1 == fromInteger(valueOf(DataCntTotal512b)) ) begin
			readDataCnt <= 0;
			readDataOn  <= False;
		end else begin
			readDataCnt <= readDataCnt + 1;
		end
	endrule

	// Read the 2-bit encoded reference genome [MEMPORT 1]
	Reg#(Bit#(32)) reqReadRefCnt <- mkReg(0);
	Reg#(Bit#(64)) reqReadRefAddr <- mkReg(0);
	Reg#(Bit#(32)) reqReadRefByte <- mkReg(0);
	rule reqReadRef( reqReadRefOn );
		if ( reqReadRefCnt == 0 ) begin
			let req <- decompressor.reqReadRef;
			readReqQs[1].enq(MemPortReq{addr:rqst.addr, bytes:64});
			
			if ( req.bytes - 64 >= 64 ) begin
				reqReadRefAddr 	<= req.addr  + 64;
				reqReadRefByte	<= req.bytes - 64;
				reqReadRefCnt	<= reqReadRefCnt + 1;
			end
		end else begin
			readReqQs[1].enq(MemPortReq{addr:reqReadRefAddr, bytes:64});

			if ( reqReadRefByte - 64 >= 64 ) begin
				reqReadRefAddr 	<= reqReadRefAddr + 64;
				reqReadRefByte 	<= reqReadRefByte - 64;
				reqReadRefCnt	<= reqReadRefCnt + 1;
			end else begin
				reqReadRefAddr 	<= 0;
				reqReadRefByte 	<= 0;
				reqReadRefCnt	<= 0;
			end
		end

		readRefOn <= True;
	endrule
	rule readRef( readRefOn );
		readWordQs[1].deq;
		let ref = readWordQs[1].first;

		decompressor.readRef(ref);
	endrule
	//------------------------------------------------------------------------------------
	// [Memory Write] & [System Finish]
	// Memory Writer is going to use HBM[1] 
	// to store the decompressed 2-bit encoded sequence
	// 268,435,456 ~ 536,870,912 
	// UPDATE!!! Cycle count number will be written at 268,435,456
	//------------------------------------------------------------------------------------
	Reg#(Bit#(32)) reqWriteResultCnt <- mkReg(0);
	rule reqWriteResult( reqWriteResultOn );
		let result <- decompressor.reqWriteResult;
		writeReqQs[0].enq(MemPortReq{addr:result.addr, bytes:result.bytes});

		if ( reqWriteResultCnt + 1 == fromInteger(valueOf(ResultCntTotal512b)) ) begin
			reqWriteResultCnt <= 0;
			reqWriteResultOn  <= False;
		end else begin
			reqWriteResultCnt <= reqWriteResultCnt + 1;
		end

		writeResultOn <= True;
	endrule
	Reg#(Bit#(32)) writeResultCnt <- mkReg(0);
	rule writeResult( writeResultOn );
		let result <- decompressor.writeResult;
		writeWordQs[0].enq(result);

		// System Finish
		if ( writeResultCnt + 1 == fromInteger(valueOf(ResultCntTotal512b)) ) begin
			writeResultCnt 	<= 0;
			writeResultOn  	<= False;
			started		<= False;
			doneQ.enq(True);
		end else begin
			writeResultCnt <= writeResultCnt + 1;
		end
	endrule
	//------------------------------------------------------------------------------------
	// Interface
	//------------------------------------------------------------------------------------
	Vector#(MemPortCnt, MemPortIfc) mem_;
	for (Integer i = 0; i < valueOf(MemPortCnt); i=i+1) begin
		mem_[i] = interface MemPortIfc;
			method ActionValue#(MemPortReq) readReq;
				readReqQs[i].deq;
				return readReqQs[i].first;
			endmethod
			method ActionValue#(MemPortReq) writeReq;
				writeReqQs[i].deq;
				return writeReqQs[i].first;
			endmethod
			method ActionValue#(Bit#(512)) writeWord;
				writeWordQs[i].deq;
				return writeWordQs[i].first;
			endmethod
			method Action readWord(Bit#(512) word);
				readWordQs[i].enq(word);
			endmethod
		endinterface;
	end
	method Action start(Bit#(32) param) if ( started == False );
		startQ.enq(True);
	endmethod
	method ActionValue#(Bool) done;
		doneQ.deq;
		return doneQ.first;
	endmethod
	interface mem = mem_;
endmodule
