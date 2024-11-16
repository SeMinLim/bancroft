import FIFO::*;
import FIFOF::*;
import Vector::*;

import Decompressor::*;


typedef 64 	 Kmer;
typedef 29111534 PcktCntTotal32b;
typedef 1819471  PcktCntTotal512b;
typedef 1 RsltCntTotal128b;
typedef 1 RsltCntTotal512b;
typedef 2 	 MemPortCnt;
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

	Reg#(Bool) started	<- mkReg(False);
	Reg#(Bool) rqstRdPcktOn <- mkReg(False);
	Reg#(Bool) readPcktOn 	<- mkReg(False);
	Reg#(Bool) rqstRdRefrOn <- mkReg(False);
	Reg#(Bool) readRefrOn 	<- mkReg(False);
	Reg#(Bool) rqstWrRsltOn <- mkReg(False);
	Reg#(Bool) wrteRsltOn 	<- mkReg(False);

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
		started      <= True;
		rqstRdPcktOn <= True;
		rqstRdRefrOn <= True;
		rqstWrRsltOn <= True;
	endrule
	//------------------------------------------------------------------------------------
	// [Memory Read]
	// 29,111,534 x 32-bit packet = 1,819,471 x 512-bit packet
	// 2,864,785,220-bp x 2-bit = 5,729,570,440-bit = 716,196,305-byte
	//------------------------------------------------------------------------------------
	Vector#(MemPortCnt, FIFO#(MemPortReq)) readReqQs <- replicateM(mkFIFO);
	Vector#(MemPortCnt, FIFO#(MemPortReq)) writeReqQs <- replicateM(mkFIFO);
	Vector#(MemPortCnt, FIFO#(Bit#(512))) writeWordQs <- replicateM(mkFIFO);
	Vector#(MemPortCnt, FIFO#(Bit#(512))) readWordQs <- replicateM(mkFIFO);

	// Read the compressed genomic data 	   [MEMPORT 0]
	Reg#(Bit#(32)) rqstRdPcktCnt <- mkReg(0);
	Reg#(Bit#(64)) addrBuf <- mkReg(0);
	rule rqstRdPckt( rqstRdPcktOn );
		readReqQs[0].enq(MemPortReq{addr:addrBuf, bytes:64});
	
		if ( rqstPcktCnt + 1 == fromInteger(valueOf(PcktCntTotal512b)) ) begin
			addrBuf       <= 0;
			rqstRdPcktCnt <= 0;
			rqstRdPcktOn  <= False;
		end else begin
			addrBuf       <= addrBuf + 64;
			rqstRdPcktCnt <= rqstRdPcktCnt + 1;
		end

		readPcktOn <= True;
	endrule
	Reg#(Bit#(32)) readPcktCnt <- mkReg(0);
	rule readPckt( readPcktOn );
		readWordQs[0].deq;
		let pckt = readWordQs[0].first;
	
		decompressor.readPckt(pckt);
		
		if ( readPcktCnt + 1 == fromInteger(valueOf(PcktCntTotal512b)) ) begin
			readPcktCnt <= 0;
			readPcktOn  <= False;
		end else begin
			readPcktCnt <= readPcktCnt + 1;
		end
	endrule

	// Read the 2-bit encoded reference genome [MEMPORT 1]
	rule rqstRdRefr( rqstRdRefrOn );
		let rqst <- decompressor.rqstRdRefr;
		readReqQs[1].enq(MemPortReq{addr:rqst.addr, bytes:rqst.bytes});

		readRefrOn <= True;
	endrule
	rule readRefr( readRefrOn );
		readWordQs[1].deq;
		let refr = readWordQs[1].first;

		decompressor.readRefr(refr);
	endrule
	//------------------------------------------------------------------------------------
	// [Memory Write] & [System Finish]
	// Memory Writer is going to use HBM[1] 
	// to store the decompressed 2-bit encoded sequence
	// 268,435,456 ~ 536,870,912 
	//------------------------------------------------------------------------------------
	Reg#(Bit#(32)) rqstWrRsltCnt <- mkReg(0);
	rule rqstWrRslt( rqstWrRsltOn );
		let rslt <- decompressor.rqstWrRslt;
		writeReqQs[0].enq(MemPortReq{addr:rslt.addr, bytes:rslt.bytes});

		if ( rqstWrRsltCnt + 1 == fromInteger(valueOf(RsltCntTotal512b)) ) begin
			rqstWrRsltCnt <= 0;
			rqstWrRsltOn  <= False;
		end else begin
			rqstWrRsltCnt <= rqstWrRsltCnt + 1;
		end

		wrteRsltOn <= True;
	endrule
	Reg#(Bit#(32)) wrteRsltCnt <- mkReg(0);
	rule wrteRslt( wrteRsltOn );
		let rslt <- decompressor.wrteRslt;
		writeWordQs[0].enq(rslt);

		// System Finish
		if ( wrteRsltCnt + 1 == fromInteger(valueOf(RsltCntTotal512b)) ) begin
			wrteRsltCnt <= 0;
			wrteRsltOn  <= False;
			started     <= False;
			doneQ.enq(True);
		end else begin
			wrteRsltCnt <= wrteRsltCnt + 1;
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
		startQ.enq(param);
	endmethod
	method ActionValue#(Bool) done;
		doneQ.deq;
		return doneQ.first;
	endmethod
	interface mem = mem_;
endmodule

