import FIFO::*;
import FIFOF::*;
import Vector::*;

import BRAM::*;
import BRAMFIFO::*;

import Serializer::*;
import BLShifter::*;


typedef 16 Head;
typedef 64 Kmer;
typedef 343525 	DataCntHead32b;
typedef TMul#(DataCntHead32b, Head) DataCntHead2b;
typedef struct {
	Bit#(64) addr;
	Bit#(32) bytes;
} MemPortReq deriving (Eq,Bits);


function Tuple3#(Bit#(64), Bit#(32), Bit#(32)) getParameters( Bit#(32) start, Bit#(32) continuous );
	Bit#(64) address = (zeroExtend(start) * 2) / 8;
	Bit#(32) pointer = (start * 2) % 8;

	Bit#(32) bytes 	  = 0;
	Bit#(32) bytesTmp = continuous / 16;
	if ( continuous % 16  == 0 ) begin
		if ( pointer > 0 ) bytes = (bytesTmp + 1) * 64;
		else 		   bytes = bytesTmp * 64;
	end else begin
		bytes = (bytesTmp + 1) * 64;
	end
	
	return tuple3(address, bytes, pointer);
endfunction


interface DecompressorIfc;
	method Action readData(Bit#(512) data);
	method ActionValue#(MemPortReq) reqReadRef;
	method Action readRef(Bit#(512) reference);
	method ActionValue#(MemPortReq) reqWriteResult;
	method ActionValue#(Bit#(512)) writeResult;
endinterface
(* synthesize *)
module mkDecompressor( DecompressorIfc );
	SerializerIfc#(512, 16) serializer512b32b 	<- mkSerializer;
	SerializerIfc#(32, 16) serializer32b2b 		<- mkSerializer;

	FIFO#(MemPortReq) reqReadRefQ 		<- mkSizedBRAMFIFO(1024);
	FIFO#(MemPortReq) reqWriteResultQ 	<- mkSizedBRAMFIFO(1024);
	FIFO#(Bit#(32))  dataQ 			<- mkSizedBRAMFIFO(2048);
	FIFO#(Bit#(32))  verbatimQ 		<- mkSizedBRAMFIFO(1024);
	FIFO#(Bit#(512)) refQ			<- mkSizedBRAMFIFO(1024);
	FIFO#(Bit#(512)) resultQ 		<- mkSizedBRAMFIFO(2048);
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
	rule serialize512b32b;
		let data <- serializer512b32b.get;
		dataQ.enq(data);
	endrule
	//------------------------------------------------------------------------------------
	// [STAGE 1]
	// Interpret the head section first
	//------------------------------------------------------------------------------------
	FIFO#(Bit#(2)) headQ <- mkSizedBRAMFIFO(1024);
	Reg#(Bit#(32)) headBuf 	   <- mkReg(0);
	Reg#(Bit#(32)) readHeadCnt <- mkReg(0);
	Reg#(Bool) readHeadOn <- mkReg(True);
	rule readHead(readHeadOn);
		dataQ.deq;
		serializer32b2b.put(dataQ.first);
		readHeadOn <= False;
	endrule

	// Bytes, Pointer, Continuous, Direction, First
	FIFO#(Tuple5#(Bit#(32), Bit#(32), Bit#(32), Bit#(8), Bool)) parameterQ <- mkSizedBRAMFIFO(1024);
	FIFO#(Bit#(8)) caseQ	      	<- mkSizedBRAMFIFO(1024);
	Reg#(Bit#(8)) direction    	<- mkReg(0);
	Reg#(Bit#(32)) start    	<- mkReg(0);
	Reg#(Bit#(32)) continuous    	<- mkReg(0);
	Reg#(Bit#(32)) interpretHeadCnt <- mkReg(0);
	rule interpretHead;
		Bit#(2) head = 0;
		if ( interpretHeadCnt == fromInteger(valueOf(DataCntHead2b)) ) begin
			head = 2'b00;
		end else begin
			head <- serializer32b2b.get;
		end
		
		if ( head == 2'b00 ) begin 	    // [UNMATCH]
			if ( direction != 0 ) begin
				let p = getParameters(start, continuous);
				let address = tpl_1(p);	
				let bytes   = tpl_2(p);
				let pointer = tpl_3(p);

				reqReadRefQ.enq(MemPortReq{addr:address, bytes:bytes});
				parameterQ.enq(tuple5(bytes, pointer, continuous, direction, True));
			
				direction <= 0;
				start <= 0;
				continuous <= 0;
			end

			dataQ.deq;
			verbatimQ.enq(dataQ.first);
			caseQ.enq(0);
		end else if ( head == 2'b01 ) begin // [FORWARD NORMAL MATCH]
			if ( direction != 0 ) begin
				let p = getParameters(start, continuous);
				let address = tpl_1(p);
				let bytes   = tpl_2(p);
				let pointer = tpl_3(p);

				reqReadRefQ.enq(MemPortReq{addr:address, bytes:bytes});
				parameterQ.enq(tuple5(bytes, pointer, continuous, direction, True));
			end

			dataQ.deq;
			start <= dataQ.first;
			direction <= 1;
			continuous <= 1;
			caseQ.enq(1);
		end else if ( head == 2'b10 ) begin // [CONTINUOUS MATCH]
			if ( direction == 1 ) begin		// [FORWARD]
				continuous <= continuous + 1;
			end else if ( direction == 2 ) begin	// [REVERSE]
				start <= start - fromInteger(valueOf(Kmer));
				continuous <= continuous + 1;
			end
		end else if ( head == 2'b11 ) begin // [REVERSE NORMAL MATCH]
			if ( direction != 0 ) begin
				let p = getParameters(start, continuous);
				let address = tpl_1(p);
				let bytes   = tpl_2(p);
				let pointer = tpl_3(p);

				reqReadRefQ.enq(MemPortReq{addr:address, bytes:bytes});
				parameterQ.enq(tuple5(bytes, pointer, continuous, direction, True));
			end

			dataQ.deq;
			start <= dataQ.first;
			direction <= 2;
			continuous <= 1;
			caseQ.enq(2);
		end
	endrule
	//------------------------------------------------------------------------------------
	// [STAGE 2]
	// Decompress the sequence
	//------------------------------------------------------------------------------------
	FIFO#(Bit#(8)) decompCaseQ <- mkSizedBRAMFIFO(512);
	FIFO#(Tuple5#(Bit#(32), Bit#(32), Bit#(32), Bit#(8), Bool)) decompParameterQ <- mkSizedBRAMFIFO(512);
	FIFO#(Tuple3#(Bit#(512), Bit#(32), Bool)) decompressSub1Q               <- mkSizedBRAMFIFO(512);
        FIFO#(Tuple4#(Bit#(512), Bit#(32), Bit#(32), Bool)) decompressSub2Q     <- mkSizedBRAMFIFO(512);
	FIFO#(Bit#(8)) remainQ <- mkSizedBRAMFIFO(512);
	Reg#(Bit#(8)) remain <- mkReg(8);
	Reg#(Bit#(64)) addr <- mkReg(268435456);
	Reg#(Bool) getNewCase <- mkReg(True);
	rule decompStart( getNewCase );
		caseQ.deq;
		let c = caseQ.first;
		decompCaseQ.enq(c);
		if ( c != 0 ) begin
			parameterQ.deq;
			let param = parameterQ.first;
			decompParameterQ.enq(param);
			if ( tpl_2(param) != 0 ) remainQ.enq(0);
			getNewCase <= False;
		end
	endrule
	rule decompressMain;
		decompCaseQ.deq;
		let c = decompCaseQ.first;
	
		if ( c == 0 ) begin 	// [UNMATCH]
			verbatimQ.deq;
			let v = verbatimQ.first;

			reqWriteResultQ.enq(MemPortReq{addr:addr, bytes:64});
			resultQ.enq(zeroExtend(v));

			addr <= addr + 4;
		end else begin 	       	// [MATCH]
			refQ.deq;
			decompParameterQ.deq;
			let value = refQ.first;
			let param = decompParameterQ.first;
			let bytes 	= tpl_1(param);
			let pointer 	= tpl_2(param);
			let continuous 	= tpl_3(param);
			let direction 	= tpl_4(param);
			let first	= tpl_5(param);

			if ( direction == 2 ) value = ~(reverseBits(value));
			
			reqWriteResultQ.enq(MemPortReq{addr:addr, bytes:64});
			resultQ.enq((value >> 2) | zeroExtend(remain));
			if ( bytes > 64 ) begin
				//decompressSub1Q.enq(tuple3(value, pointer, first));
				decompCaseQ.enq(c);
				decompParameterQ.enq(tuple5(bytes - 64, pointer, continuous - 16, direction, False));
				getNewCase <= False;
			end else begin
				//decompressSub2Q.enq(tuple4(value, pointer, continuous, first));
				getNewCase <= True;
			end
		end
	endrule

	rule decompressSub1;
		decompressSub1Q.deq;
		let param = decompressSub1Q.first;
		let value 	= tpl_1(param);
		let pointer 	= tpl_2(param);
		let first 	= tpl_3(param);
		
		Bit#(512) decomp = 0;
		if ( first ) begin
			if ( pointer > 0 ) begin
				if ( pointer == 2 ) begin
					decomp = value >> 2;
				end else if ( pointer == 4 ) begin
					decomp = value >> 4;
				end else if ( pointer == 6 ) begin
					decomp = value >> 6;
				end
			end else begin
				decomp = value;
			end
		end else begin
			if ( pointer > 0 ) begin
				if ( pointer == 2 ) begin
					decomp = (value << 6) | zeroExtend(remain);
					remain <= zeroExtend(value[511:506]);
				end else if ( pointer == 4 ) begin
					decomp = (value << 4) | zeroExtend(remain);
					remain <= zeroExtend(value[511:508]);
				end else if ( pointer == 6 ) begin
					decomp = (value << 2) | zeroExtend(remain);
					remain <= zeroExtend(value[511:510]);
				end
			end else begin
				decomp = value;
			end 
		end

		addr <= addr + 64;
		resultQ.enq(decomp);
	endrule

	rule decompressSub2;
		decompressSub2Q.deq;
		let param = decompressSub2Q.first;
		let value	= tpl_1(param);
		let pointer	= tpl_2(param);
		let continuous	= tpl_3(param);
		let first 	= tpl_4(param);

		Bit#(512) decomp = 0;
		if ( first ) begin
			if ( pointer > 0 ) begin
				if ( pointer == 2 ) begin
					decomp = value >> 2;
				end else if ( pointer == 4 ) begin
					decomp = value >> 4;
				end else if ( pointer == 6 ) begin
					decomp = value >> 6;	
				end
			end else begin
				decomp = value;
			end 
		end else begin
			if ( pointer > 0 ) begin
				if ( pointer == 2 ) begin
					decomp = (value << 6) | zeroExtend(remain);
				end else if ( pointer == 4 ) begin
					decomp = (value << 4) | zeroExtend(remain);
				end else if ( pointer == 6 ) begin
					decomp = (value << 2) | zeroExtend(remain);
				end
			end else begin
				decomp = value;
			end
		end

		addr <= addr + zeroExtend(continuous * 4);
		resultQ.enq(decomp);
	endrule
	//------------------------------------------------------------------------------------
	// Interface
	//------------------------------------------------------------------------------------
	// Read the compressed genomic data
	method Action readData(Bit#(512) data);
		serializer512b32b.put(data);
	endmethod
	// Read the 2-bit encoded reference genome
	method ActionValue#(MemPortReq) reqReadRef;
		reqReadRefQ.deq;
		return reqReadRefQ.first;
	endmethod
	method Action readRef(Bit#(512) reference);
		refQ.enq(reference);
	endmethod
	// Write the decompressed 2-bit encoded sequence
	method ActionValue#(MemPortReq) reqWriteResult;
		reqWriteResultQ.deq;
		return reqWriteResultQ.first;
	endmethod
	method ActionValue#(Bit#(512)) writeResult;
		resultQ.deq;
		return resultQ.first;
	endmethod
endmodule

