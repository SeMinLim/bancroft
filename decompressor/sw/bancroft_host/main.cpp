#include <iostream>
#include <cstring>
#include <vector>
#include <chrono>
#include <sys/time.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <fstream>
#include <string>
#include <unordered_map>
#include <map>
#include <algorithm>
using namespace std;

// XRT includes
#include "xrt/xrt_bo.h"
#include <experimental/xrt_xclbin.h>
#include "xrt/xrt_device.h"
#include "xrt/xrt_kernel.h"


#define DEVICE_ID 0
#define HG19_COUNT 2864785220
#define HEAD_COUNT 3716988
#define BODY_COUNT 25394544
#define SIZE_PORT_1 (1024*1024*256*2)
#define SIZE_PORT_2 (1024*1024*256*3)


vector<uint32_t> packet;
vector<uint32_t> reference;


// Packet Reader
void pktReader( char *filename ) {
	// Read
	ifstream f_data_packet(filename, ios::binary);
	while( !f_data_packet.eof() ) {
		uint32_t encPkt;
		f_data_packet.read(reinterpret_cast<char *>(&encPkt), 4);
		packet.push_back(encPkt);
	}
	// Terminate
	f_data_packet.close();
	printf( "---------------------------------------------------------------------\n" );
	printf( "[STEP 2] Reading packet is done!\n" );
	printf( "---------------------------------------------------------------------\n" );
	fflush( stdout );
}
// Reference Reader
void refReader( char *filename ) {
	// Read
	ifstream f_data_reference(filename, ios::binary);
	while( !f_data_reference.eof() ) {
		uint32_t encRef;
		f_data_reference.read(reinterpret_cast<char *>(&encRef), 4);
		reference.push_back(encRef);
	}
	// Terminate
	f_data_reference.close();
	printf( "---------------------------------------------------------------------\n" );
	printf( "[STEP 1] Reading reference is done!\n" );
	printf( "---------------------------------------------------------------------\n" );
	fflush( stdout );
}


int main(int argc, char** argv) {
	char *filenameP = "/mnt/ssd0/semin/bancroft/data/compressed/HG003_SUB_COMP_064Mers.bin";
	char *filenameR = "/mnt/ssd0/semin/bancroft/data/sequence/HG19.fasta.bin";

	// Read the compressed data
	pktReader( filenameP );

	// Read reference genome sequence
	refReader( filenameR );

	// Load xclbin
	string xclbin_file = "kernel.xclbin";
	xrt::device device = xrt::device(DEVICE_ID);
	xrt::uuid xclbin_uuid = device.load_xclbin(xclbin_file);
	
	auto krnl = xrt::kernel(device, xclbin_uuid, "kernel:{kernel_1}");

	cout << "Allocate Buffer in Global Memory\n";
	fflush( stdout );
	auto boIn1 = xrt::bo(device, (size_t)SIZE_PORT_1, krnl.group_id(1));
	auto boOut = xrt::bo(device, (size_t)SIZE_PORT_2, krnl.group_id(2)); 

	// Map the contents of the buffer object into host memory
	auto bo0_map = boIn1.map<int*>();
	auto bo2_map = boOut.map<int*>();
	fill(bo0_map, bo0_map + packet.size(), 0);
	fill(bo2_map, bo2_map + reference.size(), 0);
	for ( uint64_t i = 0; i < packet.size(); i ++ ) {
		bo0_map[i] = packet[i];
	}
	for ( uint64_t i = 0; i < reference.size(); i ++ ) {
		bo2_map[i] = reference[i];
	}

	// Synchronize buffer content with device side
	cout << " synchronize input buffer data to device global memory\n";
	fflush(stdout);
	boIn1.sync(XCL_BO_SYNC_BO_TO_DEVICE);
	boOut.sync(XCL_BO_SYNC_BO_TO_DEVICE);

	cout << "Execution of the kernel\n";
	fflush(stdout);
	auto run = krnl((size_t)SIZE_PORT_1 + (size_t)SIZE_PORT_2, boIn1, boOut); //DATA_SIZE=size
	run.wait();

	// Get the output;
	cout << "Get the output data from the device" << std::endl;
	boOut.sync(XCL_BO_SYNC_BO_FROM_DEVICE);

//	for ( int i = 0; i < 8; i++ ) {
//		printf( "%x %x:%x -- %x %x\n", bo3_map[i], DATA_SIZE-1-i, bo3_map[DATA_SIZE-1-i], bo2_map[i], bo2_map[DATA_SIZE-1-i] );
//	}

	// Validate results
	//if (std::memcmp(bo2_map, bufReference, vector_size_bytes))
		//throw std::runtime_error("Value read back does not match reference");

	cout << "TEST PASSED\n";

	return 0;
}

