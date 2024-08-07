#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <map>
using namespace std;


#define KMERLENGTH 256
#define ENCKMERBUFUNIT 32
#define ENCKMERBUFSIZE 8
#define BINARYRWUNIT 8
#define RDCREFSIZE 268435456


map<uint32_t, vector<uint64_t>> referenceOrg;
multimap<uint32_t, uint32_t> reference32b;
uint64_t referenceRdc[ENCKMERBUFSIZE * RDCREFSIZE];


uint32_t refSizeOrg = 2836860451;
uint32_t refSizeUsd = 2836860451;
uint32_t refSizeRdc = 268435456;
uint32_t refInstCnt = 0;


static inline double timeChecker( void ) {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)(tv.tv_sec) + (double)(tv.tv_usec) / 1000000;
}

void refReader( char *filename ) {
	ifstream f_data_reference(filename, ios::binary);
	for ( uint32_t i = 0; i < refSizeUsd; i ++ ) {
		uint64_t encKmer[ENCKMERBUFSIZE] = {0, };

		// Read
		for ( uint32_t j = 0; j < ENCKMERBUFSIZE; j ++ ) {
			f_data_reference.read(reinterpret_cast<char *>(&encKmer[j]), BINARYRWUNIT);
		}
		
		// Insert 256-Mers to Map
		referenceOrg[i] = vector<uint64_t>{encKmer[0], encKmer[1], encKmer[2], encKmer[3], 
						   encKmer[4], encKmer[5], encKmer[6], encKmer[7]};
		
		// Insert 32-bit Portion of 256-mers to MultiMap
		uint64_t encKmerTmp1 = encKmer[0] << 32;
		uint64_t encKmerTmp2 = encKmerTmp1 >> 32;
		uint32_t encKmer32Bits = (uint32_t)encKmerTmp2;
		reference32b.insert(make_pair(encKmer32Bits, i));
		
		if ( i % 1000000 == 0 ) {
			printf( "Reference: %u\n", i );
			fflush( stdout );
		}
	}

	f_data_reference.close();
}

void refShrinker( void ) {
	for ( uint32_t i = 0 i < refSizeRdc; i ++ ) {
		auto begin = referenceOrg.begin();

		uint64_t encKmerTmp1 = begin->second[0] << 32;
		uint64_t encKmerTmp2 = encKmerTmp2 >> 32;
		uint32_t encKmer32Bits = (uint32_t)encKmerTmp2;

		uint32_t min = refSizeUsd;
		for ( auto iter = reference32b.lower_bound(encKmer32Bits); 
			   iter != reference32b.upper_bound(encKmer32Bits); 
			   iter ++ ) {
			if ( iter->second < min ) min = iter->second;
		}

		referenceRdc[refInstCnt++] = 
	}
}


int main( void ) {
	char *filenameR = "/mnt/ephemeral/hg19hg38RefBook256Mers.bin";

	// Read reference file
	refReader( filenameR );

	return 0;
}
