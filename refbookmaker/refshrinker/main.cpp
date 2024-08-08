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
vector<uint64_t> referenceRdc;


uint32_t refSizeOrg = 2836860451;
uint32_t refSizeUsd = 2836860451;
uint32_t refSizeRdc = 268435456;
uint32_t refInstCnt = 0;


uint64_t verification_x[ENCKMERBUFSIZE];
uint64_t verification_y[ENCKMERBUFSIZE];


void refReader( char *filename ) {
	ifstream f_reference_original(filename, ios::binary);
	for ( uint32_t i = 0; i < refSizeUsd; i ++ ) {
		uint64_t encKmer[ENCKMERBUFSIZE] = {0, };

		// Read
		for ( uint32_t j = 0; j < ENCKMERBUFSIZE; j ++ ) {
			f_reference_original.read(reinterpret_cast<char *>(&encKmer[j]), BINARYRWUNIT);
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

	f_reference_original.close();
}

void refShrinker( void ) {
	for ( uint32_t i = 0; i < refSizeRdc; i ++ ) {
		auto begin = referenceOrg.begin();

		// Verification Purpose
		if ( i == 0 ) {
			for ( uint32_t j = 0; j < ENCKMERBUFSIZE; j ++ ) {
				verification_x[j] = begin->second[j];
			}
		}

		uint64_t encKmerTmp1 = begin->second[0] << 32;
		uint64_t encKmerTmp2 = encKmerTmp1 >> 32;
		uint32_t encKmer32Bits = (uint32_t)encKmerTmp2;

		vector<uint32_t> index;
		for ( auto iter = reference32b.lower_bound(encKmer32Bits); 
			   iter != reference32b.upper_bound(encKmer32Bits); 
			   iter ++ ) {
			index.push_back(iter->second);
		}

		for ( uint32_t j = 0; j < index.size(); j ++ ) {
			referenceOrg.erase(index[j]);
		}
		
		for ( uint32_t j = 0; j < ENCKMERBUFSIZE; j ++ ) {
			referenceRdc.push_back(begin->second[j]);
		}
	}
}

void refWriter( char *filename ) {
	ofstream f_reference_reduced(filename, ios::binary);
	for ( uint32_t i = 0; i < refSizeRdc;  i ++ ) {
		for ( uint32_t j = 0; j < ENCKMERBUFSIZE; j ++ ) {
			f_reference_reduced.write(reinterpret_cast<char *>(&referenceRdc[i*ENCKMERBUFSIZE + j]), BINARYRWUNIT);
		}
	}

	f_reference_reduced.close();
}


int main( void ) {
	char *filenameOriginal = "/mnt/ephemeral/hg19hg38RefBook256Mers.bin";
	char *filenameReduced = "/mnt/ephemeral/hg19hg38RefBook256MersReduced.bin";

	// Read Original Reference File
	refReader( filenameOriginal );

	// Shrink Original Reference File
	refShrinker();

	// Write the Reduced Reference File
	refWriter( filenameReduced );

	// Verification
	ifstream f_verification(filenameReduced, ios::binary);
	for ( uint32_t i = 0; i < ENCKMERBUFSIZE; i ++ ) {
		f_verification.read(reinterpret_cast<char *>(&verification_y[i]), BINARYRWUNIT);
	}

	uint32_t cnt = 0;
	for ( uint32_t i = 0; i < ENCKMERBUFSIZE; i ++ ) {
		if ( verification_x[i] != verification_y[i] ) cnt ++;
	}
	printf( "%u\n", cnt );

	return 0;
}
