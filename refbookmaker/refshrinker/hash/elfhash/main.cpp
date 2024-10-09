#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <algorithm>
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


unordered_map<uint32_t, uint8_t> referenceRdcX;
vector<uint64_t> referenceRdcY;


uint64_t refSizeRead = 0;
uint64_t refSizeInst = 0;
uint64_t refSizeOrg = 2836860451;
uint64_t refSizeRdc = 268435456;
//uint64_t refSizeRdc = 536870912;
//uint64_t refSizeRdc = 1073741824;
//uint64_t refSizeRdc = 2147483648;


void decoder( const uint64_t *encKmer, string &seqLine ) {
	for ( uint64_t i = 0; i < ENCKMERBUFSIZE; i ++ ) {
		for ( uint64_t j = 0; j < ENCKMERBUFUNIT; j ++ ) {
			uint64_t encCharT = encKmer[i] << (ENCKMERBUFUNIT - 1 - j) * 2;
			uint64_t encCharF = encCharT >> (ENCKMERBUFUNIT - 1) * 2;
			if ( encCharF == 0 ) seqLine.push_back('A');
			else if ( encCharF == 1 ) seqLine.push_back('C');
			else if ( encCharF == 2 ) seqLine.push_back('G');
			else if ( encCharF == 3 ) seqLine.push_back('T');
		}
	}
}

uint32_t elfHash( const char *str ) {
	uint32_t hash = 0;
	uint32_t x = 0;
	uint32_t i = 0;

	for ( i = 0; i < KMERLENGTH; ++str, ++i ) {
		hash = (hash << 4) + (*str);

		if ( (x = hash & 0xF0000000L) != 0 ) {
			hash ^= (x >> 24);
		}

		hash &= ~x;
	}

	return hash;
}

void refShrinker( char *filename ) {
	ifstream f_reference_original(filename, ios::binary);
	while ( refSizeInst < refSizeRdc ) {
		uint64_t encKmer[ENCKMERBUFSIZE] = {0, };

		// Read
		if ( refSizeRead == refSizeOrg ) break;
		for ( uint64_t cnt = 0; cnt < ENCKMERBUFSIZE; cnt ++ ) {
			f_reference_original.read(reinterpret_cast<char *>(&encKmer[cnt]), BINARYRWUNIT);
		}
		refSizeRead ++;

		// Do decoding and apply ELFHash
		string decoded;
		decoder( encKmer, decoded );
		const char *str = decoded.c_str();
		uint32_t hashed = elfHash( str );

		// Insert 32-bit hashed value to UNORDERED_MAP
		if ( referenceRdcX.insert(make_pair(hashed, 1)).second == true ) {
			for ( uint64_t cnt = 0; cnt < ENCKMERBUFSIZE; cnt ++ ) {
				referenceRdcY.push_back(encKmer[cnt]);
			}
			
			refSizeInst ++;
		}

		// Check the progress
		if ( refSizeRead % 1000000 == 0 ) {
			printf( "Read: %lu, Reduced Reference: %lu\n", refSizeRead, refSizeInst );
			fflush( stdout );
		}
	}

	f_reference_original.close();
	printf( "Reduced Reference: %lu\n", refSizeInst );
	printf( "Making Reduced Reference Book is Done!\n" );
	fflush( stdout );
}

void refWriter( char *filename ) {
	// Write the result
	ofstream f_reference_reduced(filename, ios::binary);
	for ( uint64_t i = 0; i < refSizeInst;  i ++ ) {
		for ( uint64_t j = 0; j < ENCKMERBUFSIZE; j ++ ) {
			f_reference_reduced.write(reinterpret_cast<char *>(&referenceRdcY[i*ENCKMERBUFSIZE + j]), BINARYRWUNIT);
		}
	}

	f_reference_reduced.close();

	printf( "Writing the Reduced Reference File is Done!\n" );
	fflush( stdout );
}


int main( void ) {
	char *filenameOriginal = "hg19hg38Reference256Mers.bin";
	char *filenameReduced = "hg19hg38Reference256Mers256MELFHash.bin";

	// Shrink Original Reference File
	refShrinker( filenameOriginal );

	// Write the Reduced Reference File
	refWriter( filenameReduced );

	return 0;
}
