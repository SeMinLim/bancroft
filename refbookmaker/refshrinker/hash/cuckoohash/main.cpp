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


unordered_map<uint32_t, uint8_t> referenceRdcX_RS;
unordered_map<uint32_t, uint8_t> referenceRdcX_JS;
vector<uint64_t> referenceRdcY;


uint64_t refSizeRead = 0;
uint64_t refSizeInst = 0;
uint64_t refSizeOrg = 2836860451;
//uint64_t refSizeRdc = 268435456;
//uint64_t refSizeRdc = 536870912;
//uint64_t refSizeRdc = 1073741824;
uint64_t refSizeRdc = 2147483648;


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

uint32_t rsHash( const char *str ) {
	uint32_t b = 378551;
	uint32_t a = 63689;
	uint32_t hash = 0;
	uint32_t i = 0;

	for ( i = 0; i < KMERLENGTH; ++str, ++i ) {
		hash = hash*a + (*str);
		a = a * b;
	}

	return hash;
}

uint32_t jsHash( const char *str ) {
	uint32_t hash = 1315423911;
	uint32_t i = 0;

	for ( i = 0; i < KMERLENGTH; ++str, ++i ) {
		hash ^= ((hash << 5) + (*str) + (hash >> 2));
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

		// Do decoding and apply RSHash
		string decoded_rs;
		decoder( encKmer, decoded_rs );
		const char *str_rs = decoded_rs.c_str();
		uint32_t hashed_rs = rsHash( str_rs );

		// Do decoding and apply JSHash
		string decoded_js;
		decoder( encKmer, decoded_js );
		const char *str_js = decoded_js.c_str();
		uint32_t hashed_js = rsHash( str_js );

		// Insert 32-bit hashed value to UNORDERED_MAP
		if ( referenceRdcX_RS.insert(make_pair(hashed_rs, 1)).second == true ) {
			for ( uint64_t cnt = 0; cnt < ENCKMERBUFSIZE; cnt ++ ) {
				referenceRdcY.push_back(encKmer[cnt]);
			}
			
			refSizeInst ++;
		} else {
			if ( referenceRdcX_JS.insert(make_pair(hashed_js, 1)).second == true ) {
				for ( uint64_t cnt = 0; cnt < ENCKMERBUFSIZE; cnt ++ ) {
					referenceRdcY.push_back(encKmer[cnt]);
				}

				refSizeInst ++;
			}
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
	char *filenameOriginal = "/mnt/ephemeral/hg19hg38RefBook256Mers.bin";
	char *filenameReduced = "/mnt/ephemeral/hg19hg38RefBook256Mers_2048M_CuckooHash.bin";

	// Shrink Original Reference File
	refShrinker( filenameOriginal );

	// Write the Reduced Reference File
	refWriter( filenameReduced );

	return 0;
}
