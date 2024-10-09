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


vector<uint64_t> referenceRdc;


uint64_t refSizeRead = 0;
uint64_t refSizeInst = 0;
uint64_t refSizeOrg = 2825518939;
uint64_t refSizeRdc = 268435456;
//uint64_t refSizeRdc = 536870912;
//uint64_t refSizeRdc = 1073741824;
//uint64_t refSizeRdc = 2147483648;


void refShrinker( char *filename ) {
	ifstream f_reference_original(filename, ios::binary);
	while ( refSizeInst < refSizeRdc ) {
		uint64_t encKmer[ENCKMERBUFSIZE] = {0, };

		// Read
		for ( uint64_t cnt = 0; cnt < ENCKMERBUFSIZE; cnt ++ ) {
			f_reference_original.read(reinterpret_cast<char *>(&encKmer[cnt]), BINARYRWUNIT);
		}
		refSizeRead ++;

		// Insert 256-mers to the vector
		for ( uint64_t cnt = 0; cnt < ENCKMERBUFSIZE; cnt ++ ) {
			referenceRdc.push_back(encKmer[cnt]);
		}
		refSizeInst ++;

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
			f_reference_reduced.write(reinterpret_cast<char *>(&referenceRdc[i*ENCKMERBUFSIZE + j]), BINARYRWUNIT);
		}
	}

	f_reference_reduced.close();

	printf( "Writing the Reduced Reference File is Done!\n" );
	fflush( stdout );
}


int main( void ) {
	char *filenameOriginal = "hg19hg38Reference256Mers.bin";
	char *filenameReduced = "hg19hg38Reference256Mers256MVanila.bin";

	// Shrink Original Reference File
	refShrinker( filenameOriginal );

	// Write the Reduced Reference File
	refWriter( filenameReduced );

	return 0;
}
