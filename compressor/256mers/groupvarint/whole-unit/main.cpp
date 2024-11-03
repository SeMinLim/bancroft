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


string sequence;
vector<uint32_t> index;
map<pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, uint64_t>>>>>>>, 
	 uint32_t> reference;


uint64_t seqSizeOrg = 0;
uint64_t seqSizeCmpN = 0;
uint64_t seqSizeCmpI = 0;
uint64_t seqSizeCmpP = 0;
uint64_t seqSizeRmnd = 0;


// Reference: hg19From2
//uint64_t refSizeOrg = 10976427;
//uint64_t refSizeUsd = 10976427;
// Reference: hg19From1
//uint64_t refSizeOrg = 2849207900;
//uint64_t refSizeUsd = 268435456;
//uint64_t refSizeUsd = 536870912;
//uint64_t refSizeUsd = 1073741824;
//uint64_t refSizeUsd = 2147483648;
uint64_t refSizeUsd = 2849207900;
// Reference: hg19From1 [32LSB]
//uint64_t refSizeUsd = 268435456;
//uint64_t refSizeUsd = 536870912;
//uint64_t refSizeUsd = 1073741824;
//uint64_t refSizeUsd = 1131364271;
// Reference: hg19From1 [RS Hash]
//uint64_t refSizeUsd = 268435456;
//uint64_t refSizeUsd = 536870912;
//uint64_t refSizeUsd = 1073741824;
//uint64_t refSizeUsd = 2082620638;
// Reference: hg19From1 [JS Hash]
//uint64_t refSizeUsd = 268435456;
//uint64_t refSizeUsd = 536870912;
//uint64_t refSizeUsd = 1073741824;
//uint64_t refSizeUsd = 1968745020;
// Reference: hg19From1 [Cuckoo Hash]
//uint64_t refSizeUsd = 268435456;
//uint64_t refSizeUsd = 536870912;
//uint64_t refSizeUsd = 1073741824;
//uint64_t refSizeUsd = 2147483648;
// Reference: hg19Shrinked4KB
//uint64_t refSizeOrg = 1428002113;
//uint64_t refSizeUsd = 268435456;
//uint64_t refSizeUsd = 536870912;
//uint64_t refSizeUsd = 1073741824;
//uint64_t refSizeUsd = 1428002113;
// Reference: hg19+hg38
//uint64_t refSizeOrg = 2836860451;
//uint64_t refSizeUsd = 268435456;
//uint64_t refSizeUsd = 536870912;
//uint64_t refSizeUsd = 1073741824;
//uint64_t refSizeUsd = 2147483648;
//uint64_t refSizeUsd = 2836860451;
// Reference: hg19+hg38 [32LSB]
//uint64_t refSizeUsd = 268435456;
//uint64_t refSizeUsd = 536870912;
//uint64_t refSizeUsd = 1073741824;
//uint64_t refSizeUsd = 1128110488;
// Reference: hg19+hg38 [RS Hash]
//uint64_t refSizeUsd = 268435456;
//uint64_t refSizeUsd = 536870912;
//uint64_t refSizeUsd = 1073741824;
//uint64_t refSizeUsd = 2076244715;
// Reference: hg19+hg38 [JS Hash]
//uint64_t refSizeUsd = 268435456;
//uint64_t refSizeUsd = 536870912;
//uint64_t refSizeUsd = 1073741824;
//uint64_t refSizeUsd = 1963010781;
// Reference: hg19+hg38 [Cuckoo Hash]
//uint64_t refSizeUsd = 268435456;
//uint64_t refSizeUsd = 536870912;
//uint64_t refSizeUsd = 1073741824;
//uint64_t refSizeUsd = 2147483648;
// Reference: hg19hg38
//uint64_t refSizeOrg = 2825518939;
//uint64_t refSizeUsd = 268435456;
//uint64_t refSizeUsd = 536870912;
//uint64_t refSizeUsd = 1073741824;
//uint64_t refSizeUsd = 2147483648;
//uint64_t refSizeUsd = 2825518939;
// Reference: hg19hg38 [32LSB]
//uint64_t refSizeUsd = 268435456;
//uint64_t refSizeUsd = 536870912;
//uint64_t refSizeUsd = 1073741824;
//uint64_t refSizeUsd = 1127665898;
// Reference: hg19hg38 [RS Hash]
//uint64_t refSizeUsd = 268435456;
//uint64_t refSizeUsd = 536870912;
//uint64_t refSizeUsd = 1073741824;
//uint64_t refSizeUsd = 2070382593;
// Reference: hg19hg38 [JS Hash]
//uint64_t refSizeUsd = 268435456;
//uint64_t refSizeUsd = 536870912;
//uint64_t refSizeUsd = 1073741824;
//uint64_t refSizeUsd = 1957754995;
// Reference: hg19+hg38 [Cuckoo Hash]
//uint64_t refSizeUsd = 268435456;
//uint64_t refSizeUsd = 536870912;
//uint64_t refSizeUsd = 1073741824;
//uint64_t refSizeUsd = 2147483648;


// Required Functions
// Time Checker
static inline double timeChecker( void ) {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)(tv.tv_sec) + (double)(tv.tv_usec) / 1000000;
}
// Sequence File Reader
void seqReader( char *filename ) {
	string seqLine;
	// Read
	ifstream f_data_sequence(filename);
	while ( getline(f_data_sequence, seqLine) ) {
		if ( seqLine[0] != '>' ) {
			sequence += seqLine;
			seqSizeOrg += seqLine.size();
		}
	}
	// Terminate
	f_data_sequence.close();
	printf( "---------------------------------------------------------------------\n" );
	printf( "[STEP 1] Reading sequence fasta file is done!\n" );
	printf( "---------------------------------------------------------------------\n" );
	fflush( stdout );
}
// Reference File Reader
void refReader( char *filename ) {
	ifstream f_data_reference(filename, ios::binary);
	for ( uint64_t i = 0; i < refSizeUsd; i ++ ) {
		uint64_t encKmer[ENCKMERBUFSIZE + 1] = {0, };
		// Read 256-mer and index
		for ( uint64_t j = 0; j < ENCKMERBUFSIZE + 1; j ++ ) {
			f_data_reference.read(reinterpret_cast<char *>(&encKmer[j]), BINARYRWUNIT);
		}
		// Insert 256-mer and index to Map
		if ( reference.insert(make_pair(make_pair(encKmer[0], make_pair(encKmer[1], make_pair(encKmer[2], 
				      make_pair(encKmer[3], make_pair(encKmer[4], make_pair(encKmer[5], 
				      make_pair(encKmer[6], encKmer[7]))))))), (uint32_t)encKmer[8])).second == false ) {
			printf( "There's a problem on reference code book...\n" );
			fflush( stdout );
			exit(1);
		}
		// Check the progress
		if ( i % 1000000 == 0 ) {
			printf( "[STEP 2] Reading reference is processing...[%lu/%lu]\n", i, refSizeUsd );
			fflush( stdout );
		}
	}
	// Terminate
	f_data_reference.close();
	printf( "[STEP 2] Reading reference is done!\n" );
	printf( "---------------------------------------------------------------------\n" );
	fflush( stdout );
}
// 2-bit Encoder
void encoder( string seqLine, uint64_t *encKmer ) {
	for ( uint64_t i = 0; i < ENCKMERBUFSIZE; i ++ ) {
		encKmer[i] = 0;
		for ( uint64_t j = 0; j < ENCKMERBUFUNIT; j ++ ) {
			if ( seqLine[ENCKMERBUFUNIT*i + j] == 'A' ) {
				encKmer[i] = ((uint64_t)0 << 2 * j) | encKmer[i];
			} else if ( seqLine[ENCKMERBUFUNIT*i + j] == 'C' ) {
				encKmer[i] = ((uint64_t)1 << 2 * j) | encKmer[i];
			} else if ( seqLine[ENCKMERBUFUNIT*i + j] == 'G' ) {
				encKmer[i] = ((uint64_t)2 << 2 * j) | encKmer[i];
			} else if ( seqLine[ENCKMERBUFUNIT*i + j] == 'T' ) {
				encKmer[i] = ((uint64_t)3 << 2 * j) | encKmer[i];
			}
		}
	}
}
// 2-bit Decoder
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
// Compressor
void compressor( const uint64_t stride ) {
	uint64_t start = 0;
	while ( start <= sequence.size() - KMERLENGTH ) {
		string subseq = sequence.substr(start, KMERLENGTH);
		// Encode first
		uint64_t encSubseq[ENCKMERBUFSIZE] = {0, };
		encoder(subseq, encSubseq);
		// Check possible to compress
		if ( reference.find(make_pair(encSubseq[0], make_pair(encSubseq[1], make_pair(encSubseq[2], 
				    make_pair(encSubseq[3], make_pair(encSubseq[4], make_pair(encSubseq[5], 
				    make_pair(encSubseq[6], encSubseq[7])))))))) != reference.end() ) {
			// Possible to compress, then put index to the vector first
			index.push_back(reference.at(make_pair(encSubseq[0], make_pair(encSubseq[1], make_pair(encSubseq[2], 
						     make_pair(encSubseq[3], make_pair(encSubseq[4], make_pair(encSubseq[5], 
						     make_pair(encSubseq[6], encSubseq[7])))))))));
			// Compare the current index to the previous one
			if ( seqSizeCmpP != 0 ) {
				if ( index[seqSizeCmpP] == (index[seqSizeCmpP - 1] + KMERLENGTH) ) {
					seqSizeCmpI ++;
				}
			}	
			seqSizeCmpP ++;
			start += KMERLENGTH;
		} else {
			seqSizeCmpN ++;
			start += stride;
		}
		// Check the progress
		if ( start % 1000000 == 0 ) {
			printf( "[STEP 3] Compressing the sequences is processing...[%lu/%lu]\n", start, sequence.size() );
			fflush( stdout );
		}
	}
	// Handle remainder
	uint64_t remainder = sequence.size() - start;
	if ( remainder > 0 ) seqSizeRmnd += (remainder * 2) + 1;
	// Terminate
	printf( "[STEP 3] Compressing the sequences is done!\n" );
	printf( "---------------------------------------------------------------------\n" );
	fflush( stdout );
}


int main( int argc, char **argv ) {
	char *filenameS = "/mnt/ephemeral/hg002.fasta";
	char *filenameR = "/mnt/ephemeral/hg19Reference256MersFrom1IndexIncluded.bin";

	// Read sequence file
	seqReader( filenameS );

	// Read reference file
	refReader( filenameR );

	// Compression
	for ( uint64_t stride = 1; stride < 512; stride = stride * 2 ) {
		// Variable initialization
		seqSizeCmpN = 0;
		seqSizeCmpI = 0;
		seqSizeCmpP = 0;
		seqSizeRmnd = 0;
		// Compress
		double processStart = timeChecker();
		compressor( stride );
		double processFinish = timeChecker();
		double elapsedTime = processFinish - processStart;
		// Results
		uint64_t refCompN = (2 + (stride * 2)) * seqSizeCmpN;
		uint64_t refCompP = (2 * seqSizeCmpP) + (32 * (seqSizeCmpP - seqSizeCmpI));
		printf( "REFERENCE\n" );
		printf( "The Length of K-Mer: %lu\n", KMERLENGTH );
		printf( "The Number of K-Mer: %lu\n", refSizeUsd );
		printf( "---------------------------------------------------------------------\n" );
		printf( "SEQUENCE\n" );
		printf( "The Number of Base Pair : %lu\n", seqSizeOrg );
		printf( "The Original File Size  : %0.4f MB\n", (double)seqSizeOrg / 1024 / 1024 / 4 );
		printf( "---------------------------------------------------------------------\n" );
		printf( "COMPRESSION RESULT\n" );
		printf( "Stride                  : %lu\n", stride );
		printf( "The Number of Base Pair : %lu\n", seqSizeCmpP * KMERLENGTH );
		printf( "The Compressed File Size: %0.4f MB\n", 
		     	(double)(refCompN + refCompP + seqSizeRmnd) / 8 / 1024 / 1024 );
		printf( "Sequential Percentage   : %0.4f\n", (double)((double)seqSizeCmpI / (double)seqSizeCmpP) * (double)100.00 );
		printf( "Elapsed Time: %lf\n", elapsedTime );
		printf( "---------------------------------------------------------------------\n" );
	}

	return 0;
}
