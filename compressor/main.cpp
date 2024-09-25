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
// Stride == 128
#define STRIDE 128
#define REFCOMPINDEX_N 129
// Stride == 64
//#define STRIDE 64
//#define REFCOMPINDEX_N 65
// Stride == 32
//#define STRIDE 32
//#define REFCOMPINDEX_N 33
// Stride == 16
//#define STRIDE 16
//#define REFCOMPINDEX_N 17
// Stride == 8
//#define STRIDE 8
//#define REFCOMPINDEX_N 9
// Stride == 4
//#define STRIDE 4
//#define REFCOMPINDEX_N 5
// Stride == 2
//#define STRIDE 2
//#define REFCOMPINDEX_N 3
// Stride == 1
//#define STRIDE 1
//#define REFCOMPINDEX_N 2
#define REFCOMPINDEX_P 33

#define ENCKMERBUFUNIT 32
#define ENCKMERBUFSIZE 8
#define BINARYRWUNIT 8


vector<string> sequences;
map<pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, uint64_t>>>>>>>, uint32_t> reference;


uint64_t seqSizeOrg = 0;
uint64_t seqSizeCmpN = 0;
uint64_t seqSizeCmpP = 0;

uint64_t refSizeOrg = 2836860451;
//uint64_t refSizeUsd = 268435456;
//uint64_t refSizeUsd = 536870912;
//uint64_t refSizeUsd = 1073741824;
//uint64_t refSizeUsd = 2147483648;
//uint64_t refSizeUsd = 2836860451;
//uint64_t refSizeUsd = 1128110488;
//uint64_t refSizeUsd = 2076244715;
uint64_t refSizeUsd = 1963010781;


// Time Checker
static inline double timeChecker( void ) {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)(tv.tv_sec) + (double)(tv.tv_usec) / 1000000;
}

// File Readers
void seqReader( char *filename ) {
	string seqLine;

	ifstream f_data_sequences(filename);
	while ( getline(f_data_sequences, seqLine) ) {
		if ( seqLine[0] != '>' ) {
			sequences.push_back(seqLine);
			seqSizeOrg += seqLine.size();
		}
	}

	f_data_sequences.close();
}
void refReader( char *filename ) {
	ifstream f_data_reference(filename, ios::binary);
	for ( uint64_t i = 0; i < refSizeUsd; i ++ ) {
		uint64_t encKmer[ENCKMERBUFSIZE] = {0, };
		// Read
		for ( uint64_t j = 0; j < ENCKMERBUFSIZE; j ++ ) {
			f_data_reference.read(reinterpret_cast<char *>(&encKmer[j]), BINARYRWUNIT);
		}
		// Insert 256-Mers to Map
		if ( reference.insert(make_pair(make_pair(encKmer[0], make_pair(encKmer[1], make_pair(encKmer[2], 
				      make_pair(encKmer[3], make_pair(encKmer[4], make_pair(encKmer[5], 
				      make_pair(encKmer[6], encKmer[7]))))))), i)).second == false ) {
			printf( "There's a problem on reference code book...\n" );
			fflush( stdout );
			exit(1);
		}
		// Check the progress
		if ( i % 1000000 == 0 ) {
			printf( "Reference: %lu\n", i );
			fflush( stdout );
		}
	}

	f_data_reference.close();
}

// Encoder & Decoder
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
void compressor( void ) {
	for ( uint64_t seqIdx = 0; seqIdx < sequences.size() - 1; seqIdx ++ ) {
		uint64_t start = 0;
		while ( start <= sequences[seqIdx].size() - KMERLENGTH ) {
			string subseq = sequences[seqIdx].substr(start, KMERLENGTH);

			// Encode first
			uint64_t encSubseq[ENCKMERBUFSIZE] = {0, };
			encoder(subseq, encSubseq);

			// Check possible to compress
			if ( reference.find(make_pair(encSubseq[0], make_pair(encSubseq[1], make_pair(encSubseq[2], 
					    make_pair(encSubseq[3], make_pair(encSubseq[4], make_pair(encSubseq[5], 
					    make_pair(encSubseq[6], encSubseq[7])))))))) != reference.end() ) {
				seqSizeCmpP ++;
				start += KMERLENGTH;
			} else {
				seqSizeCmpN ++;
				start += STRIDE;
			}
		}
		printf( "Compressing #%lu Sequences is Done!\n", seqIdx );
		fflush( stdout );
	}
}


int main( void ) {
	char *filenameS = "/mnt/ephemeral/hg16.fasta";
	char *filenameR = "/mnt/ephemeral/hg19hg38RefBook256Mers_2048M_JSHash.bin";

	// Read sequence file
	seqReader( filenameS );

	// Read reference file
	refReader( filenameR );

	// Compression
	double processStart = timeChecker();
	compressor();
	double processFinish = timeChecker();
	double elapsedTime = processFinish - processStart;

	printf( "--------------------------------------------\n" );
	printf( "REFERENCE\n" );
	printf( "The Length of K-Mer: %lu\n", KMERLENGTH );
	printf( "The Number of K-Mer: %lu\n", refSizeUsd );
	printf( "--------------------------------------------\n" );
	printf( "COMPRESSION RESULT\n" );
	printf( "The Original File Size  : %0.4f MB\n", (double)seqSizeOrg / 1024 / 1024 / 4 );
	printf( "The Compressed File Size: %0.4f MB\n", 
	     	(double)((seqSizeCmpN * REFCOMPINDEX_N) + (seqSizeCmpP * REFCOMPINDEX_P)) / 8 / 1024 / 1024 );
	printf( "Elapsed Time: %lf\n", elapsedTime );

	return 0;
}