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
#define STRIDE 1
#define REFINDEX 32
#define TESTSEQ 5
#define ENCKMERBUFUNIT 32
#define ENCKMERBUFSIZE 8
#define BINARYRWUNIT 8


vector<string> sequences;
unordered_map<uint64_t, uint8_t> ref_1;
unordered_map<uint64_t, uint8_t> ref_2;
unordered_map<uint64_t, uint8_t> ref_3;
unordered_map<uint64_t, uint8_t> ref_4;
unordered_map<uint64_t, uint8_t> ref_5;
unordered_map<uint64_t, uint8_t> ref_6;
unordered_map<uint64_t, uint8_t> ref_7;
unordered_map<uint64_t, uint8_t> ref_8;


uint64_t seqSizeOrg = 0;
uint64_t seqSizeCmp = 0;
uint64_t refSizeOrg = 2836860451;
uint64_t refSizeUsd = 2836860451;


static inline double timeChecker( void ) {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)(tv.tv_sec) + (double)(tv.tv_usec) / 1000000;
}

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
		// Reference 1
		f_data_reference.read(reinterpret_cast<char *>(&encKmer[0]), BINARYRWUNIT);
		ref_1.insert(make_pair(encKmer[0], 1));
		// Reference 2
		f_data_reference.read(reinterpret_cast<char *>(&encKmer[1]), BINARYRWUNIT);
		ref_2.insert(make_pair(encKmer[1], 1));
		// Reference 3
		f_data_reference.read(reinterpret_cast<char *>(&encKmer[2]), BINARYRWUNIT);
		ref_3.insert(make_pair(encKmer[2], 1));
		// Reference 4
		f_data_reference.read(reinterpret_cast<char *>(&encKmer[3]), BINARYRWUNIT);
		ref_4.insert(make_pair(encKmer[3], 1));
		// Reference 5
		f_data_reference.read(reinterpret_cast<char *>(&encKmer[4]), BINARYRWUNIT);
		ref_5.insert(make_pair(encKmer[4], 1));
		// Reference 6
		f_data_reference.read(reinterpret_cast<char *>(&encKmer[5]), BINARYRWUNIT);
		ref_6.insert(make_pair(encKmer[5], 1));
		// Reference 7
		f_data_reference.read(reinterpret_cast<char *>(&encKmer[6]), BINARYRWUNIT);
		ref_7.insert(make_pair(encKmer[6], 1));
		// Reference 8
		f_data_reference.read(reinterpret_cast<char *>(&encKmer[7]), BINARYRWUNIT);
		ref_8.insert(make_pair(encKmer[7], 1));
		if ( i % 1000000 == 0 ) {
			printf( "Reference: %ld\n", i );
			fflush( stdout );
		}
	}

	f_data_reference.close();
}

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

void compressor( void ) {
	uint64_t start = 0;
	while ( start <= sequences[TESTSEQ].size() - KMERLENGTH ) {
		string subseq = sequences[TESTSEQ].substr(start, KMERLENGTH);

		// Encode first
		uint64_t encSubseq[ENCKMERBUFSIZE] = {0, };
		encoder(subseq, encSubseq);

		// Check possible to compress
		if ( ref_1.find(encSubseq[0]) != ref_1.end() ) {
			if ( ref_2.find(encSubseq[1]) != ref_2.end() ) {
				if ( ref_3.find(encSubseq[2]) != ref_3.end() ) {
					if ( ref_4.find(encSubseq[3]) != ref_4.end() ) {
						if ( ref_5.find(encSubseq[4]) != ref_5.end() ) {
							if ( ref_6.find(encSubseq[5]) != ref_6.end() ) {
								if ( ref_7.find(encSubseq[6]) != ref_7.end() ) {
									if ( ref_8.find(encSubseq[7]) != ref_8.end() ) {
										seqSizeCmp ++;
										start += KMERLENGTH;
									} else start += STRIDE;
								} else start += STRIDE;
							} else start += STRIDE;
						} else start += STRIDE;
					} else start += STRIDE;
				} else start += STRIDE;
			} else start += STRIDE;
		} else start += STRIDE;
	}
}


int main( void ) {
	char *filenameS = "/mnt/ephemeral/hg16.fasta";
	char *filenameR = "/mnt/ephemeral/hg19hg38RefBook256Mers.bin";

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
	printf( "Reference Book [#KMER]: %ld\n", refSizeUsd );
	printf( "Reference Book [Size]: %0.4f GB\n", ((double)refSizeUsd * KMERLENGTH) / 1024 / 1024 / 1024 );
	printf( "--------------------------------------------\n" );
	printf( "SEQUENCE\n" );
	printf( "Number of Base Pairs [Original]: %ld\n", sequences[TESTSEQ].size() );
	printf( "Original File Size: %0.4f MB\n", (double)sequences[TESTSEQ].size() / 1024 / 1024 );
	printf( "--------------------------------------------\n" );
	printf( "COMPRESSION RESULT\n" );
	printf( "Number of Base Pairs [Compressed]: %ld\n", seqSizeCmp * KMERLENGTH );
	printf( "Compressed File Size [Original]: %0.4f MB\n", 
		(double)(((sequences[TESTSEQ].size() - (seqSizeCmp * KMERLENGTH)) * 8) + (seqSizeCmp * REFINDEX)) / 8 / 1024 / 1024 );
	printf( "Compressed File Size [2-b Encd]: %0.4f MB\n", 
	     	(double)(((sequences[TESTSEQ].size() - (seqSizeCmp * KMERLENGTH)) * 2) + (seqSizeCmp * REFINDEX)) / 8 / 1024 / 1024 );
	printf( "Elapsed Time: %lf\n", elapsedTime );

	return 0;
}
