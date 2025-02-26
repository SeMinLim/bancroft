#include <sys/time.h>
#include <unistd.h>
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
#include <algorithm>
using namespace std;


#define CHROMOSOMEUNIT 1


string sequence;
vector<string> sequences;


uint64_t seqSizeOrg = 0;
uint64_t seqSizeCmp = 0;
uint64_t seqSizeRmnd = 0;
uint64_t portion_1 = 0;
uint64_t portion_2 = 0;
uint64_t portion_t = 0;


// Required Functions
// Time Checker
static inline double timeChecker( void ) {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)(tv.tv_sec) + (double)(tv.tv_usec) / 1000000;
}
// Sequence FASTQ File Reader
void seqReaderFASTQ( char *filename ) {
	uint64_t counter = 0;
	string seqLine;
	// Read
	ifstream f_data_sequences(filename);
	while ( getline(f_data_sequences, seqLine) ) {
		if ( counter == 0 ) {
			counter ++;
		} else if ( counter == 1 ) {
			counter ++;
		} else if ( counter == 2 ) {
			counter ++;
		} else if ( counter == 3 ) {
			if ( CHROMOSOMEUNIT ) sequences.push_back(seqLine);
			else sequence += seqLine;
			seqSizeOrg += (seqLine.size() * 8);
			counter = 0;
		}
	}
	// Terminate
	f_data_sequences.close();
	printf( "---------------------------------------------------------------------\n" );
	printf( "[STEP 1] Reading sequence fastq file is done!\n" );
	printf( "---------------------------------------------------------------------\n" );
	fflush( stdout );
}
// Compressor [CHROMOSOME UNIT]
void compressor_unit_ch( const uint64_t stride ) {
	for ( uint64_t seqIdx = 0; seqIdx < sequences.size(); seqIdx ++ ) {
		char prevChar_1 = 'I';
		char prevChar_2 = 'I';
		string prevSubseq;
		uint64_t start = 0;
		while ( (double)start <= (double)sequences[seqIdx].size() - (double)stride ) {
			// Get subsequence
			string subseq = sequences[seqIdx].substr(start, stride);
			if ( prevSubseq.compare(subseq) == 0 ) {
				seqSizeCmp += 2;
				portion_1 ++;
			} else if ( prevChar_1 == subseq[0] ) {
				uint64_t counterCmp = 0;
				for ( uint64_t subseqIdx = 0; subseqIdx < stride; subseqIdx ++ ) {
					if ( prevChar_1 == subseq[subseqIdx] ) counterCmp ++;
				}
				if ( counterCmp == stride ) {
					seqSizeCmp += 2;
					portion_2 ++;
				}
			} else if ( prevChar_2 == subseq[0] ) {
				uint64_t counterCmp = 0;
				for ( uint64_t subseqIdx = 0; subseqIdx < stride; subseqIdx ++ ) {
					if ( prevChar_2 == subseq[subseqIdx] ) counterCmp ++;
				}
				if ( counterCmp == stride ) {
					seqSizeCmp += 2;
					portion_2 ++;
				}

			} else {
				seqSizeCmp += (stride * 6) + 2;
			}
			prevChar_1 = subseq[0];
			prevChar_2 = subseq[stride-1];
			prevSubseq = subseq;
			start += stride;
			portion_t ++;
		}
		// Handle remainder
		uint64_t remainder = sequences[seqIdx].size() - start;
		if ( remainder > 0 ) seqSizeRmnd += (remainder * 6) + 2;
		// Check the progress
		if ( seqIdx % 1000 == 0 ) {
			printf( "[STEP 3] Compressing the sequences is processing...[%lu/%lu]\n", seqIdx, sequences.size() );
			fflush( stdout );
		}
	}
	// Terminate
	printf( "[STEP 3] Compressing the sequences is done!\n" );
	printf( "---------------------------------------------------------------------\n" );
	fflush( stdout );
}
// Compressor [WHOLE UNIT]
void compressor_unit_wh( const uint64_t stride ) {
	char prevChar_1 = 'I';
	char prevChar_2 = 'I';
	string prevSubseq;
	uint64_t start = 0;
	while ( start <= sequence.size() - stride ) {
		// Get subsequence
		string subseq = sequence.substr(start, stride);
		if ( prevSubseq.compare(subseq) == 0 ) {
			seqSizeCmp += 2;
			portion_1 ++;
		} else if ( prevChar_1 == subseq[0] ) {
			uint64_t counterCmp = 0;
			for ( uint64_t subseqIdx = 0; subseqIdx < stride; subseqIdx ++ ) {
				if ( prevChar_1 == subseq[subseqIdx] ) counterCmp ++;
			}
			if ( counterCmp == stride ) {
				seqSizeCmp += 2;
				portion_2 ++;
			}
		} else if ( prevChar_2 == subseq[0] ) {
			uint64_t counterCmp = 0;
			for ( uint64_t subseqIdx = 0; subseqIdx < stride; subseqIdx ++ ) {
				if ( prevChar_2 == subseq[subseqIdx] ) counterCmp ++;
			}
			if ( counterCmp == stride ) {
				seqSizeCmp += 2;
				portion_2 ++;
			}
		} else {
			seqSizeCmp += (stride * 6) + 2;
		}
		prevChar_1 = subseq[0];
		prevChar_2 = subseq[stride-1];
		prevSubseq = subseq;
		start += stride;
		portion_t ++;
		// Check the progress
		if ( start % 1000000 == 0 ) {
			printf( "[STEP 3] Compressing the sequences is processing...[%lu/%lu]\n", start, sequence.size() );
			fflush( stdout );
		}
	}
	// Handle remainder
	uint64_t remainder = sequence.size() - start;
	if ( remainder > 0 ) seqSizeRmnd += (remainder * 6) + 2;
	// Terminate
	printf( "[STEP 3] Compressing the sequences is done!\n" );
	printf( "---------------------------------------------------------------------\n" );
	fflush( stdout );

}


int main( int argc, char **argv ) {
	char *filenameS = "/mnt/ssd0/semin/bancroft/data/sequences/fastq/HG002.fastq";

	// Read sequence file
	seqReaderFASTQ( filenameS );

	// Compression
	uint64_t base = 85;
	for ( uint64_t iter = 0; iter < 4; iter ++ ) {
		uint64_t stride = 0;
		// Set the stride
		if ( iter == 0 ) stride = base;
		else stride = (base * 4) + 1;
		// Update the base
		base = stride;
		// Variable initialization
		seqSizeCmp = 0;
		seqSizeRmnd = 0;
		// Compress
		double processStart = timeChecker();
		if ( CHROMOSOMEUNIT ) compressor_unit_ch( stride );
		else compressor_unit_wh( stride );
		double processFinish = timeChecker();
		double elapsedTime = processFinish - processStart;
		// Results
		printf( "QUALITY SCORE\n" );
		printf( "The Original File Size   : %0.4f B\n", (double)seqSizeOrg / 8.00 );
		printf( "The Original File Size   : %0.4f KB\n", (double)seqSizeOrg / 1024.00 / 8.00 );
		printf( "The Original File Size   : %0.4f MB\n", (double)seqSizeOrg / 1024.00 / 1024.00 / 8.00 );
		printf( "---------------------------------------------------------------------\n" );
		printf( "COMPRESSION RESULT\n" );
		printf( "Stride                   : %lu\n", stride );
		printf( "The Compressed File Size : %0.4f B\n", (double)seqSizeCmp / 8.00 );
		printf( "The Compressed File Size : %0.4f KB\n", (double)seqSizeCmp / 1024.00 / 8.00 );
		printf( "The Compressed File Size : %0.4f MB\n", (double)seqSizeCmp / 1024.00 / 1024.00 / 8.00 );
		printf( "---------------------------------------------------------------------\n" );
		printf( "Previous Subsequence     : %0.4f\n", (double)portion_1 / (double)portion_t );
		printf( "Previous Character       : %0.4f\n", (double)portion_2 / (double)portion_t );
		printf( "Elapsed Time: %lf\n", elapsedTime );
		printf( "---------------------------------------------------------------------\n" );
	}

	return 0;
}
