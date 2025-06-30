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


vector<string> qscores;
map<char, uint64_t> qscoreHash;


char mostFreqChar_1;
char mostFreqChar_2;
uint64_t qscoreSizeOrgn = 0;
uint64_t qscoreSizeCmpr = 0;
uint64_t qscoreSizeRmnd = 0;
uint64_t portion_1 = 0;
uint64_t portion_2 = 0;
uint64_t portion_3 = 0;
uint64_t portion_4 = 0;
uint64_t portion_5 = 0;
uint64_t portion_t = 0;


// Required Functions
// Time Checker
static inline double timeChecker( void ) {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)(tv.tv_sec) + (double)(tv.tv_usec) / 1000000;
}
// Descending Order Sorter
bool descendingOrder( pair<char, uint64_t>& x, pair<char, uint64_t>& y ) {
	return x.second > y.second;
}
// FASTQ File Reader
void readerFASTQ( char *filename ) {
	uint64_t counter = 0;
	string qscoreLine;
	// Read
	ifstream f_data_qscores(filename);
	while ( getline(f_data_qscores, qscoreLine) ) {
		if ( counter == 0 ) {
			counter ++;
		} else if ( counter == 1 ) {
			counter ++;
		} else if ( counter == 2 ) {
			counter ++;
		} else if ( counter == 3 ) {
			qscores.push_back(qscoreLine);
			qscoreSizeOrgn += qscoreLine.size();
			counter = 0;
		}
	}
	// Terminate
	f_data_qscores.close();
	printf( "---------------------------------------------------------------------\n" );
	printf( "[STEP 1] Reading FASTQ file is done!\n" );
	printf( "---------------------------------------------------------------------\n" );
	fflush( stdout );
}
// Hash Table Maker To Get Most Frequently Occurred Character
void qscoreHashTableMaker( void ) {
	for ( uint64_t idx = 0; idx < qscores.size(); idx ++ ) {
		for ( uint64_t c = 0; c < qscores[idx].size(); c ++ ) {
			if ( qscoreHash.insert(make_pair(qscores[idx][c], 1)).second == false ) {
				qscoreHash.at(qscores[idx][c]) += 1;
			}
		}
	}
	vector<pair<char, uint64_t>> qscoreHash_vector(qscoreHash.begin(), qscoreHash.end());
	sort(qscoreHash_vector.begin(), qscoreHash_vector.end(), descendingOrder);
	mostFreqChar_1 = qscoreHash_vector[0].first;
	mostFreqChar_2 = qscoreHash_vector[1].first;
	printf( "[STEP 2] Getting the most frequently occurred character is done!\n" );
	printf( "---------------------------------------------------------------------\n" );
	fflush( stdout );
}
// Compressor [CHROMOSOME UNIT]
void compressor_unit_ch( const uint64_t stride ) {
	for ( uint64_t idx = 0; idx < qscores.size(); idx ++ ) {
		char prevChar_1 = mostFreqChar_1;
		char prevChar_2 = mostFreqChar_2;
		string prevSubseq;
		uint64_t start = 0;
		while ( (double)start <= (double)qscores[idx].size() - (double)stride ) {
			// Get subsequence
			string subseq = qscores[idx].substr(start, stride);
			if ( prevSubseq.compare(subseq) == 0 ) {
				portion_1 ++;
				qscoreSizeCmpr += 2;
			} else {
				if ( prevChar_1 == subseq[0] ) {
					uint64_t counterCmp = 0;
					for ( uint64_t subseqIdx = 0; subseqIdx < stride; subseqIdx ++ ) {
						if ( prevChar_1 == subseq[subseqIdx] ) counterCmp ++;
					}
					if ( counterCmp == stride ) {
						portion_2 ++;
						qscoreSizeCmpr += 2;
					} else {
						portion_4 ++;
						qscoreSizeCmpr += (stride * 6) + 2;
					}
				} else if ( prevChar_2 == subseq[0] ) {
					uint64_t counterCmp = 0;
					for ( uint64_t subseqIdx = 0; subseqIdx < stride; subseqIdx ++ ) {
						if ( prevChar_2 == subseq[subseqIdx] ) counterCmp ++;
					}
					if ( counterCmp == stride ) {
						portion_3 ++;
						qscoreSizeCmpr += 2;
					} else {
						portion_4 ++;
						qscoreSizeCmpr += (stride * 6) + 2;
					}
				} else {
					portion_4 ++;
					qscoreSizeCmpr += (stride * 6) + 2;
				}
			}
			// Parameter update
			if ( subseq[0] == subseq[stride-1] ) {
				if ( subseq[0] == mostFreqChar_1 ) prevChar_1 = mostFreqChar_2;
				else prevChar_1 = mostFreqChar_1;
			}
			prevChar_2 = subseq[stride-1];
			prevSubseq = subseq;
			start += stride;
			portion_t ++;
		}
		// Handle remainder
		uint64_t remainder = qscores[idx].size() - start;
		if ( remainder > 0 ) {
			qscoreSizeRmnd += ((remainder * 6) + 2);
			portion_5 += remainder;
		}
		// Check the progress
		if ( idx % 1000 == 0 ) {
			printf( "[STEP 3] Compressing the sequences is processing...[%lu/%lu]\n", idx, qscores.size() );
			fflush( stdout );
		}
	}
	// Terminate
	printf( "[STEP 3] Compressing the sequences is done!\n" );
	printf( "---------------------------------------------------------------------\n" );
	fflush( stdout );
}


int main( int argc, char **argv ) {
	char *filenameS = "/mnt/ssd0/semin/bancroft/data/sequences/fastq/HG002.fastq";

	// Read sequence file
	readerFASTQ( filenameS );

	// Make hash table to get the most frequenly occurred character
	qscoreHashTableMaker();

	// Compression
	uint64_t base = 5;
	for ( uint64_t iter = 0; iter < 1; iter ++ ) {
		uint64_t stride = 0;
		// Set the stride
		if ( iter == 0 ) stride = base;
		else stride = (base * 4) + 1;
		// Update the base
		base = stride;
		// Variable initialization
		qscoreSizeCmpr = 0;
		qscoreSizeRmnd = 0;
		// Compress
		double processStart = timeChecker();
		compressor_unit_ch( stride );
		double processFinish = timeChecker();
		double elapsedTime = processFinish - processStart;
		// Results
		printf( "QUALITY SCORE\n" );
		printf( "The Original File Size   : %0.4f MB\n", (double)qscoreSizeOrgn / 1024.00 / 1024.00 );
		printf( "---------------------------------------------------------------------\n" );
		printf( "COMPRESSION RESULT\n" );
		printf( "Stride                   : %lu\n", stride );
		printf( "The Compressed File Size : %0.4f MB\n", (double)(qscoreSizeCmpr + qscoreSizeRmnd) / 8.00 / 1024.00 / 1024.00);
		printf( "---------------------------------------------------------------------\n" );
		printf( "Previous Subsequence     : %0.4f\n", ((double)portion_1 / (double)portion_t) * 100.00 );
		printf( "Most Frequence Character : %0.4f\n", ((double)portion_2 / (double)portion_t) * 100.00 );
		printf( "Previous Character       : %0.4f\n", ((double)portion_3 / (double)portion_t) * 100.00 );
		printf( "Unmatched                : %0.4f\n", ((double)portion_4 / (double)portion_t) * 100.00 );
		printf( "Sanity Check             : %lu, %lu\n", ((portion_t * stride) + portion_5), qscoreSizeOrgn );
		printf( "Elapsed Time: %lf\n", elapsedTime );
		printf( "---------------------------------------------------------------------\n" );
	}

	return 0;
}
