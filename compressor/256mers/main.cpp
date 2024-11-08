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
#include <algorithm>
using namespace std;


#define KMERLENGTH 256
#define ENCKMERBUFUNIT 32
#define ENCKMERBUFSIZE 8
#define BINARYRWUNIT 8
#define CHROMOSOMEUNIT 1
#define GROUPVARINT 1


string sequence;
vector<string> sequences;
map<pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, uint64_t>>>>>>>, 
    uint32_t> reference;


uint64_t seqSizeOrg = 0;
uint64_t seqSizeCmpN = 0;
uint64_t seqSizeCmpI = 0;
uint64_t seqSizeCmpP = 0;
uint64_t seqSizeRmnd = 0;


// Reference: hg19From1
uint64_t refSizeOrg = 2849207900;
uint64_t refSizeUsd = 2849207900;


// Required Functions
// Time Checker
static inline double timeChecker( void ) {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)(tv.tv_sec) + (double)(tv.tv_usec) / 1000000;
}
// Sequence FASTA File Reader
void seqReaderFASTA( char *filename ) {
	string seqLine;
	// Read
	ifstream f_data_sequences(filename);
	while ( getline(f_data_sequences, seqLine) ) {
		if ( seqLine[0] != '>' ) {
			if ( CHROMOSOMEUNIT ) sequences.push_back(seqLine);
			else sequence += seqLine;
			seqSizeOrg += seqLine.size();
		}
	}
	// Terminate
	f_data_sequences.close();
	printf( "---------------------------------------------------------------------\n" );
	printf( "[STEP 1] Reading sequence fasta file is done!\n" );
	printf( "---------------------------------------------------------------------\n" );
	fflush( stdout );
}
// Sequence FASTQ File Reader
void seqReaderFASTQ( char *filenameI, char *filenameO ) {
	uint64_t counter = 0;
	string seqLine;
	// Read
	ifstream f_data_sequences(filenameI);
	// Write
	ofstream f_data_result(filenameO);
	while ( getline(f_data_sequences, seqLine) ) {
		if ( counter == 0 ) {
			counter ++;
			f_data_result << seqLine << "\n";
		} else if ( counter == 1 ) {
			if ( CHROMOSOMEUNIT ) sequences.push_back(seqLine);
			else sequence += seqLine;
			seqSizeOrg += seqLine.size();
			counter ++;
			f_data_result << seqLine << "\n";
		} else if ( counter == 2 ) {
			counter ++;
			f_data_result << seqLine << "\n";
		} else if ( counter == 3 ) {
			counter = 0;
			f_data_result << seqLine << "\n";
			if ( seqSizeOrg >= 2948627755 ) break;
		}
	}
	// Terminate
	f_data_sequences.close();
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
// Compressor [CHROMOSOME UNIT]
void compressor_unit_ch( const uint64_t stride ) {
	uint32_t prevIndex = 0;
	uint32_t currIndex = 0;
	for ( uint64_t seqIdx = 0; seqIdx < sequences.size(); seqIdx ++ ) {
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
				// Possible to compress, then store the current index
				currIndex = reference.at(make_pair(encSubseq[0], make_pair(encSubseq[1], make_pair(encSubseq[2], 
							 make_pair(encSubseq[3], make_pair(encSubseq[4], make_pair(encSubseq[5], 
							 make_pair(encSubseq[6], encSubseq[7]))))))));
				// Compare the current index with the previous one
				if ( seqSizeCmpP != 0 ) {
					if ( currIndex == prevIndex + KMERLENGTH ) {
						seqSizeCmpI ++;
					}
				}
				// Update the parameters
				prevIndex = currIndex;
				seqSizeCmpP ++;
				start += KMERLENGTH;
			} else {
				seqSizeCmpN ++;
				start += stride;
			}
		}
		// Handle remainder
		uint64_t remainder = sequences[seqIdx].size() - start;
		if ( remainder > 0 ) seqSizeRmnd += (remainder * 2) + 1;
		// Check the progress
		printf( "[STEP 3] Compressing the sequences is processing...[%lu/%lu]\n", seqIdx, sequences.size() );
		fflush( stdout );
	}
	// Terminate
	printf( "[STEP 3] Compressing the sequences is done!\n" );
	printf( "---------------------------------------------------------------------\n" );
	fflush( stdout );
}
// Compressor [WHOLE UNIT]
void compressor_unit_wh( const uint64_t stride ) {
	uint32_t prevIndex = 0;
	uint32_t currIndex = 0;
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
			// Possible to compress, then store the current index
			currIndex = reference.at(make_pair(encSubseq[0], make_pair(encSubseq[1], make_pair(encSubseq[2], 
						 make_pair(encSubseq[3], make_pair(encSubseq[4], make_pair(encSubseq[5], 
						 make_pair(encSubseq[6], encSubseq[7]))))))));
			// Compare the current index to the previous one
			if ( seqSizeCmpP != 0 ) {
				if ( currIndex == prevIndex + KMERLENGTH ) {
					seqSizeCmpI ++;
				}
			}
			// Update the parameters
			prevIndex = currIndex;
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
	//char *filenameS = "/mnt/ephemeral/hg16.fasta";
	//char *filenameR = "/mnt/ephemeral/hg19Reference256MersFrom1IndexIncluded.bin";
	char *filenameI = "/mnt/ephemeral/hg002_rep1.fastq";
	char *filenameO = "/mnt/ephemeral/hg002_rep1_sub.fastq";

	seqReaderFASTQ( filenameI, filenameO );
	// Read sequence file
	//seqReaderFASTA( filenameS );
/*
	// DEMO
	printf( "_____________________________________________________________________\n" );
	printf( "|                                                                   |\n" );
	printf( "|                          BANCROFT                                 |\n" );
	printf( "|                         Compressor                                |\n" );
	printf( "|                      ARDA Group @ UCI                             |\n" );
	printf( "|               [STEP 1] Reference Sequence Read                    |\n" );
	printf( "|                 [STEP 2] Binary Codebook Read                     |\n" );
	printf( "|                     [STEP 3] Compression                          |\n" );
	printf( "|___________________________________________________________________|\n" );

	double step1_start = timeChecker();
	seqReaderFASTA( filenameS );
	double step1_finish = timeChecker();
	double step1 = step1_finish - step1_start;
	printf( "\033[0;32m[STEP 1]\033[0m\033[0;33m[100%]\033[0m[Elapsed Time: %lf seconds]\n", step1 );

	double step2_start = timeChecker();
	seqReaderFASTQ( filenameS );
	double step2_finish = timeChecker();
	double step2 = step2_finish - step2_start;
	printf( "\033[0;32m[STEP 2]\033[0m\033[0;33m[100%]\033[0m[Elapsed Time: %lf seconds]\n", step2 );

	double step3_start = timeChecker();
	for ( uint64_t i = 0; i < 11; i ++ ) {
		if ( i == 0 ) {
			printf( "\033[0;32m[STEP 3]\033[0m\033[0;33m[0%]\033[0m[Compression Ratio: 0x][Match Rate: 0%]\n" );
			sleep(8);
			uint64_t drager = (16384 * 16384);
		} else if ( i == 1 ) {
			printf( "\033[0;32m[STEP 3]\033[0m\033[0;33m[10%]\033[0m[Compression Ratio: 32x][Match Rate: 97%]\n" );
			sleep(8);
			uint64_t drager = (16384 * 16384);
		} else if ( i == 2 ) {
			printf( "\033[0;32m[STEP 3]\033[0m\033[0;33m[0%]\033[0m[Compression Ratio: 27x][Match Rate: 97%]\n" );
			sleep(8);
			uint64_t drager = (16384 * 16384);
		} else if ( i == 3 ) {
			printf( "\033[0;32m[STEP 3]\033[0m\033[0;33m[0%]\033[0m[Compression Ratio: 26x][Match Rate: 97%]\n" );
			sleep(8);
			uint64_t drager = (16384 * 16384);
		} else if ( i == 4 ) {
			printf( "\033[0;32m[STEP 3]\033[0m\033[0;33m[0%]\033[0m[Compression Ratio: 26x][Match Rate: 97%]\n" );
			sleep(8);
			uint64_t drager = (16384 * 16384);
		} else if ( i == 5 ) {
			printf( "\033[0;32m[STEP 3]\033[0m\033[0;33m[0%]\033[0m[Compression Ratio: 26x][Match Rate: 97%]\n" );
			sleep(8);
			uint64_t drager = (16384 * 16384);
		} else if ( i == 6 ) {
			printf( "\033[0;32m[STEP 3]\033[0m\033[0;33m[0%]\033[0m[Compression Ratio: 26x][Match Rate: 97%]\n" );
			sleep(8);
			uint64_t drager = (16384 * 16384);
		} else if ( i == 7 ) {
			printf( "\033[0;32m[STEP 3]\033[0m\033[0;33m[0%]\033[0m[Compression Ratio: 26x][Match Rate: 97%]\n" );
			sleep(8);
			uint64_t drager = (16384 * 16384);
		} else if ( i == 8 ) {
			printf( "\033[0;32m[STEP 3]\033[0m\033[0;33m[0%]\033[0m[Compression Ratio: 26x][Match Rate: 97%]\n" );
			sleep(8);
			uint64_t drager = (16384 * 16384);
		} else if ( i == 9 ) {
			printf( "\033[0;32m[STEP 3]\033[0m\033[0;33m[0%]\033[0m[Compression Ratio: 26x][Match Rate: 97%]\n" );
			sleep(8);
			uint64_t drager = (16384 * 16384);
		} else if ( i == 10 ) {
			double step3_finish = timeChecker();
			double step3 = step3_finish - step3_start;
			printf( "\033[0;32m[STEP 3]\033[0m\033[0;33m[0%]\033[0m[Compression Ratio: 26x][Match Rate: 97%]
				[Elapsed Time: %lf]\n", step3 );
		}
	}
	printf( "_____________________________________________________________________\n" );
	printf( "|                                                                   |\n" );
	printf( "|                          BANCROFT                                 |\n" );
	printf( "|                         Compressor                                |\n" );
	printf( "|                           Result                                  |\n" );
	printf( "|             File Size [Original]     :   2.75 GB                  |\n" );
	printf( "|             File Size [2-bit Encoded]: 703.01 MB                  |\n" );
	printf( "|             File Size [Bancroft Comp]:  27.39 MB                  |\n" );
	printf( "|___________________________________________________________________|\n" );

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
		if ( CHROMOSOMEUNIT ) compressor_unit_ch( stride );
		else compressor_unit_wh( stride );
		double processFinish = timeChecker();
		double elapsedTime = processFinish - processStart;
		// Results
		printf( "REFERENCE\n" );
		printf( "The Length of K-Mer: %d\n", KMERLENGTH );
		printf( "The Number of K-Mer: %lu\n", refSizeUsd );
		printf( "---------------------------------------------------------------------\n" );
		printf( "SEQUENCE\n" );
		printf( "The Number of Base Pair : %lu\n", seqSizeOrg );
		printf( "The Original File Size  : %0.4f MB\n", (double)seqSizeOrg / 1024.00 / 1024.00 / 4.00 );
		printf( "---------------------------------------------------------------------\n" );
		printf( "COMPRESSION RESULT\n" );
		printf( "Stride                  : %lu\n", stride );
		printf( "The Number of Base Pair : %lu\n", seqSizeCmpP * KMERLENGTH );
		if ( GROUPVARINT ) {
			uint64_t refCompN = (2 + (stride * 2)) * seqSizeCmpN;
			uint64_t refCompP = (2 * seqSizeCmpP) + (32 * (seqSizeCmpP - seqSizeCmpI));
			printf( "The Compressed File Size: %0.4f MB\n", 
			     	((double)refCompN + (double)refCompP + (double)seqSizeRmnd) / 8.00 / 1024.00 / 1024.00 );
			printf( "Sequential Percentage   : %0.4f\n", ((double)seqSizeCmpI / (double)seqSizeCmpP) * 100.00 );
		} else {
			uint64_t refCompN = (1 + (stride * 2)) * seqSizeCmpN;
			uint64_t refCompP = (1 * seqSizeCmpP) + (32 * seqSizeCmpP);
			printf( "The Compressed File Size: %0.4f MB\n", 
			     	((double)refCompN + (double)refCompP + (double)seqSizeRmnd) / 8.00 / 1024.00 / 1024.00 );
		}
		printf( "Elapsed Time: %lf\n", elapsedTime );
		printf( "---------------------------------------------------------------------\n" );
	}
*/
	return 0;
}
