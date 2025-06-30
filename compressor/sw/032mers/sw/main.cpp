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


#define KMERLENGTH 32
#define ENCKMERBUFUNIT 32
#define ENCKMERBUFSIZE 1
#define BINARYRWUNIT 8
#define CHROMOSOMEUNIT 1
#define GROUPVARINT 1
#define FASTQ 1


string sequence;
vector<string> sequences;
map<uint64_t, uint32_t> reference;


uint64_t seqSizeOrg = 0;
uint64_t seqSizeCmpN = 0;
uint64_t seqSizeCmpI = 0;
uint64_t seqSizeCmpP = 0;
uint64_t seqSizeRmnd = 0;


// Reference: hg19From1
//uint64_t refSizeOrg = 2552813890;
//uint64_t refSizeUsd = 2552813890;
// Reference: B73From1
//uint64_t refSizeOrg = 916513317;
//uint64_t refSizeUsd = 916513317;
// Reference: GRCm39
uint64_t refSizeOrg = 2231482983;
uint64_t refSizeUsd = 2231482983;


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
void seqReaderFASTQ( char *filename ) {
	uint64_t counter = 0;
	string seqLine;
	// Read
	ifstream f_data_sequences(filename);
	while ( getline(f_data_sequences, seqLine) ) {
		if ( counter == 0 ) {
			counter ++;
		} else if ( counter == 1 ) {
			if ( CHROMOSOMEUNIT ) sequences.push_back(seqLine);
			else sequence += seqLine;
			seqSizeOrg += seqLine.size();
			counter ++;
		} else if ( counter == 2 ) {
			counter ++;
		} else if ( counter == 3 ) {
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
// Reference File Reader
void refReader( char *filename ) {
	ifstream f_data_reference(filename, ios::binary);
	for ( uint64_t i = 0; i < refSizeUsd; i ++ ) {
		uint64_t encKmer[ENCKMERBUFSIZE + 1] = {0, };
		// Read 32-mer and index
		for ( uint64_t j = 0; j < ENCKMERBUFSIZE + 1; j ++ ) {
			f_data_reference.read(reinterpret_cast<char *>(&encKmer[j]), BINARYRWUNIT);
		}
		// Insert 32-mer and index to Map
		if ( reference.insert(make_pair(encKmer[0], (uint32_t)encKmer[1])).second == false ) {
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
void encoder( string seqLine, uint64_t &encKmer ) {
	encKmer = 0;
	for ( uint64_t i = 0; i < ENCKMERBUFUNIT; i ++ ) {
		if ( seqLine[i] == 'A' ) {
			encKmer = ((uint64_t)0 << 2 * i) | encKmer;
		} else if ( seqLine[i] == 'C' ) {
			encKmer = ((uint64_t)1 << 2 * i) | encKmer;
		} else if ( seqLine[i] == 'G' ) {
			encKmer = ((uint64_t)2 << 2 * i) | encKmer;
		} else if ( seqLine[i] == 'T' ) {
			encKmer = ((uint64_t)3 << 2 * i) | encKmer;
		}
	}
}
// 2-bit Decoder
void decoder( const uint64_t encKmer, string &seqLine ) {
	for ( uint64_t j = 0; j < ENCKMERBUFUNIT; j ++ ) {
		uint64_t encCharT = encKmer << (ENCKMERBUFUNIT - 1 - j) * 2;
		uint64_t encCharF = encCharT >> (ENCKMERBUFUNIT - 1) * 2;
		if ( encCharF == 0 ) seqLine.push_back('A');
		else if ( encCharF == 1 ) seqLine.push_back('C');
		else if ( encCharF == 2 ) seqLine.push_back('G');
		else if ( encCharF == 3 ) seqLine.push_back('T');
	}
}
// Reverse Complement
void complementer( string reversed, string &complemented ) {
	for ( uint64_t i = 0; i < KMERLENGTH; i ++ ) {
		if ( reversed[i] == 'A' ) complemented.push_back('T');
		else if ( reversed[i] == 'C' ) complemented.push_back('G');
		else if ( reversed[i] == 'G' ) complemented.push_back('C');
		else if ( reversed[i] == 'T' ) complemented.push_back('A');
	}
}
// Compressor [CHROMOSOME UNIT]
void compressor_unit_ch( const uint64_t stride ) {
	uint32_t prevIndex = 0;
	uint32_t currIndex = 0;
	for ( uint64_t seqIdx = 0; seqIdx < sequences.size(); seqIdx ++ ) {
		uint64_t start = 0;
		while ( (double)start <= (double)sequences[seqIdx].size() - (double)KMERLENGTH ) {
			// Get subsequence
			string subseq = sequences[seqIdx].substr(start, KMERLENGTH);
			string subseqOrg = subseq;
			reverse(subseq.begin(), subseq.end());
			string subseqRev = subseq;
			// Get reverse complement
			string subseqCom;
			complementer(subseqRev, subseqCom);
			// Encode original subsequence first
			uint64_t encSubseqOrg = 0;
			encoder(subseqOrg, encSubseqOrg);
			// Check possible to compress
			if ( reference.find(encSubseqOrg) != reference.end() ) {
				// Possible to compress, then store the current index
				currIndex = reference.at(encSubseqOrg);
				// Compare the current index with the previous one
				if ( seqSizeCmpP != 0 ) {
					if ( (currIndex == prevIndex + KMERLENGTH) || (currIndex == prevIndex - KMERLENGTH) ) {
						seqSizeCmpI ++;
					}
				}
				// Update the parameters
				prevIndex = currIndex;
				seqSizeCmpP ++;
				start += KMERLENGTH;
			} else {
				// Encode reverse complement
				uint64_t encSubseqCom = 0;
				encoder(subseqCom, encSubseqCom);
				// Check possible to compress
				if ( reference.find(encSubseqCom) != reference.end() ) {
					// Possible to compress, then store the current index
					currIndex = reference.at(encSubseqCom);
					// Compare the current index with the previous one
					if ( seqSizeCmpP != 0 ) {
						if ( (currIndex == prevIndex + KMERLENGTH) || (currIndex == prevIndex - KMERLENGTH) ) {
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
		}
		// Handle remainder
		uint64_t remainder = sequences[seqIdx].size() - start;
		if ( remainder > 0 ) seqSizeRmnd += (remainder * 2) + 1;
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
	uint32_t prevIndex = 0;
	uint32_t currIndex = 0;
	uint64_t start = 0;
	while ( start <= sequence.size() - KMERLENGTH ) {
		// Get subsequence
		string subseq = sequence.substr(start, KMERLENGTH);
		string subseqOrg = subseq;
		reverse(subseq.begin(), subseq.end());
		string subseqRev = subseq;
		// Get reverse complement
		string subseqCom;
		complementer(subseqRev, subseqCom);
		// Encode original subsequence first
		uint64_t encSubseqOrg = 0;
		encoder(subseqOrg, encSubseqOrg);
		// Check possible to compress
		if ( reference.find(encSubseqOrg) != reference.end() ) {
			// Possible to compress, then store the current index
			currIndex = reference.at(encSubseqOrg);
			// Compare the current index to the previous one
			if ( seqSizeCmpP != 0 ) {
				if ( (currIndex == prevIndex + KMERLENGTH) || (currIndex == prevIndex - KMERLENGTH) ) {
					seqSizeCmpI ++;
				}
			}
			// Update the parameters
			prevIndex = currIndex;
			seqSizeCmpP ++;
			start += KMERLENGTH;
		} else {
			// Encode reverse complement
			uint64_t encSubseqCom = 0;
			encoder(subseqCom, encSubseqCom);
			// Check possible to compress
			if ( reference.find(encSubseqCom) != reference.end() ) {
				// Possible to compress, then store the current index
				currIndex = reference.at(encSubseqCom);
				// Compare the current index with the previous one
				if ( seqSizeCmpP != 0 ) {
					if ( (currIndex == prevIndex + KMERLENGTH) || (currIndex == prevIndex + KMERLENGTH) ) {
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
	char *filenameS = "/mnt/ephemeral/sequence/longread/fastq/mouse-rep2.fastq";
	char *filenameR = "/mnt/ephemeral/reference/GRCm39Reference032MersFrom1IndexIncluded.bin";

	// Read sequence file
	if ( FASTQ ) seqReaderFASTQ( filenameS );
	else seqReaderFASTA( filenameS );

	// Read reference file
	refReader( filenameR );

	// Compression
	for ( uint64_t stride = 1; stride < 64; stride = stride * 2 ) {
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

	return 0;
}
