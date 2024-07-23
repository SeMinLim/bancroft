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
#define REFINDEX 32
#define TESTSEQ 5


vector<string> sequences;
unordered_map<string, uint32_t> reference;


uint64_t seqSizeOrg = 0;
uint64_t seqSizeCmp = 0;
uint32_t refIdx = 0;
//size_t usedReference = 4194304;


static inline double timeChecker( void ) {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)(tv.tv_sec) + (double)(tv.tv_usec) / 1000000;
}

void fastaReader( char *filename ) {
	string seqLine;

	ifstream f_data_sequences(filename);
	if ( !f_data_sequences.is_open() ) {
		printf( "File not found: %s\n", filename );
		exit(1);
	}

	while ( getline(f_data_sequences, seqLine) ) {
		if ( seqLine[0] != '>' ) {
			sequences.push_back(seqLine);
			seqSizeOrg += seqLine.size();
		}
	}

	f_data_sequences.close();
}

void refBookReader_1( char *filename ) {
	reference.clear();
	string refLine;

	ifstream f_data_reference(filename);
	if ( !f_data_reference.is_open() ) {
		printf( "File not found: %s\n", filename );
		exit(1);
	}

	while ( getline(f_data_reference, refLine) ) {
		string kmer;
		kmer.reserve(KMERLENGTH);

		for ( size_t i = 0; i < KMERLENGTH; i ++ ) {
			kmer.push_back(refLine[i]);
		}
		if ( reference.insert(make_pair(kmer, refIdx)).second == false ) {
			printf( "There's an issue on reference code book\n" );
			exit(1);
		}
		refIdx++;
	}

	f_data_reference.close();
}

bool subseqFinder( string subseq ) {
	if ( reference.find(subseq) != reference.end() ) {
		seqSizeCmp++;
		return true;
	} else {
		return false;
	}
}

void compressor( void ) {
	char *filenameRef_1 = "/mnt/smartssd0/semin/hg19hg38RefBook256Mers_1.txt";
	char *filenameRef_2 = "/mnt/smartssd0/semin/hg19hg38RefBook256Mers_2.txt";

	size_t sp = 0;
	bool success = false;	
	while ( sp <= sequences[TESTSEQ].size() - KMERLENGTH ) {
		string subseq = sequences[TESTSEQ].substr(sp, KMERLENGTH);
		
		// Reference 1
		refIdx = 0;
		refBookReader_1(filenameRef_1);
		success = subseqFinder(subseq);
		if ( success == true ) {
			sp += KMERLENGTH;
		} else {
			// Reference 2
			string refLine;
			ifstream f_data_reference(filenameRef_2);
			while ( true ) {
				reference.clear();
				while ( getline(f_data_reference, refLine) ) {
					string kmer;
					kmer.reserve(KMERLENGTH);
	
					for ( size_t i = 0; i < KMERLENGTH; i ++ ) {
						kmer.push_back(refLine[i]);
					}
					if ( reference.insert(make_pair(kmer, refIdx)).second == false ) {
						printf( "There's an issue on reference code book\n" );
						exit(1);
					}
					refIdx++;
					if ( refIdx % 21923549 == 0  ) break;
				}
				printf( "%d\n", refIdx );

				success = subseqFinder(subseq);
				if ( success == true ) {
					sp += KMERLENGTH;
					break;
				} else {
					if ( refIdx == 2836860451 ) {
						sp += 1;
						break;
					}
				}
			}

			f_data_reference.close();
		}
	}
}


int main( void ) {
	char *filenameS = "../../data/sequences/hg16.fasta";

	// Read sequence file
	fastaReader( filenameS );

	// Compression
	double processStart = timeChecker();
	compressor();
	double processFinish = timeChecker();
	double elapsedTime = processFinish - processStart;

	printf( "--------------------------------------------\n" );
	printf( "REFERENCE\n" );
	printf( "Reference Book [#KMER]: %ld\n", reference.size() );
	printf( "Reference Book [Size]: %0.4f MB\n", ((double)reference.size() * KMERLENGTH) / 1024 / 1024 );
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
