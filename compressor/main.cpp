#include <sys/resource.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;


typedef struct RefBook {
	string kmer;
	int cnt;
}RefBook;


vector<string> sequences;
vector<RefBook> reference;
uint64_t seqSizeOrg = 0;
uint64_t seqSizeCmp = 0;


static inline double timeCheckerCPU(void) {
        struct rusage ru;
        getrusage(RUSAGE_SELF, &ru);
        return (double)ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000;
}

void fastaReader( char *filename ) {
	sequences.clear();
	sequences.shrink_to_fit();
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

void refBookReader( char *filename ) {
	reference.clear();
	reference.shrink_to_fit();
	string refLine;

	ifstream f_data_reference(filename);
	if ( !f_data_reference.is_open() ) {
		printf( "File not found: %s\n", filename );
		exit(1);
	}

	while ( getline(f_data_reference, refLine) ) {
		RefBook element;
		string kmer;
		int cnt;

		kmer.reserve(32);		
		for ( size_t i = 0; i < 32; i ++ ) {
			kmer.push_back(refLine[i]);
		}
		
		refLine.erase(0, 33);
		refLine.shrink_to_fit();
		cnt = atoi(refLine.c_str());

		element.kmer = kmer;
		element.cnt = cnt;

		reference.push_back(element);

		refLine.clear();
		refLine.shrink_to_fit();
	}

	f_data_reference.close();
}

void compressor(void) {
	for ( size_t i = 0; i < 1; i ++ ) {
		//printf( "%ldth sequence size: %ld\n", i, sequences[i].size() );
		size_t sp = 0;
		while ( sp + 32 < sequences[i].size() ) {
			//printf( "status: %ld\n", sp );
			string subseq = sequences[i].substr(sp, 32);
			for ( size_t j = 0; j < 131072; j ++ ) {
				if ( subseq.compare(reference[j].kmer) == 0 ) {
					seqSizeCmp++;
					break;
				}
			}
			sp += 32;
		}
		//printf( "%ldth sequence done!\n", i );
	}
}


int main() {
	char *filenameS = "../data/references/hg38.fasta";
	char *filenameR = "../data/references/hg19refbook.txt";

	printf( "Reading sequence file is started!\n" );
	fflush( stdout );
	fastaReader( filenameS );
	printf( "Reading sequence file is finished!\n" );
	fflush( stdout );

	printf( "Reading reference file is started!\n" );
	fflush( stdout );
	refBookReader( filenameR );
	printf( "Reading reference file is finished!\n" );
	fflush( stdout );

	printf( "Compression is started!\n" );
	fflush( stdout );
	double processStart = timeCheckerCPU();
	compressor();
	double processFinish = timeCheckerCPU();
	double elapsedTime = processFinish - processStart;
	printf( "Compression is finished!\n" );
	fflush( stdout );

	printf( "Reference book size: %0.4f GB\n", ((double)reference.size() * 32) / 1024 / 1024 / 1024 );
	printf( "Original file size: %0.4f MB\n", (double)sequences[0].size() / 1024 / 1024 );
	printf( "Compressed size: %0.4f MB\n", (double)((sequences[0].size() - (seqSizeCmp * 32)) + seqSizeCmp) / 1024 / 1024 );
	printf( "# of compressed kmer: %ld\n", seqSizeCmp );
	printf( "Elapsed Time: %f\n", elapsedTime );

	return 0;
}
