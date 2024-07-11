#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <pthread.h>
#include <time.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;


#define NUMTHREAD 16
#define KMERLENGTH 10
#define TESTSEQ 24


typedef struct RefBook {
	string kmer;
	int cnt;
}RefBook;
typedef struct PthreadArg {
	size_t seqIdx;
	size_t thrIdx;
	size_t subSeq;
}PthreadArg;


vector<string> sequences;
vector<RefBook> reference;
uint64_t seqSizeOrg = 0;
uint64_t seqSizeCmp = 0;
//size_t usedReference = 131072;

pthread_mutex_t mutex;


static inline double timeChecker( void ) {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)(tv.tv_sec) + (double)(tv.tv_usec) / 1000000;
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

		kmer.reserve(KMERLENGTH);		
		for ( size_t i = 0; i < KMERLENGTH; i ++ ) {
			kmer.push_back(refLine[i]);
		}
		
		refLine.erase(0, KMERLENGTH+1);
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

void *compressor(void *pthreadArg) {
	PthreadArg *arg = (PthreadArg *)pthreadArg;
	size_t sp = arg->thrIdx * arg->subSeq;
	size_t fp = 0;
	if ( arg->thrIdx + 1 == NUMTHREAD ) fp = sequences[arg->seqIdx].size();
	else fp = (arg->thrIdx + 1) * arg->subSeq;
	
	size_t flag = 0;
	while ( sp <= fp - KMERLENGTH) {
		string subseq = sequences[arg->seqIdx].substr(sp, KMERLENGTH);
		for ( size_t j = 0; j < reference.size(); j ++ ) {
			if ( subseq.compare(reference[j].kmer) == 0 ) {
				pthread_mutex_lock(&mutex);
				seqSizeCmp++;
				pthread_mutex_unlock(&mutex);
				flag = 1;
				break;
			} else flag = 0;
		}
		
		if ( flag == 1 ) sp += KMERLENGTH;
		else sp += 1;
	}
	return 0;
}


int main() {
	char *filenameS = "../../data/sequences/hg38.fasta";
	char *filenameR = "../../data/references/hg19RefBook10Mers.txt";

	// Read sequence file
	fastaReader( filenameS );

	// Read reference file
	refBookReader( filenameR );

	// Compression
	pthread_t pthread[NUMTHREAD];
	pthread_mutex_init(&mutex, NULL);
	PthreadArg pthreadArg[NUMTHREAD];

	double processStart = timeChecker();
	size_t subSeq = sequences[TESTSEQ].size() / NUMTHREAD;
	// Thread Create
	for ( size_t i = 0; i < NUMTHREAD; i ++ ) {
		pthreadArg[i].seqIdx = TESTSEQ;
		pthreadArg[i].thrIdx = i;
		pthreadArg[i].subSeq = subSeq;

		pthread_create(&pthread[i], NULL, compressor, (void *)&pthreadArg[i]);
	}
	// Thread Join
	for ( size_t j = 0; j < NUMTHREAD; j ++ ) {
		pthread_join(pthread[j], NULL);
	}
	// Mutex Destroy
	pthread_mutex_destroy(&mutex);
	double processFinish = timeChecker();
	double elapsedTime = processFinish - processStart;

	printf( "--------------------------------------------\n" );
	printf( "REFERENCE\n" );
	printf( "Original Reference Book [#KMER]: %ld\n", reference.size() );
	printf( "Original Reference Book [Size]: %0.4f MB\n", ((double)reference.size() * KMERLENGTH) / 1024 / 1024 );
	printf( "Used Reference Book Size: %0.4f MB\n", ((double)reference.size() * KMERLENGTH) / 1024 / 1024 );
	printf( "--------------------------------------------\n" );
	printf( "SEQUENCE\n" );
	printf( "Original File Size: %0.4f MB\n", (double)sequences[TESTSEQ].size() / 1024 / 1024 );
	printf( "Number of Base Pairs [Original]: %ld\n", sequences[TESTSEQ].size() );
	printf( "--------------------------------------------\n" );
	printf( "COMPRESSION RESULT\n" );
	printf( "Compressed File Size: %0.4f MB\n", (double)((sequences[TESTSEQ].size() - (seqSizeCmp * KMERLENGTH)) + (seqSizeCmp * 3)) / 1024 / 1024 );		
	printf( "Number of Base Pairs [Compressed]: %ld\n", seqSizeCmp * KMERLENGTH );
	printf( "Elapsed Time: %lf\n", elapsedTime );

	return 0;
}
