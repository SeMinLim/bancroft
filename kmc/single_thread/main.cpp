#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <unordered_map>
using namespace std;


#define KMERLENGTH 32


unordered_map<string, int> reference;
vector<string> sequences;
size_t seqSizeOrg = 0;


static inline double timeChecker( void ) {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)(tv.tv_sec) + (double)(tv.tv_usec) / 1000000;
}

void fastaReader( char *filename ) {
	string seqLine;
	sequences.clear();
	sequences.shrink_to_fit();
	seqLine.clear();
	seqLine.shrink_to_fit();

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
		seqLine.clear();
		seqLine.shrink_to_fit();
	}

	f_data_sequences.close();
}

void fastaWriter( char *filename ) {
	ofstream f_data_result(filename);
	if ( !f_data_result.is_open() ) {
		printf( "File not found: %s\n", filename );
		exit(1);
	}

	for ( auto it = reference.begin(); it != reference.end(); it ++ ) {
		f_data_result << it->first << "\t" << it->second << "\n";
	}
	f_data_result << endl;

	f_data_result.close();
}

void kmc( void ) {
	for ( size_t i = 0; i < sequences.size(); i ++ ) {
		for ( size_t j = 0; j < sequences[i].size() - KMERLENGTH + 1; j ++ ) {
			if ( j % 1000000 == 0 ) printf( "%ld\n", j );
			string subseq = sequences[i].substr(j, 32);
			if ( reference.insert(make_pair(subseq, 1)).second == false ) {
				reference.at(subseq) += 1;
			}
		}
		printf( "Processing %ld sequence is finished!\n", i );
	}
}

void refOrganizer( void ) {
	for ( auto it = reference.begin(); it != reference.end(); ) {
		if ( it->second == 1 ) {
			reference.erase(it);
		} else it++;
	}
}


int main() {
	char *filenameS = "../../../data/references/hg19.fasta";
	char *filenameR = "../../../data/references/hg19RefBookARDA.txt";
	
	// Read sequence file
	fastaReader( filenameS );
	printf( "Read sequence file is finished!\n" );
	fflush( stdout );

	// Kmer counting
	double processStartKmc = timeChecker();
	kmc();
	double processFinishKmc = timeChecker();
	double elapsedTimeKmc = processFinishKmc - processStartKmc;
	printf( "KMC is finished!\n" );
	fflush( stdout );
	
	// Sorting and remove useless kmer
	double processStartOrg = timeChecker();
	refOrganizer();
	double processFinishOrg = timeChecker();
	double elapsedTimeOrg = processFinishOrg - processStartOrg;
	printf( "Sorting reference book is finished!\n" );
	fflush( stdout );

	// Write reference book
	fastaWriter( filenameR );
	printf( "Writing reference book is finished!\n" );
	fflush( stdout );

	printf( "--------------------------------------------\n" );
	printf( "KMC RESULT\n" );
	printf( "Reference Book [KMER]: %ld\n", reference.size() );
	printf( "Reference Book [Size]: %0.4f GB\n", (double)(reference.size() * 32) / 1024 / 1024 / 1024 );
	printf( "Elapsed Time: %lf\n", elapsedTimeKmc + elapsedTimeOrg );
	printf( "--------------------------------------------\n" );
	
	return 0;
}
