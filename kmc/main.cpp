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
#include <map>
#include <unordered_map>
using namespace std;


#define KMERLENGTH 34


unordered_map<string, int> reference;
vector<string> sequences;
size_t seqSize = 0;
size_t refSize = 0;


static inline double timeChecker( void ) {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)(tv.tv_sec) + (double)(tv.tv_usec) / 1000000;
}

bool decendingOrder( pair<string, int>& x, pair<string, int>& y ) {
	return x.second > y.second;
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
			seqSize += seqLine.size();
		}
		seqLine.clear();
		seqLine.shrink_to_fit();
	}

	f_data_sequences.close();
	printf( "Reading Sequence File is Finished!\n" );
}

void fastaWriter( char *filename ) {
	// Sort first before writing the result
	vector<pair<string, int>> reference_vec(reference.begin(), reference.end());
	sort(reference_vec.begin(), reference_vec.end(), decendingOrder);
	printf( "Sorting Task is Finished!\n" );

	ofstream f_data_result(filename);
	if ( !f_data_result.is_open() ) {
		printf( "File not found: %s\n", filename );
		exit(1);
	}

	for ( size_t i = 0; i < reference_vec.size(); i ++ ) {
		f_data_result << reference_vec[i].first << "\t" << reference_vec[i].second << "\n";
	}
	f_data_result << endl;

	f_data_result.close();
	printf( "Writing Reference Book is Finished!\n" );
}

void kmc( void ) {
	// K-mer counting
	for ( size_t i = 0; i < sequences.size(); i ++ ) {
		for ( size_t j = 0; j < sequences[i].size() - KMERLENGTH + 1; j ++ ) {
			if ( j % 1000000 == 0 ) printf( "%ld\n", j );
			string subseq = sequences[i].substr(j, KMERLENGTH);
			if ( reference.insert(make_pair(subseq, 1)).second == false ) {
				reference.at(subseq) += 1;
			}
		}
		printf( "Processing %ld sequence is finished!\n", i );
	}
	printf( "KMC is Finished!\n" );

	// Erase the elements that have a value as 1
	for ( auto it = reference.begin(); it != reference.end(); ) {
		if ( it->second == 1 ) {
			it = reference.erase(it);
		} else it++;
	}
	refSize = reference.size();
	printf( "Erasing Task is Finished!\n" );
}


int main() {
	char *filenameS = "/dfs6/pub/seminl1/DNACompressor/hg19.fasta";
	char *filenameR = "/dfs6/pub/seminl1/DNACompressor/hg19RefBookARDA.txt";
	
	// Read sequence file
	fastaReader( filenameS );

	// Kmer counting
	double processStartKmc = timeChecker();
	kmc();
	double processFinishKmc = timeChecker();
	double elapsedTimeKmc = processFinishKmc - processStartKmc;
	
	// Sorting, erasing, and writing
	fastaWriter( filenameR );

	printf( "--------------------------------------------\n" );
	printf( "KMC RESULT\n" );
	printf( "KMER [Total]: %ld\n", seqSize );
	printf( "KMER [Count]: %ld\n", refSize );
	printf( "KMER [Percentage]: %0.8f\n", (double)refSize / seqSize * 100 );
	printf( "Reference Book [Size]: %0.4f GB\n", (double)(refSize * 32) / 1024 / 1024 / 1024 );
	printf( "Elapsed Time: %lf\n", elapsedTimeKmc );
	printf( "--------------------------------------------\n" );
	
	
	return 0;
}
