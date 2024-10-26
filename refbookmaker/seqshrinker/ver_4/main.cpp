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


#define KMERLENGTH 256
#define SHRINKLENGTH 1073741824
#define ENCKMERBUFUNIT 32
#define ENCKMERBUFSIZE 8
#define BINARYRWUNIT 8


string sequences;
vector<uint32_t> index_1;
vector<vector<uint32_t>> index_2;
vector<uint32_t> index_size;
map<pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, uint64_t>>>>>>>, 
    uint32_t> 
    reference_kmer;
map<uint32_t, 
    pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, uint64_t>>>>>>>> 
    reference_index;


uint64_t seqSizeOrg = 0;
uint64_t refSizeOrg = 2849207900;
uint64_t refSizeUsd = 536870912;


//// Required Functions
// Time Checker
static inline double timeChecker( void ) {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)(tv.tv_sec) + (double)(tv.tv_usec) / 1000000;
}
// Descending Order Sorter
bool descendingOrder( vector<vector<uint32_t>> &x, vector<vector<uint32_t>> &y ) {
	return x.size() > y.size();
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
// Assembled Sequence File Reader
void seqReader( char *filename ) {
	string seqLine;

	ifstream f_data_sequences(filename);
	while ( getline(f_data_sequences, seqLine) ) {
		if ( seqLine[0] != '>' ) {
			sequences += seqLine;
			seqSizeOrg += seqLine.size();
		}
	}

	f_data_sequences.close();
	printf( "[STEP 1] Reading sequence fasta file is done!\n" );
	fflush( stdout );
}
// Pre-Made Reference File Reader
void refReader( char *filename ) {
	ifstream f_data_reference(filename, ios::binary);
	for ( uint64_t i = 0; i < refSizeUsd; i ++ ) {
		uint64_t encKmer[ENCKMERBUFSIZE] = {0, };
		// Read
		for ( uint64_t j = 0; j < ENCKMERBUFSIZE; j ++ ) {
			f_data_reference.read(reinterpret_cast<char *>(&encKmer[j]), BINARYRWUNIT);
		}
		// Insert 256-Mers to Map_1
		if ( reference_kmer.insert(make_pair(make_pair(encKmer[0], make_pair(encKmer[1], make_pair(encKmer[2], 
				           make_pair(encKmer[3], make_pair(encKmer[4], make_pair(encKmer[5], 
				           make_pair(encKmer[6], encKmer[7]))))))), i)).second == false ) {
			printf( "There's a problem on reference code book...\n" );
			fflush( stdout );
			exit(1);
		}
		// Insert 256-Mers to Map_2
		if ( reference_index.insert(make_pair(i, make_pair(encKmer[0], make_pair(encKmer[1], make_pair(encKmer[2], 
				            make_pair(encKmer[3], make_pair(encKmer[4], make_pair(encKmer[5], 
				            make_pair(encKmer[6], encKmer[7]))))))))).second == false ) {
			printf( "There's a problem on reference code book...\n" );
			fflush( stdout );
			exit(1);
		}
		// Check the progress
		if ( i % 1000000 == 0 ) {
			printf( "[STEP 2] Reading pre-made reference file is processing...[%lu]\n", i );
			fflush( stdout );
		}
	}

	f_data_reference.close();
	printf( "[STEP 2] Reading pre-made reference file is done!\n" );
	fflush( stdout );
}
// Sequence Shrinker
void seqShrinker( char *filename ) {
	// Get the index of sequence for each reference slot first
	uint32_t start = 0;
	while ( start <= sequences.size() - KMERLENGTH ) {
		string subseq = sequences.substr(start, KMERLENGTH);

		// Encode subsequence first
		uint64_t encSubseq[ENCKMERBUFSIZE] = {0, };
		encoder(subseq, encSubseq);

		// Store the index to vector if there is the same subsequence in the pre-made reference book
		if ( reference.find(make_pair(encSubseq[0], make_pair(encSubseq[1], make_pair(encSubseq[2], 
				    make_pair(encSubseq[3], make_pair(encSubseq[4], make_pair(encSubseq[5], 
				    make_pair(encSubseq[6], encSubseq[7])))))))) != reference.end() ) {
			index_1.push_back(start);
		}
		
		// Move forward to the next index
		start += 1;

		// Check the progress
		if ( start % 1000000 == 0 ) {
			printf( "[STEP 3] Getting the index of the sequence is processing...[%lu]\n", start );
			fflush( stdout );
		}
	}
	printf( "[STEP 3] Getting the index of the sequence is done!\n" );
	fflush( stdout );
	
	// Do sorting index_1
	sort(index_1.begin(), index_1.end());
	printf( "[STEP 4] Sorting the index is done!\n" );
	fflush( stdout );

	// Make the groups of index consisting of the sequential index
	uint32_t group = 0;
	for ( uint32_t i = 0 i < index_1.size(); i ++ ) {
		// Initialization
		if ( i == 0 ) {
			index_2[group].push_back(index_1[i]);
		} else {
			if ( index_1[i] = index_1[i-1] + 1 ) {
				index_2[group].push_back(index_1[i]);
			} else {
				index_2[++group].push_back(index_1[i]);
			}
		}
	}
	sort(index_2.begin(), index_2.end(), descendingOrder);
	printf( "[STEP 5] Grouping the index is done!\n" );
	fflush( stdout );
	
	// Check each group that has the same subsequence with other group
	for ( uint32_t i = index_2.size() - 1; i = 0; ) {
		uint32_t flag = 0;
		for ( uint32_t j = 0; j < index_2[i].size(); j ++ ) {
			for ( uint32_t k = 0; k < index_2[i-1].size(); k ++ ) {
				if ( reference_index.at(index_2[i][j]) == reference_index.at(index_2[i-1][k]) ) {
					index_2.erase(index_2.begin() + i);
					flag = 1;
					break;
				}
			}
			if ( flag == 1 ) break;
		}
		if (flag == 0) i --;
	}
	printf( "[STEP 6] Erasing a smaller group is done!\n" );
	fflush( stdout );

	// Check the length of new sequence
	uint32_t size = 0;
	for ( uint32_t i = 0; i < index_2.size(); i ++ ) {
		size = 256 + index_2[i].size() - 1;
	}
	printf( "size: %u\n", size );
	fflush( stdout );
/*
	// Write the reference as binary file
	ofstream f_data_result(filename, ios::binary);
	for ( uint64_t i = 0; i < reference_vector.size(); i ++ ) {
		uint64_t encKmer[ENCKMERBUFSIZE] = {0, };
		encKmer[0] = reference_vector[i].first.first;
		encKmer[1] = reference_vector[i].first.second.first;
		encKmer[2] = reference_vector[i].first.second.second.first;
		encKmer[3] = reference_vector[i].first.second.second.second.first;
		encKmer[4] = reference_vector[i].first.second.second.second.second.first;
		encKmer[5] = reference_vector[i].first.second.second.second.second.second.first;
		encKmer[6] = reference_vector[i].first.second.second.second.second.second.second.first;
		encKmer[7] = reference_vector[i].first.second.second.second.second.second.second.second;
		for ( uint64_t j = 0; j < ENCKMERBUFSIZE; j ++ ) {
			f_data_result.write(reinterpret_cast<char *>(&encKmer[j]), BINARYRWUNIT);
		}
	}
	f_data_result.close();
	printf( "[STEP 4] Writing the k-mers as binary is done\n" );
	fflush( stdout );
*/
}


int main() {
	char *filenameIn = "/mnt/ephemeral/hg19.fasta";
	char *filenameOut = "/home/jovyan/hg19ReferenceShrinked256Mers256M.bin";
	
	// Read sequence file
	seqReader( filenameIn );

	// Kmer counting
	kmc( filenameOut);
	
	return 0;
}
