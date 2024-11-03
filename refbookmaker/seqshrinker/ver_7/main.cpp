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
#define ENCKMERBUFUNIT 32
#define ENCKMERBUFSIZE 8
#define BINARYRWUNIT 8
#define BLOCKLENGTH 1073741824


// Sequence
string sequence;
// Occurrence-Ordered Index
map<uint32_t, uint32_t> reference_index;
// <Index, K-Mer>
map<uint32_t, 
    pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, uint64_t>>>>>>>> 
    reference_kmer;
// <Index, 0>
map<uint32_t, uint32_t> 
    reference_final_index;
// <Index, K-Mer>
map<uint32_t, 
    pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, uint64_t>>>>>>>> 
    reference_final_kmer;


uint64_t seqSizeOrg = 0;
uint64_t seqSizeRdc = 0;
uint64_t seqSizeAln = 0;
uint64_t refSizeOrg = 2849207900;
uint64_t refSizeUsd = 2849207900;


// Required Functions
// Time Checker
static inline double timeChecker( void ) {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)(tv.tv_sec) + (double)(tv.tv_usec) / 1000000;
}
// Sequence File Reader
void seqReader( char *filename ) {
	string seqLine;
	// Read
	ifstream f_data_sequences(filename);
	while ( getline(f_data_sequences, seqLine) ) {
		if ( seqLine[0] != '>' ) {
			sequence += seqLine;
			seqSizeOrg += seqLine.size();
		}
	}
	// Terminate
	f_data_sequences.close();
	printf( "---------------------------------------------------------------------\n" );
	printf( "[STEP 1] The Length of Sequence: %lu\n", sequence.size() );
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
		// [1st] Insert index to the occurrence-ordered map
		if ( reference_index.insert(make_pair((uint32_t)i, encKmer[8])).second == false ) {
			printf( "There's a problem on reference code book...\n" );
			fflush( stdout );
			exit(1);
		}
		// [2nd] Insert index and k-mer to map<index, k-mer>
		if ( reference_kmer.insert(make_pair((uint32_t)encKmer[8], 
				      make_pair(encKmer[0], make_pair(encKmer[1], 
				      make_pair(encKmer[2], make_pair(encKmer[3], 
				      make_pair(encKmer[4], make_pair(encKmer[5], 
				      make_pair(encKmer[6], encKmer[7]))))))))).second == false ) {
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
// Variable Taker
void varTaker( uint64_t *encKmer, uint32_t index ) {
	encKmer[0] = reference_kmer[index].first;
	encKmer[1] = reference_kmer[index].second.first;
	encKmer[2] = reference_kmer[index].second.second.first;
	encKmer[3] = reference_kmer[index].second.second.second.first;
	encKmer[4] = reference_kmer[index].second.second.second.second.first;
	encKmer[5] = reference_kmer[index].second.second.second.second.second.first;
	encKmer[6] = reference_kmer[index].second.second.second.second.second.second.first;
	encKmer[7] = reference_kmer[index].second.second.second.second.second.second.second;
}
// Sequence Shrinker by De Bruijn
void seqShrinker( char *filename ) {
	// Align the indices first
	while ( seqSizeAln < BLOCKLENGTH ) {
		uint32_t start = 0;
		//// Phase 1
		// Pick the most frequently occurred k-mer as a start point
		for ( auto iter = reference_index.begin(); iter != reference_index.end(); ++ iter ) {
			if ( reference_final_index.find(iter->second) == reference_final_index.end() ) {
				start = iter->second;
				reference_final_index.insert(make_pair(start, 0));
				// Check the progress
				if ( seqSizeAln % 1000 == 0 ) {
					printf( "[STEP 3] Reducing the sequence is processing...[%lu/%lu]\n", seqSizeAln, BLOCKLENGTH );
					fflush( stdout );
				}
				seqSizeAln ++;
				break;
			}
		}
		// Decide to terminate or not
		if ( seqSizeAln >= BLOCKLENGTH ) break;
		//// Phase 2
		// Try to find the right next index 
		while ( true ) {
			if ( reference_kmer.find(start + 1) != reference_kmer.end() ) {
				if ( reference_final_index.insert(make_pair(start + 1, 0)).second == false ) {
					break;
				}
				// Check the progress
				if ( seqSizeAln % 1000 == 0 ) {
					printf( "[STEP 3] Reducing the sequences is processing...[%lu/%lu]\n", seqSizeAln, BLOCKLENGTH);
					fflush( stdout );
				}
				// Decide to terminate or not
				if ( ++seqSizeAln >= BLOCKLENGTH ) break;
			} else break;
		}
	}
	// Make the final reference based on the alignment
	uint32_t previous = 0;
	bool init = true;
	for ( auto iter = reference_final_index.begin(); iter != reference_final_index.end(); ++ iter ) {
		// About sequence length
		if ( !init ) {
			if ( iter->first == (previous + 1) ) {
				seqSizeRdc += 1;
			} else {
				string prev = sequence.substr(previous + 1, KMERLENGTH - 1);
				string curr = sequence.substr(iter->first, KMERLENGTH - 1);
				if ( prev.compare(curr) == 0 ) seqSizeRdc += 1;
				else seqSizeRdc += 256;
			}
		} else {
			seqSizeRdc += 256;
			init = false;
		}
		// About reference
		uint64_t encKmer[ENCKMERBUFSIZE] = {0, };
		varTaker(encKmer, iter->first);
		if ( reference_final_kmer.insert(make_pair(iter->first,
					         make_pair(encKmer[0], make_pair(encKmer[1],
					    	 make_pair(encKmer[2], make_pair(encKmer[3],
					    	 make_pair(encKmer[4], make_pair(encKmer[5],
			      		    	 make_pair(encKmer[6], encKmer[7]))))))))).second == false ) {
			printf( "There's a problem on the system...\n" );
			fflush( stdout );
			exit(1);
		}
		// Finalization
		previous = iter->first;
		if ( seqSizeRdc >= BLOCKLENGTH ) break;
	}
	printf( "[STEP 3] The Length of New Sequence : %lu\n", seqSizeRdc );
	printf( "[STEP 3] The Number of New Reference: %lu\n", reference_final_kmer.size() );
	printf( "[STEP 3] Reducing the sequences is done!\n" );
	printf( "---------------------------------------------------------------------\n" );
	fflush( stdout );
	// Write the new reference as binary file
	vector<pair<uint32_t,
	       pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t,
	       pair<uint64_t, pair<uint64_t, pair<uint64_t, uint64_t>>>>>>>>> 
	       reference_vector(reference_final_kmer.begin(), reference_final_kmer.end());
	ofstream f_data_result(filename, ios::binary);
	for ( uint64_t i = 0; i < reference_vector.size(); i ++ ) {
		uint64_t encKmer[ENCKMERBUFSIZE + 1] = {0, };
		encKmer[0] = reference_vector[i].second.first;
		encKmer[1] = reference_vector[i].second.second.first;
		encKmer[2] = reference_vector[i].second.second.second.first;
		encKmer[3] = reference_vector[i].second.second.second.second.first;
		encKmer[4] = reference_vector[i].second.second.second.second.second.first;
		encKmer[5] = reference_vector[i].second.second.second.second.second.second.first;
		encKmer[6] = reference_vector[i].second.second.second.second.second.second.second.first;
		encKmer[7] = reference_vector[i].second.second.second.second.second.second.second.second;
		encKmer[8] = (uint64_t)reference_vector[i].first;
		for ( uint64_t j = 0; j < ENCKMERBUFSIZE + 1; j ++ ) {
			f_data_result.write(reinterpret_cast<char *>(&encKmer[j]), BINARYRWUNIT);
		}
	}
	f_data_result.close();
	printf( "[STEP 5] Writing the k-mers as binary is done\n" );
	printf( "----------------------------------------------------------------\n" );
	fflush( stdout );
}


int main( int argc, char **argv ) {
	char *filenameS = "/mnt/ephemeral/hg19.fasta";
	char *filenameR = "/mnt/ephemeral/hg19Reference256MersFrom1IndexIncluded.bin";
	char *filenameF = "/mnt/ephemeral/hg19Reference256MersFrom1256MBVer7.bin";

	// Read sequence file
	seqReader( filenameS );

	// Read reference file
	refReader( filenameR );

	// Reduce the sequence
	seqShrinker( filenameF );
	return 0;
}
