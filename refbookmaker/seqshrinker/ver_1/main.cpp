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
#define ENCKMERBUFUNIT 32
#define ENCKMERBUFSIZE 8
#define BINARYRWUNIT 8
#define BLOCKLENGTH 16384


uint64_t seqSizeOrg = 0;
uint64_t seqSizeRdc = 1073741824;
uint64_t refSizeRdc = 0;


string sequenceOrg;
string sequenceRdc;
map<pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, uint64_t>>>>>>>, 
    pair<uint32_t, uint32_t>> reference;


//// Required Functions
// Time Checker
static inline double timeChecker( void ) {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)(tv.tv_sec) + (double)(tv.tv_usec) / 1000000;
}
// Descending Order Sorter
bool descendingOrder( pair<
		     pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, 
		     pair<uint64_t, pair<uint64_t, pair<uint64_t, uint64_t>>>>>>>, pair<uint32_t, uint32_t>>& x, 
		     pair<
		     pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, 
		     pair<uint64_t, pair<uint64_t, pair<uint64_t, uint64_t>>>>>>>, pair<uint32_t, uint32_t>>& y ) {
	return x.second.second > y.second.second;
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
	// Read
	ifstream f_data_sequences(filename);
	while ( getline(f_data_sequences, seqLine) ) {
		if ( seqLine[0] != '>' ) {
			sequenceOrg += seqLine;
			seqSizeOrg += seqLine.size();
		}
	}
	// Terminate
	f_data_sequences.close();
	printf( "----------------------------------------------------------------\n" );
	printf( "[STEP 1] The Original Sequence Length: %lu\n", sequenceOrg.size() );
	printf( "[STEP 1] Reading sequence fasta file is done!\n" );
	printf( "----------------------------------------------------------------\n" );
	fflush( stdout );
}
// Sequence Shrinker
void seqShrinker( void ) {
	for ( uint64_t cnt = 0; cnt < (seqSizeRdc / BLOCKLENGTH); cnt ++ ) {
		// Set the block
		uint64_t start = cnt * (BLOCKLENGTH * 2);
		string block = sequenceOrg.substr(start, BLOCKLENGTH);
		// Store the block to the new sequence
		sequenceRdc += block;
	}
	// Terminate
	printf( "[STEP 2] The Reduced Sequence Length : %lu\n", sequenceRdc.size() );
	printf( "[STEP 2] Reducing sequence is done!\n" );
	printf( "----------------------------------------------------------------\n" );
	fflush( stdout );
}
// KMC
void kmc( char *filename ) {
	// Run KMC
	uint64_t start = 0;
	while ( start <= sequenceRdc.size() - KMERLENGTH ) {
		string subseq = sequenceRdc.substr(start, KMERLENGTH);
		// Encode first
		uint64_t encSubseq[ENCKMERBUFSIZE] = {0, };
		encoder(subseq, encSubseq);
		// Put 256-mers and index & Update the number of occurrence
		if ( reference.insert(make_pair(make_pair(encSubseq[0], make_pair(encSubseq[1], make_pair(encSubseq[2], 
				      make_pair(encSubseq[3], make_pair(encSubseq[4], make_pair(encSubseq[5], 
				      make_pair(encSubseq[6], encSubseq[7]))))))), make_pair(start, 1))).second == false ) {
			reference.at(make_pair(encSubseq[0], make_pair(encSubseq[1], make_pair(encSubseq[2], 
				     make_pair(encSubseq[3], make_pair(encSubseq[4], make_pair(encSubseq[5],
				     make_pair(encSubseq[6], encSubseq[7])))))))).second += 1;
		} else refSizeRdc ++;
		start += 1;
		// Check the progress
		if ( refSizeRdc % 1000000 == 0 ) {
			printf( "[STEP 3] Generating 2-bit encoded k-mer table...[%lu]\n", refSizeRdc );
			fflush( stdout );
		}
	}
	printf( "[STEP 3] Generating 2-bit encoded k-mer table is done!\n" );		
	printf( "----------------------------------------------------------------\n" );
	fflush( stdout );
	// Do sorting
	vector<pair<
	       pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, 
	       pair<uint64_t, pair<uint64_t, pair<uint64_t, uint64_t>>>>>>>, 
	       pair<uint32_t, uint32_t>>> reference_vector(reference.begin(), reference.end());
	sort(reference_vector.begin(), reference_vector.end(), descendingOrder);
	printf( "[STEP 4] Sorting the k-mers through decending order is done!\n" );
	printf( "----------------------------------------------------------------\n" );
	fflush( stdout );
	// Write the reference as binary file
	ofstream f_data_result(filename, ios::binary);
	for ( uint64_t i = 0; i < reference_vector.size(); i ++ ) {
		uint64_t encKmer[ENCKMERBUFSIZE + 1] = {0, };
		encKmer[0] = reference_vector[i].first.first;
		encKmer[1] = reference_vector[i].first.second.first;
		encKmer[2] = reference_vector[i].first.second.second.first;
		encKmer[3] = reference_vector[i].first.second.second.second.first;
		encKmer[4] = reference_vector[i].first.second.second.second.second.first;
		encKmer[5] = reference_vector[i].first.second.second.second.second.second.first;
		encKmer[6] = reference_vector[i].first.second.second.second.second.second.second.first;
		encKmer[7] = reference_vector[i].first.second.second.second.second.second.second.second;
		encKmer[8] = (uint64_t)reference_vector[i].second.first;
		for ( uint64_t j = 0; j < ENCKMERBUFSIZE + 1; j ++ ) {
			f_data_result.write(reinterpret_cast<char *>(&encKmer[j]), BINARYRWUNIT);
		}
	}
	f_data_result.close();
	printf( "[STEP 5] Writing the k-mers as binary is done\n" );
	printf( "----------------------------------------------------------------\n" );
	fflush( stdout );
}


int main() {
	char *filenameIn = "/mnt/ephemeral/hg19.fasta";
	char *filenameOut = "/mnt/ephemeral/hg19Reference256MersFrom14KBVer1.bin";
	
	// Read sequence file
	seqReader( filenameIn );

	// Reduce the sequence 
	seqShrinker();

	// Kmer counting
	double processStartKmc = timeChecker();
	kmc( filenameOut );
	double processFinishKmc = timeChecker();
	double elapsedTimeKmc = processFinishKmc - processStartKmc;
	
	printf( "KMC RESULT\n" );
	printf( "KMER [Total]: %ld\n", seqSizeRdc );
	printf( "KMER [Count]: %ld\n", refSizeRdc );
	printf( "KMER [Percentage]: %0.4f\n", ((double)refSizeRdc / (double)seqSizeRdc) * (double)100.00 );
	printf( "Reference Book [Size]: %0.4f GB\n", (double)((refSizeRdc * 512) + (refSizeRdc * 30)) / 8 / 1024 / 1024 / 1024 );
	printf( "Elapsed Time: %lf\n", elapsedTimeKmc );
	printf( "----------------------------------------------------------------\n" );
	
	return 0;
}
