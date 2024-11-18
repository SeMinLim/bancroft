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


#define KMERLENGTH 32
#define ENCKMERBUFUNIT 32
#define ENCKMERBUFSIZE 1
#define BINARYRWUNIT 8


uint64_t seqSizeOrg = 0;
uint64_t refSizeOrg = 0;


vector<string> sequences;
map<uint64_t, pair<uint32_t, uint32_t>> reference;


//// Required Functions
// Time Checker
static inline double timeChecker( void ) {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)(tv.tv_sec) + (double)(tv.tv_usec) / 1000000;
}
// Descending Order Sorter
bool descendingOrder( pair<
		      uint64_t, pair<uint32_t, uint32_t>>& x, 
		      pair<
		      uint64_t, pair<uint32_t, uint32_t>>& y ) {
	return x.second.second > y.second.second;
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
// Sequence File Reader
void fastaReader( char *filename ) {
	string seqLine;
	// Read
	ifstream f_data_sequences(filename);
	while ( getline(f_data_sequences, seqLine) ) {
		if ( seqLine[0] != '>' ) {
			sequences.push_back(seqLine);
			seqSizeOrg += seqLine.size();
		}
	}
	// Terminate
	f_data_sequences.close();
	printf( "[STEP 1] Reading sequence fasta file is done!\n" );
	fflush( stdout );
}
// KMC
void kmc( char *filename ) {
	uint32_t index = 0;
	for ( uint64_t seqIdx = 0; seqIdx < sequences.size(); seqIdx ++ ) {
		uint64_t start = 0;
		uint32_t remainder = 0;
		while ( start <= sequences[seqIdx].size() - KMERLENGTH ) {
			string subseq = sequences[seqIdx].substr(start, KMERLENGTH);
			// Encode subsequence first
			uint64_t encSubseq = 0;
			encoder(subseq, encSubseq);
			// Put encoded subsequence, index, and occurrence to the reference
			if ( reference.insert(make_pair(encSubseq, make_pair(index, 1))).second == false ) {
				reference.at(encSubseq).second += 1;
			} else refSizeOrg ++;
			start += 1;
			index += 1;
			// Check the progress
			if ( refSizeOrg % 1000000 == 0 ) {
				printf( "[STEP 2] Generating 2-bit encoded k-mer table...[%lu]\n", refSizeOrg );
				fflush( stdout );
			}
		}
		// Update the index
		remainder = sequences[seqIdx].size() - start;
		index += remainder;
	}
	printf( "[STEP 2] Generating 2-bit encoded k-mer table is done!\n" );
	fflush( stdout );
	// Do sorting
	vector<pair<
	       uint64_t, 
	       pair<uint32_t, uint32_t>>> reference_vector(reference.begin(), reference.end());
	sort(reference_vector.begin(), reference_vector.end(), descendingOrder);
	printf( "[STEP 3] Sorting the k-mers through decending order is done!\n" );
	fflush( stdout );
	// Write the reference as binary file
	ofstream f_data_result(filename, ios::binary);
	for ( uint64_t i = 0; i < reference_vector.size(); i ++ ) {
		uint64_t encKmer[ENCKMERBUFSIZE + 1] = {0, };
		encKmer[0] = reference_vector[i].first;
		encKmer[1] = (uint64_t)reference_vector[i].second.first;
		for ( uint64_t j = 0; j < ENCKMERBUFSIZE + 1; j ++ ) {
			f_data_result.write(reinterpret_cast<char *>(&encKmer[j]), BINARYRWUNIT);
		}
	}
	f_data_result.close();
	printf( "[STEP 4] Writing the k-mers as binary is done\n" );
	fflush( stdout );
}


int main() {
	char *filenameIn = "/mnt/ephemeral/sequence/HG19.fasta";
	char *filenameOut = "/mnt/ephemeral/reference/HG19Reference032MersFrom1IndexIncluded.bin";
	
	// Read sequence file
	fastaReader( filenameIn );

	// Kmer counting
	double processStartKmc = timeChecker();
	kmc( filenameOut );
	double processFinishKmc = timeChecker();
	double elapsedTimeKmc = processFinishKmc - processStartKmc;
	
	printf( "--------------------------------------------------------------------------------\n" );
	printf( "KMC RESULT\n" );
	printf( "KMER [Total]: %ld\n", seqSizeOrg );
	printf( "KMER [Count]: %ld\n", refSizeOrg );
	printf( "KMER [Percentage]: %0.4f\n", (double)refSizeOrg / (double)seqSizeOrg * 100.00 );
	printf( "Reference Book [Size]: %0.4f GB\n", 
		(((double)refSizeOrg * 512.00) + ((double)refSizeOrg * 32.00)) / 8.00 / 1024.00 / 1024.00 / 1024.00 );
	printf( "Elapsed Time: %lf\n", elapsedTimeKmc );
	printf( "--------------------------------------------------------------------------------\n" );
	
	return 0;
}
