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
#include <map>
using namespace std;


#define KMERLENGTH 256
#define ENCKMERBUFUNIT 32
#define ENCKMERBUFSIZE 8
#define BINARYRWUNIT 8


vector<string> hg38;
map<pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, uint64_t>>>>>>>, uint64_t> hg19;


uint64_t numHG19 = 2838231473;
uint64_t numHG19HG38 = 0;


// Function for Sorting
bool decendingOrder( pair<
		     pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, 
		     pair<uint64_t, pair<uint64_t, pair<uint64_t, uint64_t>>>>>>>, uint64_t>& x, 
		     pair<
		     pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, 
		     pair<uint64_t, pair<uint64_t, pair<uint64_t, uint64_t>>>>>>>, uint64_t>& y ) {
	return x.second > y.second;
}

// File Readers
void hg38Reader( char *filename ) {
	string hg38Line;

	ifstream f_data_hg38(filename);
	while ( getline(f_data_hg38, hg38Line) ) {
		if ( hg38Line[0] != '>' ) {
			hg38.push_back(hg38Line);
		}
	}

	f_data_hg38.close();
}
void hg19Reader( char *filename ) {
	ifstream f_data_hg19(filename, ios::binary);
	for ( uint64_t i = 0; i < numHG19; i ++ ) {
		uint64_t encKmer[ENCKMERBUFSIZE] = {0, };
		// Read
		for ( uint64_t j = 0; j < ENCKMERBUFSIZE; j ++ ) {
			f_data_hg19.read(reinterpret_cast<char *>(&encKmer[j]), BINARYRWUNIT);
		}
		// Insert 256-Mers to Map
		if ( hg19.insert(make_pair(make_pair(encKmer[0], make_pair(encKmer[1], make_pair(encKmer[2], 
				 make_pair(encKmer[3], make_pair(encKmer[4], make_pair(encKmer[5], 
				 make_pair(encKmer[6], encKmer[7]))))))), 1)).second == false ) {
			printf( "There's a problem on reference code book...\n" );
			fflush( stdout );
			exit(1);
		}
		// Check the progress
		if ( i % 1000000 == 0 ) {
			printf( "Reference: %lu\n", i );
			fflush( stdout );
		}
	}

	f_data_hg19.close();
}

// Encoder & Decoder
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

// HG19HG38 Maker
void refMaker( char *filename ) {
	// Generate the draft of hg19hg38_2
	for ( uint64_t seqIdx = 0; seqIdx < hg38.size(); seqIdx ++ ) {
		uint64_t start = 0;
		while ( start <= hg38[seqIdx].size() - KMERLENGTH ) {
			string subseq = hg38[seqIdx].substr(start, KMERLENGTH);

			// Encode first
			uint64_t encSubseq[ENCKMERBUFSIZE] = {0, };
			encoder(subseq, encSubseq);

			// Check if there is the same k-mer in hg38
			if ( hg19.find(make_pair(encSubseq[0], make_pair(encSubseq[1], make_pair(encSubseq[2], 
				       make_pair(encSubseq[3], make_pair(encSubseq[4], make_pair(encSubseq[5], 
				       make_pair(encSubseq[6], encSubseq[7])))))))) != hg19.end() ) {
				hg19.at(make_pair(encSubseq[0], make_pair(encSubseq[1], make_pair(encSubseq[2],
					make_pair(encSubseq[3], make_pair(encSubseq[4], make_pair(encSubseq[5],
					make_pair(encSubseq[6], encSubseq[7])))))))) += 1;
			}
			
			start += 1;
		}
		
		printf( "Scanning #%lu Sequences is Done!\n", seqIdx );
		fflush( stdout );
	}
	printf( "Generating hg19hg38_2 is Done!\n" );
	fflush( stdout );

	// Delete the elements that has 1 occurrence
	for ( auto iter = hg19.begin(); iter != hg19.end(); ) {
		if ( iter->second == 1 ) hg19.erase(iter++);
		else ++iter;
	}

	// Do sorting
	vector<pair<
	       pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, 
	       pair<uint64_t, pair<uint64_t, pair<uint64_t, uint64_t>>>>>>>, 
	       uint64_t>> hg19_vector(hg19.begin(), hg19.end());
	sort(hg19_vector.begin(), hg19_vector.end(), decendingOrder);
	numHG19HG38 = hg19_vector.size();
	printf( "Sorting hg19hg38_2 is Done!\n" );
	fflush( stdout );

	// Write the result
	ofstream f_data_result(filename, ios::binary);
	for ( uint64_t i = 0; i < hg19_vector.size(); i ++ ) {
		uint64_t encKmer[ENCKMERBUFSIZE+1] = {0, };
		encKmer[0] = hg19_vector[i].first.first;
		encKmer[1] = hg19_vector[i].first.second.first;
		encKmer[2] = hg19_vector[i].first.second.second.first;
		encKmer[3] = hg19_vector[i].first.second.second.second.first;
		encKmer[4] = hg19_vector[i].first.second.second.second.second.first;
		encKmer[5] = hg19_vector[i].first.second.second.second.second.second.first;
		encKmer[6] = hg19_vector[i].first.second.second.second.second.second.second.first;
		encKmer[7] = hg19_vector[i].first.second.second.second.second.second.second.second;
		encKmer[8] = hg19_vector[i].second;
		for ( uint64_t j = 0; j < ENCKMERBUFSIZE + 1; j ++ ) {
			f_data_result.write(reinterpret_cast<char *>(&encKmer[j]), BINARYRWUNIT);
		}
	}
	f_data_result.close();
	printf( "Writing hg19hg38_2 is Done\n" );
	fflush( stdout );
}


int main( int argc, char **argv ) {
	char *filenamehg38 = "/mnt/ephemeral/hg38.fasta";
	char *filenamehg19 = "/mnt/ephemeral/hg19RefBook256Mers_2.bin";
	char *filenamehg19hg38 = "/mnt/ephemeral/hg19hg38Reference256Mers_2.bin";

	// Read hg38
	hg38Reader( filenamehg38 );
	printf( "Reading hg38 is Done!\n" );
	fflush( stdout );

	// Read hg19
	hg19Reader( filenamehg19 );
	printf( "Reading hg19 is Done!\n" );
	fflush( stdout );

	// Make hg19hg38
	refMaker( filenamehg19hg38 ); 
	printf( "--------------------------------------------\n" );
	printf( "The Number of K-Mer in HG19    : %lu\n", numHG19 );
	printf( "The Number of K-Mer in HG19HG38: %lu\n", numHG19HG38 );

	// Verification	
	string decKmer;
	decKmer.reserve(KMERLENGTH);
	uint64_t encKmer[ENCKMERBUFSIZE+1] = {0, };
	
	ifstream f_data_hg19hg38(filenamehg19hg38, ios::binary);
	for ( uint64_t i = 0; i < ENCKMERBUFSIZE + 1; i ++ ) {
		f_data_hg19hg38.read(reinterpret_cast<char *>(&encKmer[i]), BINARYRWUNIT);
	}
	
	decoder(encKmer, decKmer);

	printf( "%s", decKmer.c_str() );
	printf( "\t" );
	printf( "%lu\n", encKmer[ENCKMERBUFSIZE] );
	
	f_data_hg19hg38.close();

	return 0;
}
