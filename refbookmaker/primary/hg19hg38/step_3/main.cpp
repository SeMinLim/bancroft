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


vector<pair<pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, uint64_t>>>>>>>, uint16_t>> hg19hg38;


uint64_t numHG19HG38Expect = 2825518939;
uint64_t numHG19HG38Actual = 0;


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

// Function for Sorting
bool decendingOrder( pair<
		     pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, 
		     pair<uint64_t, pair<uint64_t, pair<uint64_t, uint64_t>>>>>>>, uint16_t>& x, 
		     pair<
		     pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, 
		     pair<uint64_t, pair<uint64_t, pair<uint64_t, uint64_t>>>>>>>, uint16_t>& y ) {
	return x.second > y.second;
}

// File Readers
void hg19hg38FirstReader( char *filename ) {
	string hg19hg38Line;

	ifstream f_data_hg19hg38(filename);
	while ( getline(f_data_hg19hg38, hg19hg38Line) ) {
		// K-Mer
		string kmer;
		kmer.reserve(KMERLENGTH);
		uint64_t encKmer[ENCKMERBUFSIZE] = {0, };
		for ( uint64_t i = 0; i < KMERLENGTH; i ++ ) {
			kmer.push_back(hg19hg38Line[i]);
		}
		encoder(kmer, encKmer);

		// Count
		string cntS;
		for ( uint64_t i = KMERLENGTH + 1; i < hg19hg38Line.size(); i ++ ) {
			cntS.push_back(hg19hg38Line[i]);
		}
		uint16_t cntI = atoi(cntS.c_str());

		// Vector
		hg19hg38.push_back(make_pair(
				   make_pair(encKmer[0], make_pair(encKmer[1], make_pair(encKmer[2], make_pair(encKmer[3], 
				   make_pair(encKmer[4], make_pair(encKmer[5], make_pair(encKmer[6], encKmer[7]))))))), 
				   cntI));

		// Check the progress
		if ( numHG19HG38Actual % 1000000 == 0 ) {
			printf( "[FIRST] Reference: %lu\n", numHG19HG38Actual );
			fflush( stdout );
		}
		numHG19HG38Actual ++;		
	}

	f_data_hg19hg38.close();
}
void hg19hg38SecondReader( char *filename ) {
	ifstream f_data_hg19hg38(filename, ios::binary);
	while ( !f_data_hg19hg38.eof() ) {
		// Read
		uint64_t encoded[ENCKMERBUFSIZE+1] = {0, };
		for ( uint64_t i = 0; i < ENCKMERBUFSIZE + 1; i ++ ) {
			f_data_hg19hg38.read(reinterpret_cast<char *>(&encoded[i]), BINARYRWUNIT);
		}

		// Vector
		hg19hg38.push_back(make_pair(
				   make_pair(encoded[0], make_pair(encoded[1], make_pair(encoded[2], make_pair(encoded[3], 
				   make_pair(encoded[4], make_pair(encoded[5], make_pair(encoded[6], encoded[7]))))))), 
				   (uint16_t)encoded[8]));

		// Check the progress
		if ( numHG19HG38Actual % 1000000 == 0 ) {
			printf( "[SECOND] Reference: %lu\n", numHG19HG38Actual );
			fflush( stdout );
		}
		numHG19HG38Actual ++;
	}

	f_data_hg19hg38.close();
}

// HG19HG38 Maker
void refMaker( char *filename ) {
	// Do sorting
	sort(hg19hg38.begin(), hg19hg38.end(), decendingOrder);
	printf( "Sorting hg19hg38 is Done!\n" );
	fflush( stdout );

	// Write the result
	ofstream f_data_result(filename, ios::binary);
	for ( uint64_t i = 0; i < hg19hg38.size(); i ++ ) {
		uint64_t encKmer[ENCKMERBUFSIZE] = {0, };
		encKmer[0] = hg19hg38[i].first.first;
		encKmer[1] = hg19hg38[i].first.second.first;
		encKmer[2] = hg19hg38[i].first.second.second.first;
		encKmer[3] = hg19hg38[i].first.second.second.second.first;
		encKmer[4] = hg19hg38[i].first.second.second.second.second.first;
		encKmer[5] = hg19hg38[i].first.second.second.second.second.second.first;
		encKmer[6] = hg19hg38[i].first.second.second.second.second.second.second.first;
		encKmer[7] = hg19hg38[i].first.second.second.second.second.second.second.second;
		
		uint16_t cnt = hg19hg38[i].second;
		
		for ( uint64_t j = 0; j < ENCKMERBUFSIZE; j ++ ) {
			f_data_result.write(reinterpret_cast<char *>(&encKmer[j]), BINARYRWUNIT);
		}
		
		if ( i == 0 ) {
			string decKmer;
			decKmer.reserve(KMERLENGTH);
			decoder(encKmer, decKmer);
			printf( "%s", decKmer.c_str() );
			printf( "\t" );
			printf( "%u\n", cnt );
		}
	}
	f_data_result.close();
	printf( "Writing hg19hg38_2 is Done\n" );
	fflush( stdout );
}


int main( int argc, char **argv ) {
	char *filenamehg19hg38First = "/home/jovyan/dna_compressor/data/references/hg19hg38/hg19hg38Reference256Mers_1.txt";
	char *filenamehg19hg38Second = "/home/jovyan/dna_compressor/data/references/hg19hg38/hg19hg38Reference256Mers_2.bin";
	char *filenamehg19hg38Final = "/mnt/ephemeral/hg19hg38Reference256Mers.bin";

	// Read hg38
	hg19hg38FirstReader( filenamehg19hg38First );

	// Read hg19
	hg19hg38SecondReader( filenamehg19hg38Second );

	// Make hg19hg38
	refMaker( filenamehg19hg38Final ); 
	printf( "--------------------------------------------\n" );
	printf( "The Number of K-Mer in HG19HG38 [Expect]: %lu\n", numHG19HG38Expect );
	printf( "The Number of K-Mer in HG19HG38 [Actual]: %lu\n", numHG19HG38Actual );

	// Verification	
	string decKmer;
	decKmer.reserve(KMERLENGTH);
	uint64_t encKmer[ENCKMERBUFSIZE] = {0, };
	
	ifstream f_data_hg19hg38(filenamehg19hg38Final, ios::binary);
	for ( uint64_t i = 0; i < ENCKMERBUFSIZE; i ++ ) {
		f_data_hg19hg38.read(reinterpret_cast<char *>(&encKmer[i]), BINARYRWUNIT);
	}
	
	decoder(encKmer, decKmer);

	printf( "%s\n", decKmer.c_str() );
	
	f_data_hg19hg38.close();

	return 0;
}
