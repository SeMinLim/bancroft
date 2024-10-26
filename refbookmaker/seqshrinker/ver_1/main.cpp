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
#define SHRINKLENGTH 4096
#define ENCKMERBUFUNIT 32
#define ENCKMERBUFSIZE 8
#define BINARYRWUNIT 8


uint64_t seqSizeOrg = 0;
uint64_t refSizeOrg = 0;


vector<string> sequences;
map<pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, uint64_t>>>>>>>, 
    uint32_t> reference;


// Required Functions
static inline double timeChecker( void ) {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)(tv.tv_sec) + (double)(tv.tv_usec) / 1000000;
}
bool decendingOrder( pair<
		     pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, 
		     pair<uint64_t, pair<uint64_t, pair<uint64_t, uint64_t>>>>>>>, uint32_t>& x, 
		     pair<
		     pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, 
		     pair<uint64_t, pair<uint64_t, pair<uint64_t, uint64_t>>>>>>>, uint32_t>& y ) {
	return x.second > y.second;
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

// File Reader + Sequence Shrinker
void fastaReader( char *filename ) {
	uint64_t seqCounter = 0;
	string seqLine;

	ifstream f_data_sequences(filename);
	while ( getline(f_data_sequences, seqLine) ) {
		if ( seqLine[0] != '>' ) {
			// Shrink a sequence first
			uint64_t pointer = 0;
			string seqShrinked;
			while ( true ) {
				if ( (pointer + SHRINKLENGTH) > seqLine.size()  ) {
					break;
				}
				else {
					for ( uint64_t i = pointer; i < pointer + SHRINKLENGTH; i ++ ) {
						seqShrinked += seqLine[i];
					}
					pointer += (SHRINKLENGTH * 2);
				}
			}
			
			// Store the shrinked sequence to the vector then
			sequences.push_back(seqShrinked);
			seqSizeOrg += seqShrinked.size();

			// Check the progress
			printf( "[STEP 1] Shrinking a sequence #%lu is done!...[%lu]\n", seqCounter, seqShrinked.size() );
			fflush( stdout );
			seqCounter ++;
		}
	}

	f_data_sequences.close();
	printf( "[STEP 1] Reading sequence fasta file is done!\n" );
	fflush( stdout );
}

void kmc( char *filename ) {
	// Do 2-bit Encoding & Put all 256-mers to Map & Count the number of occurrence
	for ( uint64_t seqIdx = 0; seqIdx < sequences.size(); seqIdx ++ ) {
		uint64_t start = 0;
		while ( start <= sequences[seqIdx].size() - KMERLENGTH ) {
			string subseq = sequences[seqIdx].substr(start, KMERLENGTH);

			// Encode first
			uint64_t encSubseq[ENCKMERBUFSIZE] = {0, };
			encoder(subseq, encSubseq);

			// Put 256-mers and update the number of occurrence
			if ( reference.insert(make_pair(make_pair(encSubseq[0], make_pair(encSubseq[1], make_pair(encSubseq[2], 
					      make_pair(encSubseq[3], make_pair(encSubseq[4], make_pair(encSubseq[5], 
					      make_pair(encSubseq[6], encSubseq[7]))))))), 1)).second == false ) {
				reference.at(make_pair(encSubseq[0], make_pair(encSubseq[1], make_pair(encSubseq[2], 
					     make_pair(encSubseq[3], make_pair(encSubseq[4], make_pair(encSubseq[5],
					     make_pair(encSubseq[6], encSubseq[7])))))))) += 1;
			} else refSizeOrg ++;
			start += 1;

			// Check the progress
			if ( refSizeOrg % 1000000 == 0 ) {
				printf( "[STEP 2] Generating 2-bit encoded k-mer table...[%lu]\n", refSizeOrg );
				fflush( stdout );
			}
		}
	}
	printf( "[STEP 2] Generating 2-bit encoded k-mer table is done!\n" );
	fflush( stdout );
	
	// Do sorting
	vector<pair<
	       pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, 
	       pair<uint64_t, pair<uint64_t, pair<uint64_t, uint64_t>>>>>>>, 
	       uint32_t>> reference_vector(reference.begin(), reference.end());
	sort(reference_vector.begin(), reference_vector.end(), decendingOrder);
	printf( "[STEP 3] Sorting the k-mers through decending order is done!\n" );
	fflush( stdout );

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
}


int main() {
	char *filenameIn = "/mnt/smartssd0/semin/hg19.fasta";
	char *filenameOut = "/mnt/smartssd0/semin/hg19ShrinkedReference256Mers.bin";
	
	// Read sequence file
	fastaReader( filenameIn );

	// Kmer counting
	double processStartKmc = timeChecker();
	kmc( filenameOut);
	double processFinishKmc = timeChecker();
	double elapsedTimeKmc = processFinishKmc - processStartKmc;
	
	printf( "--------------------------------------------\n" );
	printf( "KMC RESULT\n" );
	printf( "KMER [Total]: %ld\n", seqSizeOrg );
	printf( "KMER [Count]: %ld\n", refSizeOrg );
	printf( "KMER [Percentage]: %0.8f\n", (double)(refSizeOrg / seqSizeOrg) * 100 );
	printf( "Reference Book [Size]: %0.4f GB\n", (double)((refSizeOrg * 512) + (refSizeOrg * 32)) / 8 / 1024 / 1024 / 1024 );
	printf( "Elapsed Time: %lf\n", elapsedTimeKmc );
	printf( "--------------------------------------------\n" );
	
	
	return 0;
}
