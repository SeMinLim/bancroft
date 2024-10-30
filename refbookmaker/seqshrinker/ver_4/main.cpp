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
#define GROUPUNIT 16384
#define ENCKMERBUFUNIT 32
#define ENCKMERBUFSIZE 8
#define BINARYRWUNIT 8


uint64_t seqSizeOrg = 0;
uint64_t refSizeOrg = 2849207900;
uint64_t refSizeOrgNew = 0;


vector<string> sequences;
vector<pair<uint32_t, uint64_t>> groups;
map<pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, uint64_t>>>>>>>, 
    uint32_t> reference;
map<pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, uint64_t>>>>>>>, 
    pair<uint32_t, uint32_t>> reference_final;


//// Required Functions
// Time Checker
static inline double timeChecker( void ) {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)(tv.tv_sec) + (double)(tv.tv_usec) / 1000000;
}
// Descending Order Sorter #1
bool descendingOrder_1( pair<uint32_t, uint64_t>& x, pair<uint32_t, uint64_t>& y ) {
	return x.second > y.second;
}
// Descending Order Sorter #2
bool descendingOrder_2( pair<
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
			uint64_t numGroups = seqLine.size() / GROUPUNIT;
			if ( numGroups > 0 ) {
				for ( uint64_t idx = 0; idx < numGroups; idx ++ ) {
					uint64_t start = GROUPUNIT * idx;
					string subseq = seqLine.substr(start, GROUPUNIT);
					sequences.push_back(subseq);
				}
			}
			seqSizeOrg += seqLine.size();
		}
	}
	// Terminate
	f_data_sequences.close();
	printf( "[STEP 1] Reading sequence fasta file is done!\n" );
	fflush( stdout );
}
// Pre-Made Reference File Reader
void refReader( char *filename ) {
	ifstream f_data_reference(filename, ios::binary);
	for ( uint64_t i = 0; i < refSizeOrg; i ++ ) {
		uint64_t encKmer[ENCKMERBUFSIZE + 1] = {0, };
		// Read kmer and occurrence at the same time
		for ( uint64_t j = 0; j < ENCKMERBUFSIZE + 1; j ++ ) {
			f_data_reference.read(reinterpret_cast<char *>(&encKmer[j]), BINARYRWUNIT);
		}
		// Insert 256-Mers and occurrence to Map
		if ( reference.insert(make_pair(make_pair(encKmer[0], make_pair(encKmer[1], make_pair(encKmer[2], 
				      		make_pair(encKmer[3], make_pair(encKmer[4], make_pair(encKmer[5], 
				      		make_pair(encKmer[6], encKmer[7]))))))), (uint32_t)encKmer[8])).second == false ) {
			printf( "There's a problem on reference code book...\n" );
			fflush( stdout );
			exit(1);
		}
		// Check the progress
		if ( i % 1000000 == 0 ) {
			printf( "[STEP 2] Reading pre-made reference file is processing...[%lu/%lu]\n", i, refSizeOrg );
			fflush( stdout );
		}
	}
	// Terminate
	f_data_reference.close();
	printf( "[STEP 2] Reading pre-made reference file is done!\n" );
	fflush( stdout );
}
// Group Selector [4KB Unit]
void groupSelector( void ) {
	for ( uint64_t seqIdx = 0; seqIdx < sequences.size(); seqIdx ++ ) {
		uint64_t start = 0;
		while ( start <= sequences[seqIdx].size() - KMERLENGTH ) {
			string subseq = sequences[seqIdx].substr(start, KMERLENGTH);
			// Encode subsequence first
			uint64_t encSubseq[ENCKMERBUFSIZE] = {0, };
			encoder(subseq, encSubseq);
			// Get the occurrence of the subsequence
			if ( reference.find(make_pair(encSubseq[0], make_pair(encSubseq[1], make_pair(encSubseq[2], 
					    make_pair(encSubseq[3], make_pair(encSubseq[4], make_pair(encSubseq[5], 
					    make_pair(encSubseq[6], encSubseq[7])))))))) != reference.end() ) {
				if ( start == 0 ) {
					groups.push_back(make_pair(seqIdx, reference.at(make_pair(encSubseq[0], make_pair(encSubseq[1], 
								      			make_pair(encSubseq[2], make_pair(encSubseq[3], 
								      			make_pair(encSubseq[4], make_pair(encSubseq[5],
					     	     		      			make_pair(encSubseq[6], encSubseq[7]))))))))));
				} else { 
					groups[seqIdx].second += reference.at(make_pair(encSubseq[0], make_pair(encSubseq[1], 
							 		      make_pair(encSubseq[2], make_pair(encSubseq[3], 
									      make_pair(encSubseq[4], make_pair(encSubseq[5], 
									      make_pair(encSubseq[6], encSubseq[7]))))))));
				}
			}
			start += 1;
		}
		// Check the progress
		printf( "[STEP 3] Couting occurrence is processing...[%lu/%lu]\n", seqIdx, sequences.size() );
		fflush( stdout );
	}
	printf( "[STEP 3] Counting occurrence is done!\n" );
	fflush( stdout );
	// Do sorting by descending order
	sort(groups.begin(), groups.end(), descendingOrder_1);
	printf( "[STEP 4] Sorting the groups is done!\n" );
	fflush( stdout );
	// Delete the primary reference
	reference.clear();
}
// KMC
void kmc( char *filename ) {
	// Run KMC
	uint32_t index = 0;
	for ( uint64_t seqIdx = 0; seqIdx < 65536; seqIdx ++ ) {
		uint64_t start = 0;
		uint32_t remainder = 0;
		while ( start <= sequences[groups[seqIdx].first].size() - KMERLENGTH ) {
			string subseq = sequences[groups[seqIdx].first].substr(start, KMERLENGTH);
			// Encode first
			uint64_t encSubseq[ENCKMERBUFSIZE] = {0, };
			encoder(subseq, encSubseq);
			// Put 256-mers, index, and occurrence to the reference
			if ( reference_final.insert(make_pair(make_pair(encSubseq[0], make_pair(encSubseq[1], make_pair(encSubseq[2], 
					            	      make_pair(encSubseq[3], make_pair(encSubseq[4], make_pair(encSubseq[5], 
					            	      make_pair(encSubseq[6], encSubseq[7]))))))), 
							      make_pair(index, 1))).second == false ) {
				reference_final.at(make_pair(encSubseq[0], make_pair(encSubseq[1], make_pair(encSubseq[2], 
					           make_pair(encSubseq[3], make_pair(encSubseq[4], make_pair(encSubseq[5], 
					           make_pair(encSubseq[6], encSubseq[7])))))))).second += 1;
			} else refSizeOrgNew ++;
			start += 1;
			index += 1;
			// Check the progress
			if ( refSizeOrgNew % 1000000 == 0 ) {
				printf( "[STEP 5] Generating 2-bit encoded k-mer table...[%lu]\n", refSizeOrgNew );
				fflush( stdout );
			}
		}
		// Update the index
		remainder = sequences[groups[seqIdx].first].size() - start;
		index += remainder;
	}
	printf( "[STEP 5] Generating 2-bit encoded k-mer table is done!\n" );
	fflush( stdout );
	// Do sorting
	vector<pair<
	       pair<uint64_t, pair<uint64_t, pair<uint64_t, pair<uint64_t, 
	       pair<uint64_t, pair<uint64_t, pair<uint64_t, uint64_t>>>>>>>, 
	       pair<uint32_t, uint32_t>>> reference_vector(reference_final.begin(), reference_final.end());
	sort(reference_vector.begin(), reference_vector.end(), descendingOrder_2);
	printf( "[STEP 6] Sorting the k-mers through decending order is done!\n" );
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
	printf( "[STEP 7] Writing the k-mers as binary is done\n" );
	fflush( stdout );
}


int main() {
	char *filenameSeq = "/mnt/ephemeral/hg19.fasta";
	char *filenameRef = "/mnt/ephemeral/hg19Reference256MersFrom1OccurIncluded.bin";
	char *filenameOut = "/home/jovyan/hg19Reference256Mers4KBUnit.bin";
	
	// Read sequence file
	seqReader( filenameSeq );

	// Read pre-made reference file
	refReader( filenameRef );

	// Select the groups
	groupSelector();

	// Kmer counting
	double processStartKmc = timeChecker();
	kmc( filenameOut );
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
