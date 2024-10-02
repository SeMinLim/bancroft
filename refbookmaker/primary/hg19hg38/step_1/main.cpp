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


vector<string> hg38;
unordered_map<string, uint64_t> hg19;
unordered_map<string, uint64_t> hg19hg38;


uint64_t numHG19 = 0;
uint64_t numHG19HG38 = 0;


// Function for Sorting
bool decendingOrder( pair<string, uint64_t>& x, pair<string, uint64_t>& y ) {
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
	string hg19Line;

	ifstream f_data_hg19(filename);
	while ( getline(f_data_hg19, hg19Line) ) {
		// K-Mer
		string kmer;
		kmer.reserve(KMERLENGTH);
		for ( uint64_t i = 0; i < KMERLENGTH; i ++ ) {
			kmer.push_back(hg19Line[i]);
		}

		// Count
		string cntS;
		for ( uint64_t i = KMERLENGTH + 1; i < hg19Line.size(); i ++ ) {
			cntS.push_back(hg19Line[i]);
		}
		uint64_t cntI = atoi(cntS.c_str());

		// Unordered Map
		if ( hg19.insert(make_pair(kmer, cntI)).second == false ) {
			printf( "There's an issue on reference code book\n" );
			exit(1);
		}

		numHG19 ++;
	}

	f_data_hg19.close();
}

// HG19HG38 Maker
void refMaker( char *filename ) {
	// Generate the draft of hg19hg38_1
	for ( uint64_t seqIdx = 0; seqIdx < hg38.size() - 1; seqIdx ++ ) {
		uint64_t start = 0;
		while ( start <= hg38[seqIdx].size() - KMERLENGTH ) {
			string subseq = hg38[seqIdx].substr(start, KMERLENGTH);

			// Check if there is the same k-mer in hg38
			if ( hg19.find(subseq) != hg19.end() ) {
				if ( hg19hg38.insert(make_pair(subseq, 1)).second == false ) {
					hg19hg38.at(subseq) += 1;
				}
			}

			start += 1;
		}
		
		printf( "Scanning #%lu Sequences is Done!\n", seqIdx );
		fflush( stdout );
	}
	printf( "Generating hg19hg38_1 is Done!\n" );
	fflush( stdout );

	// Update the number of count of hg19hg38
	for ( auto iter = hg19hg38.begin(); iter != hg19hg38.end(); ++iter ) {
		if ( hg19.find(iter->first) != hg19.end() ) {
			iter->second = iter->second + hg19.at(iter->first);
		}
	}
	printf( "Composing hg19hg38_1 is Done!\n" );
	fflush( stdout );

	// Do sorting
	vector<pair<string, uint64_t>> hg19hg38_vector(hg19hg38.begin(), hg19hg38.end());
	sort(hg19hg38_vector.begin(), hg19hg38_vector.end(), decendingOrder);
	numHG19HG38 = hg19hg38_vector.size();
	printf( "Sorting hg19hg38_1 is Done!\n" );
	fflush( stdout );

	// Write the result
	ofstream f_data_result(filename);
	for ( uint64_t i = 0; i < hg19hg38_vector.size(); i ++ ) {
		f_data_result << hg19hg38_vector[i].first << "\t" << hg19hg38_vector[i].second << "\n";
	}
	f_data_result << endl;
	f_data_result.close();
	printf( "Writing hg19hg38_1 is Done\n" );
	fflush( stdout );
}


int main( int argc, char **argv ) {
	char *filenamehg38 = "/home/semin/dna_compressor/data/references/hg38/hg38.fasta";
	char *filenamehg19 = "/mnt/smartssd0/semin/hg19RefBook256Mers_1.txt";
	char *filenamehg19hg38 = "/mnt/smartssd0/semin/hg19hg38Reference256Mers_1.txt";

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

	return 0;
}
