#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <unordered_map>
#include <map>
using namespace std;


#define KMERLENGTH 12


unordered_map<string, uint32_t> referenceOG_Unordered;
map<uint32_t, string> referenceOG_Ordered;
map<uint32_t, string> referenceDB;


uint64_t cntOG = 0;
uint64_t cntDB = 0;
uint64_t independent = 0;


static inline double timeChecker( void ) {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)(tv.tv_sec) + (double)(tv.tv_usec) / 1000000;
}

void refBookReader( char *filename ) {
	string refLine;

	ifstream f_data_reference(filename);
	if ( !f_data_reference.is_open() ) {
		printf( "File not found: %s\n", filename );
		exit(1);
	}
	
	uint32_t index = 0;
	while ( getline(f_data_reference, refLine) ) {
		// Kmer
		string kmer;
		kmer.reserve(KMERLENGTH);		
		for ( size_t i = 0; i < KMERLENGTH; i ++ ) {
			kmer.push_back(refLine[i]);
		}
		
		// Reference updating
		referenceOG_Unordered.insert(make_pair(kmer, index));
		referenceOG_Ordered.insert(make_pair(index, kmer));
		
		// Etc
		index ++;
		refLine.clear();
		refLine.shrink_to_fit();
	}

	f_data_reference.close();

}

void fastaWriter( char *filename ) {
	ofstream f_data_result(filename);
	if ( !f_data_result.is_open() ) {
		printf( "File not found: %s\n", filename );
		exit(1);
	}

	for ( auto iter = referenceDB.begin(); iter !=referenceDB.end(); iter ++ ) {
		f_data_result << iter->second << "\n";
	}
	f_data_result << endl;

	f_data_result.close();
}

void debruijn( void ) {
	cntOG = referenceOG_Ordered.size();
	
	bool init = true;
	string dna = "TGCA";
	string kmer;
	string subkmer;
	kmer.reserve(KMERLENGTH);
	subkmer.reserve(KMERLENGTH-1);
	while ( cntDB < cntOG ) {
		if ( init ) {
			auto begin = referenceOG_Ordered.begin();
			
			kmer = begin->second;
			subkmer = kmer.substr(1, KMERLENGTH-1);
			
			referenceOG_Unordered.erase(kmer);
			begin = referenceOG_Ordered.erase(begin);
			if ( referenceDB.insert(make_pair(cntDB, kmer)).second == false) {
				printf( "Something is wrong\n" );
			}
			
			cntDB ++;
			independent ++;
		}
		
		uint8_t notFound = 0;
		uint8_t dnaIdx = 0;
		uint32_t mapIdx = cntOG;
		for ( size_t i = 0; i < 4; i ++ ) {
			string kmerTMP = subkmer + dna[i];
			if ( referenceOG_Unordered.find(kmerTMP) != referenceOG_Unordered.end() ) {
				uint32_t idx = referenceOG_Unordered.at(kmerTMP);
				if ( mapIdx > idx ) {
					mapIdx = idx;
					dnaIdx = i;
				}
			} else notFound ++;
		}

		if ( notFound < 4 ) {
			kmer = subkmer + dna[dnaIdx];
			subkmer = kmer.substr(1, KMERLENGTH-1);
			
			referenceOG_Unordered.erase(kmer);
			referenceOG_Ordered.erase(mapIdx);
			if ( referenceDB.insert(make_pair(cntDB, kmer)).second == false ) {
				printf( "Something is wrong\n" );
			}

			cntDB ++;
			init = false;
		} else init = true;

		if ( cntDB % 1000000 == 0 ) {
			printf( "Independent: %ld, De Bruijn: %ld\n", independent, cntDB );
		}
	}
}


int main() {
	char *filenameR = "../../data/references/hg19RefBook12Mers.txt";
	char *filenameD = "../../data/references/hg19RefBook12MersDeBruijn.txt";

	// Read reference file
	refBookReader( filenameR );
	printf( "Reading Reference File is Finished!\n" );
	fflush( stdout );

	// Make a De Bruijn style reference code book
	debruijn();
	printf( "Changing Reference File as De Bruijn Style is Finished!\n" );
	fflush( stdout );

	// Write the De Bruijn style reference code book
	fastaWriter( filenameD );
	printf( "Writing the De Bruijn Style Reference Code Book is Finished!\n" );
	fflush( stdout );

	// Results
	printf( "------------------------------------------------------------\n" );
	printf( "# Origin Kmer:    %ld\n", independent );
	printf( "# De Bruijn Kmer: %ld\n", cntDB );
	printf( "------------------------------------------------------------\n" );
	fflush( stdout );

	return 0;
}
