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
using namespace std;


vector<uint64_t> length;


uint64_t counter = 0;
uint64_t sequences = 0;
uint64_t totalLength = 0;


void fastaReader( char *filename ) {
	string seqLine;

	ifstream f_data_sequences(filename);
	if ( !f_data_sequences.is_open() ) {
		printf( "File not found: %s\n", filename );
		exit(1);
	}

	while ( getline(f_data_sequences, seqLine) ) {
		if ( seqLine[0] != '>' ) {
			if ( counter == 0 ) length.push_back(seqLine.size());
			else length.back() += seqLine.size();
			totalLength += seqLine.size();
			counter ++;
		} else {
			if ( counter != 0 ) {
				printf( "Length[%lu]: %lu\n", sequences, length[sequences] );
				sequences ++;
			}
			counter = 0;
			if ( sequences == 3 ) break;
			printf( "%s\n", seqLine.c_str() );
		}
	}

	f_data_sequences.close();
	
	printf( "----------------------------------------------------------------------------\n" );
	printf( "Total Length: %lu\n", totalLength );
}


int main( void ) {
	char *filenameI = "/mnt/ssd0/semin/bancroft/data/sequences/hg002_rep1.fasta";
	
	// Read assembled sequence file 
	fastaReader( filenameI );

	return 0;
}
