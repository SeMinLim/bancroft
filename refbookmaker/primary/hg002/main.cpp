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


vector<string> sequences;


uint64_t seqSizeOrg = 0;


void seqAligner( char *filenameI, char *filenameO ) {
	// Sequence Aligner
	uint64_t flag = 0;
	string seqLine;
	ifstream f_data_in(filenameI);
	ofstream f_data_out(filenameO);
	while ( getline(f_data_in, seqLine) ) {
		if ( seqLine[0] != '>' ) {
			f_data_out << seqLine;
			flag = 1;
		} else {
			if ( flag == 1 ) f_data_out << "\n";
			f_data_out << seqLine << "\n";
		}
	}
	f_data_in.close();
	f_data_out.close();
	printf( "----------------------------------------------------------------------------\n" );
	printf( "[STEP 1] Sequence aligning is done!\n" );
	fflush( stdout );
	// Sanity Check
	ifstream f_data_sequence(filenameO);
	while ( getline(f_data_sequence, seqLine) ) {
		if ( seqLine[0] != '>' ) {
			sequences.push_back(seqLine);
			seqSizeOrg += seqLine.size();
		}
	}
	f_data_sequence.close();
	printf( "----------------------------------------------------------------------------\n" );
	printf( "[STEP 2] The Number of Long-Read Sequences: %lu\n", sequences.size() );
	printf( "[STEP 2] The Length of Long-Read Sequences: %lu\n", seqSizeOrg );
	printf( "[STEP 2] Sanity checking is done!\n" );
	printf( "----------------------------------------------------------------------------\n" );
	fflush( stdout );
}


int main( void ) {
	char *filenameI = "/mnt/ephemeral/hg002_rep1.fasta";
	char *filenameO = "/mnt/ephemeral/hg002_rep1_Modified.fasta";
	
	// Read assembled sequence file 
	seqAligner( filenameI, filenameO );

	return 0;
}
