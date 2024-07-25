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


vector<string> annotation;
vector<string> chromosomes;
vector<uint64_t> length;


void fastaReader( char *filename ) {
	string seqLine;

	ifstream f_data_sequences(filename);
	if ( !f_data_sequences.is_open() ) {
		printf( "File not found: %s\n", filename );
		exit(1);
	}

	while ( getline(f_data_sequences, seqLine) ) {
		if ( seqLine[0] != '>' ) {
			chromosomes.push_back(seqLine);
			length.push_back(seqLine.size());
		} else {
			annotation.push_back(seqLine);
		}
	}

	f_data_sequences.close();
}

void fastaWriter( char *filename ) {
	ofstream f_data_result(filename);
	if ( !f_data_result.is_open() ) {
		printf( "File not found: %s\n", filename );
		exit(1);
	}

	for ( size_t i = 0; i < length.size(); i ++ ) {
		f_data_result << annotation[i] << " : " << length[i] << "\n";
		f_data_result << chromosomes[i] << "\n";
	}
	f_data_result << endl;

	f_data_result.close();
}


int main( void ) {
	char *filenameI = "/home/semin/hg16_tmp.fasta";
	char *filenameO = "/home/semin/hg16.fasta";

	// Read assembled sequence file 
	fastaReader( filenameI );
	printf( "Read Done!\n" );
	fflush( stdout );

	// Attach the each chromosome's length to each annotation
	fastaWriter( filenameO );
	printf( "Write Done!\n" );
	fflush( stdout );

	return 0;
}
