#include <sys/resource.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <boost/algorithm/string.hpp>
using namespace std;


vector<string> sequences;


void sequenceOrganizer( char *filename ) {
	string seqLine;

	sequences.clear();
	sequences.shrink_to_fit();
	seqLine.clear();
	seqLine.shrink_to_fit();

	ifstream f_data_sequences(filename);
	if ( !f_data_sequences.is_open() ) {
		printf( "File not found: %s\n", filename );
		exit(1);
	}

	while ( getline(f_data_sequences, seqLine) ) {
		sequences.push_back(seqLine);
		seqLine.clear();
		seqLine.shrink_to_fit();
	}

	f_data_sequences.close();
}

void fastaMaker( char *filename ) {
	ofstream f_data_result(filename);
	if ( !f_data_result.is_open() ) {
		printf( "File not found: %s\n", filename );
		exit(1);
	}

	int flag = 0;
	for ( size_t i = 0; i < sequences.size(); i ++ ) {
		if ( sequences[i][0] == '>' ) {
			if ( flag == 1 ) f_data_result << "\n";
			f_data_result << sequences[i] << "\n";
		} else {
			boost::to_upper(sequences[i]);
			f_data_result << sequences[i];
			flag = 1;
		}
	}
	f_data_result << endl;

	f_data_result.close();
}


int main() {
	char *filenameIn = "../../data/references/grch38.fasta";
	char *filenameOut = "../../data/references/hg38.fasta";

	printf( "Sequence organizing is started!\n" );
	fflush( stdout );
	sequenceOrganizer( filenameIn );
	printf( "Sequence organizing is finished!\n" );
	fflush( stdout );

	printf( "Writing the organized sequence is started!\n" );
	fflush( stdout );
	fastaMaker( filenameOut );
	printf( "Writing the organized sequence is finished!\n" );
	fflush( stdout );

	return 0;
}
