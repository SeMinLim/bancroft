#include <sys/time.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <map>
#include <algorithm>
using namespace std;


#define FASTA 0


uint64_t seqSizeOrg = 0;;


// Required Functions
// Time Checker
static inline double timeChecker( void ) {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)(tv.tv_sec) + (double)(tv.tv_usec) / 1000000;
}
// Sequence FASTA File Reader / Writer
void seqFASTA( char *filenameI, char *filenameO ) {
	string seqLine;
	// Read / Write
	ifstream f_data_sequenceI(filenameI);
	ofstream f_data_sequenceO(filenameO);
	while ( getline(f_data_sequenceI, seqLine) ) {
		if ( seqLine[0] != '>' ) {
			f_data_sequenceO << seqLine << "\n";
			seqSizeOrg += seqLine.size();
			if ( seqSizeOrg >= 2948627755 ) break;
		} else {
			f_data_sequenceO << seqLine << "\n";
		}
	}
	// Terminate
	f_data_sequenceI.close();
	f_data_sequenceO.close();
}
// Sequence FASTQ File Reader / Writer
void seqFASTQ( char *filenameI, char *filenameO ) {
	uint64_t counter = 0;
	string seqLine;
	// Read / Write
	ifstream f_data_sequenceI(filenameI);
	ofstream f_data_sequenceO(filenameO);
	while ( getline(f_data_sequenceI, seqLine) ) {
		if ( counter == 0 ) {
			f_data_sequenceO << seqLine << "\n";
			counter ++;
		} else if ( counter == 1 ) {
			f_data_sequenceO << seqLine << "\n";
			seqSizeOrg += seqLine.size();
			counter ++;
		} else if ( counter == 2 ) {
			f_data_sequenceO << seqLine << "\n";
			counter ++;
		} else if ( counter == 3 ) {
			f_data_sequenceO << seqLine << "\n";
			counter = 0;
			if ( seqSizeOrg >= 2948627755 ) break;
		}
	}
	// Terminate
	f_data_sequenceI.close();
	f_data_sequenceO.close();
}


int main( int argc, char **argv ) {
	char *filenameI = argv[1];
	char *filenameO = argv[2];
	
	if ( FASTA ) seqFASTA( filenameI, filenameO );
	else seqFASTQ( filenameI, filenameO );

	return 0;
}
