#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <boost/multiprecision/cpp_int.hpp>
using namespace std;
using namespace boost::multiprecision;


#define KMERLENGTH 256
#define ENCKMERBUFUNIT 32
#define ENCKMERBUFSIZE 8
#define BINARYRWUNIT 8


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


int main( void ) {
	string seqLine;
	uint64_t cnt = 0;
	
	char *filenameI = "/mnt/smartssd0/semin/hg19RefBook256Mers_2.txt";
	char *filenameO = "/mnt/smartssd0/semin/hg19RefBook256Mers_2.bin";
	
	string verificationX;
	string verificationY;

	ofstream f_output(filenameO, ios::binary);

	// 1st Reference Book File
	ifstream f_input(filenameI);
	while ( getline(f_input, seqLine) ) {
		if ( cnt == 0 ) verificationX = seqLine.substr(0, 256);
		uint64_t encKmer[ENCKMERBUFSIZE];
		encoder(seqLine, encKmer);
		for ( size_t i = 0; i < ENCKMERBUFSIZE; i ++ ) {
			f_output.write(reinterpret_cast<char *>(&encKmer[i]), BINARYRWUNIT);
		}
		cnt ++;
	}
	f_input.close();
	f_output.close();
	printf( "Writing the Reference Book File as Binary is Done\n" );
	fflush( stdout );

	// Verification
	uint64_t verificationEncoded[ENCKMERBUFSIZE];

	ifstream f_verification(filenameO, ios::binary);
	for ( size_t i = 0; i < ENCKMERBUFSIZE; i ++ ) {
		f_verification.read(reinterpret_cast<char *>(&verificationEncoded[i]), BINARYRWUNIT);
	}
	decoder(verificationEncoded, verificationY);
	f_verification.close();
	cout << verificationX << "\n";
	cout << verificationY << "\n";
	cout << verificationX.compare(verificationY) << "\n";

	return 0;
}
