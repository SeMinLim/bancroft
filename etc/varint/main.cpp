#include <stdio.h>
#include <stdint.h>


#define NUMKMER 16333049


uint64_t totalBytes = 0;


uint64_t varint_encode( uint64_t value, unsigned char *buf ) {
	uint64_t encoded_bytes = 0;
	do {
		uint8_t byte = value & 0x7F; 	// Get the lower 7 bits
		value >>= 7; 			// Shift right by 7 bits
		if (value > 0) byte |= 0x80; 	// Set the MSB if more bytes follow
		buf[encoded_bytes++] = byte;
	} while (value > 0);
	
	return encoded_bytes;
}

uint64_t varint_decode( const unsigned char *buf, int *decoded_bytes ) {
	uint64_t value = 0;
	uint8_t byte = 0;
	int shift = 0;
	int i = 0;
	do {
		byte = buf[i++];
		value |= (uint64_t)(byte & 0x7F) << shift;
		shift += 7;
	} while (byte & 0x80);
	
	*decoded_bytes = i;
	
	return value;
}


int main() {
	// Index
	for ( uint64_t i = 0; i < NUMKMER; i ++ ) {
		unsigned char buf[10];
		uint64_t encoded_bytes = varint_encode(i, buf);
		if ( i + 1 == NUMKMER ) printf( "%ld\n", encoded_bytes );
		totalBytes += encoded_bytes;
	}

	printf( "Encoded Size: %ld\n", totalBytes );
	
	//int decoded_bytes = 0;
	//uint64_t decoded_num = varint_decode(buf, &decoded_bytes);
	//printf( "Decoded number: %lu\n", decoded_num );
	//printf( "Decoded Bytes: %d\n", decoded_bytes );
	
	return 0;
}
