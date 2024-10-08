/*
Copyright (C) 2017 Instituto Superior Tecnico

This file is part of the ZPIC Educational code suite

The ZPIC Educational code suite is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

The ZPIC Educational code suite is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with the ZPIC Educational code suite. If not, see <http://www.gnu.org/licenses/>.
*/


/**
 * ZDF version 1
 *
 * This ZDF version is totally self contained. It does not depend on XDR.
 * Current implementation should also work on big endian systems (untested)
 * This version is not compatible is version 0 (it has the opposite endianess)
 * 
 */


#include "zdf.h"


#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string.h>
#include <errno.h>

/**
 * On Windows we cannot use the POSIX.1 mkdir command, so use _mkdir instead
 */
#if defined(_MSC_VER) || defined(_WIN32) || defined(_WIN64)
#include <direct.h>
#define mkdir(path,mode) _mkdir(path)
#endif

/**
 * Number of bytes in each ZDF file unit
 */
#define BYTES_PER_ZDF_UNIT	(4)


/**
 * Number of bytes required for writing data to ZDF file
 * (round up to multiple of BYTES_PER_ZDF_UNIT)
 */
#define RNDUP(x)  ((((x) + BYTES_PER_ZDF_UNIT - 1) / BYTES_PER_ZDF_UNIT) \
		    * BYTES_PER_ZDF_UNIT)


/**
 * Magic byte sequence identifying ZDF files
 */
#define ZDF_MAGIC_LENGTH BYTES_PER_ZDF_UNIT
const char zdf_magic[ZDF_MAGIC_LENGTH] = {'Z','D','F','1'};

/**
 * ZDF datatypes - Currently only float32 and float64 are supported
 */
enum zdf_data_type{ zdf_null,
	                zdf_int8,  zdf_uint8,  zdf_int16, zdf_uint16, 
	                zdf_int32, zdf_uint32, zdf_int64, zdf_uint64,
	                zdf_float32, zdf_float64 };

/**
 * Buffer size for converting big endian to little endian data
 */
#define ENDIAN_CONV_BUF_SIZE 1024

/**
 * Sizes of datatypes
 */
const unsigned size_zdf_int32   = 4;
const unsigned size_zdf_uint32  = 4;
const unsigned size_zdf_uint64  = 8;
const unsigned size_zdf_double  = 8;
const unsigned size_zdf_float   = 4;


/**
 * IDs of ZDF records
 */
#define ZDF_INT32_ID     0x00010000
#define ZDF_DOUBLE_ID    0x00020000
#define ZDF_STRING_ID    0x00030000

#define ZDF_DATASET_ID   0x00100000

#define ZDF_ITERATION_ID 0x00200000
#define ZDF_GRID_INFO_ID 0x00210000
#define ZDF_PART_INFO_ID 0x00220000

/* -----------------------------------------------------------------------------------------------
  recursively create path if required
-------------------------------------------------------------------------------------------------- */

/**
 * Recursively creates a path if required
 * @param  path path to create
 * @return      returns 0 on success, otherwise returns the value returned by mkdir
 */
int create_path( const char path[] )
{
	char uppath[256], *p;
	int ierr = 0;
	
	if (mkdir(path,S_IRWXU | (S_IRGRP | S_IXGRP ) | (S_IROTH | S_IXOTH) )) {
		switch (errno) {
			case ENOENT : // A component of the path does not exist 
				
				// get upper path
				strncpy(uppath, path, 256);
				p = uppath + strlen(uppath);
				while(*p!='/') p--;
				*p=0;
				
				//recursively build the path
				if ( !create_path( uppath ) ) ierr = create_path( path );
				break;
				
			case EEXIST : /* if directory already exists ignore the error */
				ierr = 0;
				break;
			default: ierr = errno;
		}
	}
	return ierr;
}

/* -----------------------------------------------------------------------------------------------
  Open / Close ZDF file
-------------------------------------------------------------------------------------------------- */

/**
 * Closes ZDF file
 * @param  zdf ZDF file to close
 * @return     Returns 0 on success, -1 otherwise
 */
int zdf_close_file( t_zdf_file* zdf ) {

	if ( fclose( zdf->fp ) ) {
		perror("(*error*) Unable to close ZDF file");
		return(-1);
	}
	zdf -> fp = NULL;
	return(0);
}

/**
 * Opens ZDF file
 * @param  zdf      ZDF file to open
 * @param  filename Filename of the ZDF file to open, including path
 * @param  mode     Can be one of ZDF_WRITE (for writing) or ZDF_READ (for reading)
 * @return          Returns 0 on success, -1 otherwise
 */
int zdf_open_file( t_zdf_file* zdf, char* filename, enum zdf_file_access_mode mode ){

	zdf -> mode = mode;

	switch( mode ) {
		case ZDF_WRITE :
			// Open file for writing
			// The "wb" mode must be used for compatibility with Windows
			if (!(zdf->fp = fopen( filename, "wb"))) {
				perror("(*error*) Unable to open ZDF file for writing");
				return(-1);
			}
			
			// Write magic number
			if (!fwrite( (void *) &zdf_magic, sizeof(char), 4, zdf->fp )) {
				fprintf(stderr, "(*error*) Unable to write magic number to ZDF file.");
				zdf_close_file( zdf );
				return(-1);
			}

			break;

		case ZDF_READ :
			// Open file for reading
			if (!(zdf->fp = fopen( filename, "r"))) {
				perror("(*error*) Unable to open ZDF file for reading");
				return(-1);
			}

			// Read magic number
			char test_magic[4];
			if (! fread( (void *) &test_magic, sizeof(char), 4, zdf->fp )) {
				fprintf(stderr, "(*error*) Unable to read magic number from ZDF file.");
				zdf_close_file( zdf );
				return(-1);
			}

            // Check magic number
            for( int i = 0; i < ZDF_MAGIC_LENGTH; i++) {
	            if ( test_magic[i] != zdf_magic[i] ) {
					fprintf(stderr, "(*error*) Invalid magic number, file is not a proper ZDF file.");
					zdf_close_file( zdf );
					return(-1);
				}
			}
			break;

		default:
			fprintf(stderr, "(*error*) zdf_open_file: unsupported mode\n");
			return(-1);
	}

	return(0);
}

/* -----------------------------------------------------------------------------------------------
  Elemental types
-------------------------------------------------------------------------------------------------- */


#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__

/**
 * Implementation for little endian systems (e.g. x86)
 *
 * For these systems just writing the data to disk is sufficient
 */

/**
 * Writes int32 value to file
 * @param  zdf ZDF file descriptor
 * @param  i   int32_t value to write
 * @return     Returns 0 on success, other value on error
 */
int zdf_int32_write( t_zdf_file* zdf, const int32_t i ){
	return ( fwrite( (void *) &i, sizeof(int32_t), 1, zdf -> fp ) == 1 );
}

/**
 * Writes uint32 value to file
 * @param  zdf ZDF file descriptor
 * @param  u   uint32_t value to write
 * @return     Returns 0 on success, other value on error
 */
int zdf_uint32_write( t_zdf_file* zdf, const uint32_t u ){
	return ( fwrite( (void *) &u, sizeof(uint32_t), 1, zdf -> fp ) == 1 );
}

/**
 * Writes int64 value to file
 * @param  zdf ZDF file descriptor
 * @param  i   int64_t value to write
 * @return     Returns 0 on success, other value on error
 */
int zdf_int64_write( t_zdf_file* zdf, const int64_t i ){
	return ( fwrite( (void *) &i, sizeof(int64_t), 1, zdf -> fp ) == 1 );
}

/**
 * Writes uint64_t value to file
 * @param  zdf ZDF file descriptor
 * @param  u   uint64_t value to write
 * @return     Returns 0 on success, other value on error
 */
int zdf_uint64_write( t_zdf_file* zdf, const uint64_t u ){
	return ( fwrite( (void *) &u, sizeof(uint64_t), 1, zdf -> fp ) == 1 );
}

/**
 * Writes double (float64) value to file
 * @param  zdf ZDF file descriptor
 * @param  d   double value to write
 * @return     Returns 0 on success, other value on error
 */
int zdf_double_write( t_zdf_file* zdf, const double d ){
	return ( fwrite( (void *) &d, sizeof(double), 1, zdf -> fp ) == 1 );
}

/**
 * Write float (float32) vector to file
 * @param  zdf  ZDF file descriptor
 * @param  data Pointer to float (float32) data to write
 * @param  len  Number of vector elements
 * @return      Returns 0 on success, other value on error
 */
int zdf_float_vector_write( t_zdf_file* zdf,  float const * const data, size_t len ) {
	return( fwrite( (void *) data, sizeof(float), len, zdf -> fp ) != len );
}

/**
 * Write double (float64) vector to file
 * @param  zdf  ZDF file descriptor
 * @param  data Pointer to float (float64) data to write
 * @param  len  Number of vector elements
 * @return      Returns 0 on success, other value on error
 */
int zdf_double_vector_write( t_zdf_file* zdf,  double const * const data, size_t len ) {
	return( fwrite( (void *) data, sizeof(double), len, zdf -> fp ) != len );
}

#elif __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__

/**
 * Implementation for big endian systems (e.g. PowerPC)
 *
 * For these systems we need to swap the bytes before writing. For vectors, conversion
 * is done in chunks of ENDIAN_CONV_BUF_SIZE values.
 *
 * The system may provide hardware optimized bswap_* routines, these are usually found in
 * the <byteswap.h> header, but they are missing on many platforms.
 */


/*
// Not required
#define bswap_16(x) \
({ \
	uint16_t __x = (x); \
	((uint16_t)( \
		(((uint16_t)(__x) & (uint16_t)0x00ffU) << 8) | \
		(((uint16_t)(__x) & (uint16_t)0xff00U) >> 8) )); \
})
*/

/**
 * Swap bytes for 32 bit number
 * @param  x 32 bit number to convert
 * @return   byte-swapped version of x
 */
#define bswap_32(x) \
({ \
	uint32_t __x = (x); \
	((uint32_t)( \
		(((uint32_t)(__x) & (uint32_t)0x000000ffUL) << 24) | \
		(((uint32_t)(__x) & (uint32_t)0x0000ff00UL) <<  8) | \
		(((uint32_t)(__x) & (uint32_t)0x00ff0000UL) >>  8) | \
		(((uint32_t)(__x) & (uint32_t)0xff000000UL) >> 24) )); \
})

/**
 * Swap bytes for 64 bit number
 * @param  x 64 bit number to convert
 * @return   byte-swapped version of x
 */
#define bswap_64(x) \
({ \
	uint64_t __x = (x); \
	((uint64_t)( \
		(uint64_t)(((uint64_t)(__x) & (uint64_t)0x00000000000000ffULL) << 56) | \
		(uint64_t)(((uint64_t)(__x) & (uint64_t)0x000000000000ff00ULL) << 40) | \
		(uint64_t)(((uint64_t)(__x) & (uint64_t)0x0000000000ff0000ULL) << 24) | \
		(uint64_t)(((uint64_t)(__x) & (uint64_t)0x00000000ff000000ULL) <<  8) | \
	    (uint64_t)(((uint64_t)(__x) & (uint64_t)0x000000ff00000000ULL) >>  8) | \
		(uint64_t)(((uint64_t)(__x) & (uint64_t)0x0000ff0000000000ULL) >> 24) | \
		(uint64_t)(((uint64_t)(__x) & (uint64_t)0x00ff000000000000ULL) >> 40) | \
		(uint64_t)(((uint64_t)(__x) & (uint64_t)0xff00000000000000ULL) >> 56) )); \
})

/**
 * Writes int32 value to file
 * @param  zdf ZDF file descriptor
 * @param  i   int32_t value to write
 * @return     Returns 0 on success, other value on error
 */
int zdf_int32_write( t_zdf_file* zdf, const int32_t i ){
	uint32_t tmp = bswap_32((uint32_t)i);
	return ( fwrite( (void *) &tmp, sizeof(uint32_t), 1, zdf -> fp ) == 1 );
}

/**
 * Writes uint32 value to file
 * @param  zdf ZDF file descriptor
 * @param  u   uint32_t value to write
 * @return     Returns 0 on success, other value on error
 */
int zdf_uint32_write( t_zdf_file* zdf, const uint32_t u ){
	uint32_t tmp = bswap_32(u);
	return ( fwrite( (void *) &tmp, sizeof(uint32_t), 1, zdf -> fp ) == 1 );
}

/**
 * Writes int64 value to file
 * @param  zdf ZDF file descriptor
 * @param  i   int64_t value to write
 * @return     Returns 0 on success, other value on error
 */
int zdf_int64_write( t_zdf_file* zdf, const int64_t i ){
	uint64_t tmp = bswap_64((uint64_t)i);
	return ( fwrite( (void *) &tmp, sizeof(uint64_t), 1, zdf -> fp ) == 1 );
}

/**
 * Writes uint64_t value to file
 * @param  zdf ZDF file descriptor
 * @param  u   uint64_t value to write
 * @return     Returns 0 on success, other value on error
 */
int zdf_uint64_write( t_zdf_file* zdf, const uint64_t u ){
	uint64_t tmp = bswap_64(u);
	return ( fwrite( (void *) &tmp, sizeof(uint64_t), 1, zdf -> fp ) == 1 );
}

/**
 * Writes double (float64) value to file
 * @param  zdf ZDF file descriptor
 * @param  d   double value to write
 * @return     Returns 0 on success, other value on error
 */
int zdf_double_write( t_zdf_file* zdf, const double d ){
	uint64_t tmp = bswap_64((uint64_t)d);
	return ( fwrite( (void *) &tmp, sizeof(uint64_t), 1, zdf -> fp ) == 1 );
}

/**
 * Write float (float32) vector to file
 * @param  zdf  ZDF file descriptor
 * @param  data Pointer to float (float32) data to write
 * @param  len  Number of vector elements
 * @return      Returns 0 on success, other value on error
 */
int zdf_float_vector_write( t_zdf_file* zdf,  float const * const data, size_t len ) {

	uint32_t buffer[ENDIAN_CONV_BUF_SIZE];

	for( size_t offset = 0; offset < len; offset += ENDIAN_CONV_BUF_SIZE ) {

		// Number of values in chunk
		size_t chunk_len = (offset + ENDIAN_CONV_BUF_SIZE < len ) ? ENDIAN_CONV_BUF_SIZE : len - offset;
		
		// Convert chunk to little endian
		for( size_t i = 0; i < chunk_len; i++) buffer[i] = bswap_32( (uint32_t) data[offset+i]);
		
		// Write chunk
		if ( fwrite( (void *) buffer, sizeof(uint32_t), chunk_len, zdf -> fp) != chunk_len ) 
			return(0);
    }

	return( 1 );
}

/**
 * Write double (float64) vector to file
 * @param  zdf  ZDF file descriptor
 * @param  data Pointer to float (float64) data to write
 * @param  len  Number of vector elements
 * @return      Returns 0 on success, other value on error
 */
int zdf_double_vector_write( t_zdf_file* zdf,  double const * const data, size_t len ) {
	uint64_t buffer[ENDIAN_CONV_BUF_SIZE];

	for( size_t offset = 0; offset < len; offset += ENDIAN_CONV_BUF_SIZE ) {

		// Number of values in chunk
		size_t chunk_len = (offset + ENDIAN_CONV_BUF_SIZE < len ) ? ENDIAN_CONV_BUF_SIZE : len - offset;
		
		// Convert chunk to little endian
		for( size_t i = 0; i < chunk_len; i++) buffer[i] = bswap_64( (uint64_t) data[offset+i]);
		
		// Write chunk
		if ( fwrite( (void *) buffer, sizeof(uint64_t), chunk_len, zdf -> fp) != chunk_len ) 
			return(0);
    }

	return( 1 );
}


#else

#error "System is neither little endian nor big endian, aborting."

#endif

/**
 * Write a sequence of bytes to file
 * @param  zdf ZDF file descriptor
 * @param  u   Pointer to byte data
 * @param  len Number of bytes to write
 * @return     Returns 0 on success, other value on error
 */
int zdf_bytes_write( t_zdf_file* zdf, const uint8_t *u, size_t len ){
	
	if ( fwrite( (void *) u, sizeof(uint8_t), len, zdf -> fp ) != len ) {
		return(0);
	}

	size_t npad = RNDUP(len) - len;

	if ( npad > 0 ) {
		const uint8_t pad[BYTES_PER_ZDF_UNIT] = {0};
		if ( fwrite( (void *) pad, sizeof(uint8_t), npad, zdf -> fp ) != npad ) {
			return(0);
		}
	}

	return( RNDUP(len) );
}


/* -----------------------------------------------------------------------------------------------
  zdf_string
-------------------------------------------------------------------------------------------------- */

/**
 * Write a string to file
 * @param  zdf ZDF file descriptor
 * @param  str C string (char *) to write
 * @return     Returns 0 on success, 1 on error
 */
int zdf_string_write( t_zdf_file* zdf, const char * str ){

	uint32_t len;

	len = ( str ) ? strlen( str ) : 0;
	if ( !zdf_uint32_write( zdf, len ) ) return(0);

	if ( len > 0 ) {
		if ( !zdf_bytes_write( zdf, (const uint8_t *) str, len ) ) return(0);
	}

	return(1);
}

/**
 * Calculate size of string record
 * @param  s C string (char *) to write
 * @return   The number of bytes required for writing the string data
 */
uint64_t size_zdf_string( const char *s )
{
	unsigned len = ( s ) ? strlen( s ) : 0;
	return size_zdf_uint32 + ((len > 0) ? RNDUP(len) : 0);
}


/* -----------------------------------------------------------------------------------------------
  zdf records
-------------------------------------------------------------------------------------------------- */

typedef struct {
	uint32_t id_version;
	char*    name;
	uint64_t length;
} t_zdf_record;


int zdf_record_write( t_zdf_file* zdf, const t_zdf_record* rec ){

    if ( !zdf_uint32_write( zdf, rec -> id_version ) ) return(0);
 	if ( !zdf_string_write( zdf, rec -> name ) ) return(0);
	if ( !zdf_uint64_write( zdf, rec -> length ) ) return(0);

 	return(1);
}


/* -----------------------------------------------------------------------------------------------
  zdf basic data tags
-------------------------------------------------------------------------------------------------- */
int zdf_add_string( t_zdf_file* zdf, char* name, const char* str ){
	
	t_zdf_record rec;
	rec.id_version = ZDF_STRING_ID;
	rec.name = name;
	rec.length = size_zdf_string( str );
	

    if ( !zdf_record_write( zdf, &rec) ) return(-1);
    if ( !zdf_string_write( zdf, str) ) return(-1);

    return(0);
}

int zdf_add_int32( t_zdf_file* zdf, char* name, const int32_t value ){
	
	t_zdf_record rec = {
		.id_version = ZDF_INT32_ID,
		.name = name,
		.length = size_zdf_int32
	};

    if ( !zdf_record_write(zdf, &rec) ) return(-1);
    if ( !zdf_int32_write(zdf, value) ) return(-1);

    return(0);
}

int zdf_add_double( t_zdf_file* zdf, char* name, const double value )
{
	t_zdf_record rec = {
		.id_version = ZDF_DOUBLE_ID,
		.name = name,
		.length = size_zdf_double
	};

    if ( !zdf_record_write(zdf, &rec) ) return(-1);
    if ( !zdf_double_write(zdf, value) ) return(-1);

    return(0);
}


/* -----------------------------------------------------------------------------------------------
  zdf compound metadata tags
-------------------------------------------------------------------------------------------------- */

int zdf_add_iteration( t_zdf_file* zdf, const char* name, const t_zdf_iteration* iter ){

	t_zdf_record rec;
	rec.id_version = ZDF_ITERATION_ID;
	rec.name = (char *) name;
	rec.length = size_zdf_uint32 +
	             size_zdf_double +
	             size_zdf_string( iter -> time_units );
	

  if ( !zdf_record_write( zdf, &rec) ) return(-1);

 	if ( !zdf_int32_write( zdf, iter ->n ) ) return(-1);
	if ( !zdf_double_write( zdf, iter->t ) ) return(-1);
 	if ( !zdf_string_write( zdf, iter->time_units ) ) return(-1);

 	return(0);
}


uint64_t size_xdr_zdf_grid_info(const t_zdf_grid_info* grid) {
	uint64_t size;

	size = size_zdf_uint32 + grid -> ndims * size_zdf_uint64 +
	        size_zdf_string(grid->label) +  size_zdf_string(grid->units);

   	// Includes axis information
   	size += size_zdf_int32;
    if ( grid -> axis ) {

    	for(unsigned int i=0; i<grid -> ndims; i++) 
    		size += size_zdf_int32 + 
    	            2 * size_zdf_double + 
    	            size_zdf_string( grid -> axis[i].label ) +
    	            size_zdf_string( grid -> axis[i].units );
    }

	return size;
}

int zdf_add_grid_info( t_zdf_file* zdf, const char* name, const t_zdf_grid_info* grid ){

	t_zdf_record rec;
	rec.id_version = ZDF_GRID_INFO_ID;
	rec.name = (char *) name;
	rec.length = size_xdr_zdf_grid_info( grid );
	

    if ( !zdf_record_write( zdf, &rec) ) return(-1);

 	if ( !zdf_uint32_write( zdf, grid -> ndims ) ) return(-1);

    unsigned int i;
    for( i=0; i < grid -> ndims; i++) {
	 	if ( !zdf_uint64_write( zdf, grid -> nx[i] ) ) return(-1);
    }
 	
 	if ( !zdf_string_write( zdf, grid->label ) ) return(-1);
 	if ( !zdf_string_write( zdf, grid->units ) ) return(-1);

    int32_t has_axis = ( grid -> axis != NULL);
 	if ( !zdf_int32_write( zdf, has_axis ) ) return(-1);

 	if ( has_axis ) {

	    for( i=0; i < grid -> ndims; i++) {

		    if ( !zdf_int32_write( zdf, grid -> axis[i].type ) ) return(-1);
		    if ( !zdf_double_write( zdf, grid -> axis[i].min ) ) return(-1);
		    if ( !zdf_double_write( zdf, grid -> axis[i].max ) ) return(-1);

		    if ( !zdf_string_write( zdf, grid -> axis[i].label ) ) return(-1);
		    if ( !zdf_string_write( zdf, grid -> axis[i].units ) ) return(-1);
	    }
 	}

 	return(0);
}


uint64_t size_xdr_zdf_part_info(const t_zdf_part_info* part) {
	uint64_t size = size_zdf_string(part->name) +	// name
			size_zdf_uint32; 								// nquants

	unsigned int i;
	for( i = 0; i < part -> nquants; i++) {
		size += size_zdf_string( part -> quants[i] ) +
				size_zdf_string( part -> units[i] );
	}

	size += size_zdf_uint64; 					// np

	return size;
}

int zdf_add_part_info( t_zdf_file* zdf, char* name, t_zdf_part_info* part ){

	t_zdf_record rec;
	rec.id_version = ZDF_PART_INFO_ID;
	rec.name = name;
	rec.length = size_xdr_zdf_part_info( part );
	

    if ( !zdf_record_write( zdf, &rec) ) return(-1);
	if ( !zdf_string_write( zdf, part->name ) ) return(-1);
 	if ( !zdf_uint32_write( zdf, part -> nquants ) ) return(-1);

    unsigned int i;
    for( i=0; i < part -> nquants; i++) {
		if ( !zdf_string_write( zdf, part->quants[i] ) ) return(-1);
    }
    for( i=0; i < part -> nquants; i++) {
		if ( !zdf_string_write( zdf, part->units[i] ) ) return(-1);
    }
 	
 	if ( !zdf_uint64_write( zdf, part -> np ) ) return(-1);

 	return(0);
}

/* -----------------------------------------------------------------------------------------------
  zdf dataset
-------------------------------------------------------------------------------------------------- */

typedef struct {
	enum zdf_data_type data_type;
	uint32_t ndims;
	uint64_t nx[zdf_max_dims];

	uint8_t* data;
} t_zdf_dataset;


uint64_t zdf_datatype_size( enum zdf_data_type data_type ) {
    uint64_t size;

    switch ( data_type ) {
    	case zdf_float32: size = size_zdf_float; break;
    	case zdf_float64: size = size_zdf_double; break;
    	default : size = 0;
    }

    return(size);
}

uint64_t size_zdf_dataset(const t_zdf_dataset* dataset) {
	
	unsigned int i;
    uint64_t data_size;

    data_size = zdf_datatype_size( dataset -> data_type );
    for( i = 0; i < dataset -> ndims; i++ ) {
    	data_size *= dataset -> nx[i];
    }

	return size_zdf_int32 + size_zdf_uint32 + 
	       dataset -> ndims * size_zdf_uint64 +
	       data_size ;
}

int zdf_add_dataset( t_zdf_file* zdf, char* name, t_zdf_dataset* dataset ){

	t_zdf_record rec;
	rec.id_version = ZDF_DATASET_ID;
	rec.name       = name;
	rec.length     = size_zdf_dataset( dataset );

    if ( !zdf_record_write( zdf, &rec) ) return(-1);

 	if ( !zdf_int32_write( zdf, dataset -> data_type ) ) return(-1);
	if ( !zdf_uint32_write( zdf, dataset -> ndims ) ) return(-1);

    unsigned int i;
 	unsigned int count;
    for( i=0, count = 1; i < dataset -> ndims; i++) {
    	count *= dataset -> nx[i];
	 	if ( !zdf_uint64_write( zdf, dataset -> nx[i] ) ) return(-1);
    }

 	switch ( dataset -> data_type ) {
    	case zdf_float32: 
    	    if ( !zdf_float_vector_write( zdf, (float *) dataset -> data, count ) ) return(-1);
    	    break;
    	case zdf_float64:
    	    if ( !zdf_double_vector_write( zdf, (double *) dataset -> data, count ) ) return(-1);
    	    break;
    	default: 
    		fprintf(stderr,"(*error*) zdf_add_dataset: Unsupported datatype.");
    		return(-1);
    		break;
 	}
 	
 	return(0);
}


/* -----------------------------------------------------------------------------------------------
  zdf high level interface
-------------------------------------------------------------------------------------------------- */

int zdf_save_grid( const float* data, const t_zdf_grid_info *_info, 
	const t_zdf_iteration *_iteration, char const path[] )
{

	unsigned int i;
	char filename[1024];

    // Set iteration info
    t_zdf_iteration iteration;
    iteration.n = _iteration->n;
    iteration.t = _iteration->t;
    iteration.time_units = _iteration->time_units;

    // Set axis info
    t_zdf_grid_axis axis[zdf_max_dims];

    if ( _info -> axis ) {
	    // Copy axis info
	    for( i = 0; i < _info->ndims; i ++) {
	    	axis[i].type = zdf_linear;
	    	axis[i].min = _info -> axis[i].min;
	    	axis[i].max = _info -> axis[i].max;
	    	axis[i].label = _info -> axis[i].label;
	    	axis[i].units = _info -> axis[i].units;
	    }

    } else {
    	// set default axis info
	    for( i = 0; i < _info->ndims; i ++) {
	    	axis[i].type = zdf_linear;
	    	axis[i].min = 0;
	    	axis[i].max = _info -> nx[i];
	    	axis[i].label = "";
	    	axis[i].units = "" ;
	    }
    }

    // Set grid info
    t_zdf_grid_info   grid_info;
    grid_info.ndims = _info->  ndims;
    grid_info.label = _info -> label;
    grid_info.units = _info -> units;
    grid_info.axis  = axis;

    for( i = 0; i < _info->ndims; i ++) grid_info.nx[i] = _info->nx[i];

    // Set data
    t_zdf_dataset dataset;
    dataset.data_type = zdf_float32;
    dataset.ndims = _info->ndims;
    dataset.data = (uint8_t *) data;

    for( i = 0; i < _info->ndims; i ++) dataset.nx[i] = _info->nx[i];

	// Ensure that the path is available
	create_path( path );
	
	// Build filename
	sprintf( filename, "%s/%s-%06u.zdf", path, _info->label, _iteration -> n );
    // printf("Saving filename %s\n", filename );
	
    // Create ZDF file
    t_zdf_file zdf;
    if ( zdf_open_file( &zdf, filename, ZDF_WRITE ) ) {
    	fprintf(stderr,"(*error*) Unable to open ZDF file, aborting.");
    	return(-1);
    }

    // Add file type
    zdf_add_string( &zdf, "TYPE", "grid");

    // Add grid info
    zdf_add_grid_info( &zdf, "GRID", &grid_info );

    // Add iteration info
    zdf_add_iteration( &zdf, "ITERATION", &iteration );

    // Add dataset
    zdf_add_dataset( &zdf, "DATA", &dataset );

    // Close ZDF file and return
    return( zdf_close_file( &zdf ) );
}

int zdf_part_file_open( t_zdf_file *zdf, t_zdf_part_info *_info, 
	const t_zdf_iteration *_iteration, char const path[] ){

	char filename[1024];

    t_zdf_iteration iteration = *_iteration;
    t_zdf_part_info info = *_info;

	// Ensure that the path is available
	create_path( path );
	
	// Build filename
	sprintf( filename, "%s/%s-%s-%06u.zdf", path, "particles", _info->name, _iteration -> n );
    //printf("Saving filename %s\n", filename );
	
    // Create ZDF file
    if ( zdf_open_file( zdf, filename, ZDF_WRITE ) ) {
    	fprintf(stderr,"(*error*) Unable to open ZDF file, aborting.");
    	return(-1);
    }

    // Add file type
    zdf_add_string( zdf, "TYPE", "particles");

    // Add particle info
    zdf_add_part_info( zdf, "PARTICLES", &info );

    // Add iteration info
    zdf_add_iteration( zdf, "ITERATION", &iteration );

    return(0);
}

int zdf_part_file_add_quant( t_zdf_file *zdf, const char *name, const float* data, 
	const unsigned np ) {

    t_zdf_dataset dataset = {
    	.data_type = zdf_float32,
    	.ndims = 1,
    	.data = (uint8_t *) data
    };

    dataset.nx[0] = np;

    zdf_add_dataset( zdf, (char *) name, &dataset );

    return(0);

}

#ifdef __TEST_ZDF__

#include <math.h>

int main (int argc, const char * argv[]) {

	const unsigned NX = 128;
	float buf[NX];

	for(unsigned i = 0; i < NX; i++) {
		float x = 8 * (M_PI/NX) * (i+1);
		buf[i] = sin(x)/x;
	}

    t_zdf_grid_axis axis[1];
    axis[0] = (t_zdf_grid_axis) {
    	.min = -1.0,
    	.max =  1.0,
    	.label = "axis label",
    	.units = "axis units"
    };

    t_zdf_grid_info info = {
    	.ndims = 1,
    	.label = "data label",
    	.units = "data units",
    	.axis = axis
    };

    info.nx[0] = NX;

    t_zdf_iteration iter = {
    	.n = 123,
    	.t = 12.3,
    	.time_units = "time units"
    };

	zdf_save_grid( buf, &info, &iter, "V1TEST" );

	return 0;
}

#endif
