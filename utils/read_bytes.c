/***********************************************************************
* Shimon Panfil: Industrial Physics and Simulations                   **
* http://industrialphys.com                                           **
* THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           **
***********************************************************************/
#include "my_cdefs.h"
char * read_bytes(char* filename) {
    FILE* infile;	
    int len;
    char* source;
#ifdef WIN32
    fopen_s( &infile, filename, "rb" );
#else
    infile = fopen( filename, "rb" );
#endif // WIN32
    if ( !infile ) {
	fprintf(stderr,"can't open file %s\n",filename); 
	return NULL;
    }
    fseek( infile, 0, SEEK_END );
    len = ftell( infile );
    fseek( infile, 0, SEEK_SET );
    source = ALLOC_MEM(char,len+1);
    if(fread( source, 1, len, infile )!=len) {
	fprintf(stderr,"error reading file %s\n",filename);
	FREE_MEM(source);
	return NULL;
    }
    fclose( infile );
    source[len] = 0;
    return source;
}

