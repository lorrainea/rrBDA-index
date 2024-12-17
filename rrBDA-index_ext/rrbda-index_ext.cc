/**
    rrBDA-index_II_ext: Randomized Reduced Bi-directional Anchors (using external memory)
    Copyright (C) 2024 Lorraine A. K. Ayad, Grigorios Loukides, Solon P. Pissis

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

#include "stream.h"
#include "uint40.h"
#include <math.h>
#include "rrbda-index_ext.h"
#include <sdsl/bit_vectors.hpp>                                   // include header for bit vectors
#include <sdsl/rmq_support.hpp>	
#include <sdsl/io.hpp>
#include "krfp.h"

#ifdef _USE_64
#include <divsufsort64.h>                                         // include header for suffix sort
#endif

#ifdef _USE_32
#include <divsufsort.h>                                       	  // include header for suffix sort
#endif

int main(int argc, char **argv)
{
	unordered_set<unsigned char> alphabet;

	if( argc < 8 )
 	{
        	cout<<"Wrong arguments!\n";
 		cout<<"./rrbda-index_ext <text_file> <ell> <pattern_file> <block_size> <ram_use> <output_filename> <index_filename>\n";
 		exit(-1);
 	}
	
	// Input text file
 	ifstream is_text;
 	char * arg1 = argv[1];
 	char * arg0 = argv[0];
 	
 	is_text.open (arg1, ios::in | ios::binary);
 	
 	ifstream in_file(arg1, ios::binary);
   	in_file.seekg(0, ios::end);
   	INT text_file_size = in_file.tellg();

	// Input ell
 	std::string str_ell(argv[2]);
 	
 	INT ell;
 	std::stringstream(str_ell)>>ell;
	
	// Input block size
 	std::string str_block(argv[4]);
 	
 	INT block;
 	std::stringstream(str_block)>>block;
 	
 	// Input ram use
 	std::string str_ram(argv[5]);
 	
 	INT ram_use;
 	std::stringstream(str_ram)>>ram_use;
 	
 	// Input patterns
 	ifstream is_patterns;
 	is_patterns.open (argv[3], ios::in | ios::binary);
 	
 	// Input output file
 	string output_filename = argv[6]; 
 	
 	// Input index file
 	string index_name = argv[7];
    	ifstream is_index;
 	is_index.open (argv[7], ios::in | ios::binary);

 	unordered_set<INT> text_anchors;
   	
   	unsigned char c = 0;
  	INT text_size = 0;
	for (INT i = 0; i < text_file_size; i++)
	{	
		is_text.read(reinterpret_cast<char*>(&c), 1);
		
		alphabet.insert( (unsigned char) c );
		text_size++;
	}
	is_text.close();
	
	if( text_size < ell )
	{
	
		fprintf( stderr, " Error: Window size (ell) cannot be larger than sequence length!\n");
		return ( 1 );
	}
	
	if( block < ell )
	{
		fprintf( stderr, " Error: Window size (ell) cannot be larger than the block size!\n");
		return ( 1 );
	}
	
	INT k  = ceil(4*log2(ell)/log2(alphabet.size()));
	if( ell - k - 1 < 0 )
		k = 2;
	
	uint64_t hash = karp_rabin_hashing::init();
	uint64_t power = karp_rabin_hashing::pow_mod_mersenne(hash, k, 61);

	/* Compute bd-anchors */
	std::chrono::steady_clock::time_point  start_bd = std::chrono::steady_clock::now();

    	compute_anchors(arg1, text_anchors, text_size, block, ell, k, power );
    	
    	INT g = text_anchors.size();
    	INT n = text_size;
    	
    	std::chrono::steady_clock::time_point  end_bd = std::chrono::steady_clock::now();
	std::cout <<"bd construction took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_bd - start_bd).count() << " [ms]" << std::endl;
	cout<<"The text is of length "<< n << ", its alphabet size is "<< alphabet.size()<<", and it has "<<g<<" bd-anchors of order "<<ell<<endl;
	cout<<"The density is "<<(double) g / n<<endl;
    	
	unsigned char * text_string = ( unsigned char * ) malloc (  ( text_size + 1 ) * sizeof ( unsigned char ) );
	ifstream is_full;
 	is_full.open (arg1, ios::in | ios::binary);
 	
	c = 0;
	for (INT i = 0; i < text_size; i++)
	{	
		is_full.read(reinterpret_cast<char*>(&c), 1);
		text_string[i] = (unsigned char) c;
	}
	is_full.close();
	
	std::chrono::steady_clock::time_point  start_index = std::chrono::steady_clock::now();
	
	INT * RSA;
	INT * RLCP;

	RSA = ( INT * ) malloc( ( g ) * sizeof( INT ) );
	if( ( RSA == NULL) )
	{
	 	fprintf(stderr, " Error: Cannot allocate memory for RSA.\n" );
		return ( 0 );
	}

	RLCP = ( INT * ) malloc( ( g+1 ) * sizeof( INT ) );
	if( ( RLCP == NULL) )
	{
	 	fprintf(stderr, " Error: Cannot allocate memory for RLCP.\n" );
		return ( 0 );
	}
	
	INT * LSA;
  	INT * LLCP;

  	LSA = ( INT * ) malloc( ( g ) * sizeof( INT ) );
  	if( ( LSA == NULL) )
  	{
  		fprintf(stderr, " Error: Cannot allocate memory for LSA.\n" );
        	return ( 0 );
  	}
  	LLCP = ( INT * ) malloc( ( g+1 ) * sizeof( INT ) );
  	if( ( LLCP == NULL) )
  	{
  		fprintf(stderr, " Error: Cannot allocate memory for LLCP.\n" );
        	return ( 0 );
  	}
  	 
  	rmq_succinct_sct<> lrmq;
  	rmq_succinct_sct<> rrmq;
  	
  	compute_index( hash, index_name, text_size, g, text_string, text_anchors, RSA, RLCP, LSA, LLCP, lrmq, rrmq, arg0, arg1, ram_use );
	std::chrono::steady_clock::time_point  end_index = std::chrono::steady_clock::now();
	std::cout <<"Index construction took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_index- start_index + end_bd - start_bd).count() << " [ms]" << std::endl;
  
  	/* Query */
	std::chrono::steady_clock::time_point  begin_pt = std::chrono::steady_clock::now();

	INT hits = query(argv[3], text_string, output_filename, text_size, LSA, LLCP, RSA, RLCP, lrmq, rrmq, g, ell, power, k );
	
 	std::chrono::steady_clock::time_point  end_pt = std::chrono::steady_clock::now();
	std::cout <<"Pattern matching took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_pt - begin_pt).count() << " [ms]" << std::endl;
	std::cout <<"Occurrences: "<< hits <<endl;
 
	free ( RSA );
  	free ( RLCP );
  	free ( LSA );
  	free ( LLCP );
  	free( text_string );
	return 0;
  	
  	return 0;
}
