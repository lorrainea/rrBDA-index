/**
    rrBDA-index_II_int: Randomized Reduced Bi-directional Anchors (using internal memory)
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


#include "utils.h"
#include "stream.h"
#include "uint40.h"
#include <math.h>
#include "rrbda-index_int.h"
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

double vm, vm0, rss, rss0;

void reverse( unsigned char * &s)
{
	INT length = strlen( (char*) s ) ;
	
	unsigned char temp = ' ';
	INT i = 0, j = length-1;
	
	while( i < j )
	{
		temp = s[j];
		s[j] = s[i];
		s[i] = temp;
		
		i++, j--;
	}
}

void process_mem_usage(double& vm_usage, double& resident_set)
{
    vm_usage     = 0.0;
    resident_set = 0.0;

    // the two fields we want
    unsigned long vsize;
    long rss;
    {
        std::string ignore;
        std::ifstream ifs("/proc/self/stat", std::ios_base::in);
        ifs >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
                >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
                >> ignore >> ignore >> vsize >> rss;
    }

    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    vm_usage = vsize / 1024.0;
    resident_set = rss * page_size_kb;
}



/* Sorting comparison */
bool sort_sa(const pair<INT,INT> &a,const pair<INT,INT> &b)
{
       return a.first<b.first;

}

int main(int argc, char **argv)
{
	unordered_set<unsigned char> alphabet;

	if( argc < 7 )
 	{
        	cout<<"Wrong arguments!\n";
 		cout<<"./rbda-index_II <text_file> <ell> <pattern_file> <block_size> <output_filename> <index_filename>\n";
 		exit(-1);
 	}
	
	// Input text file
 	ifstream is_text;
 	is_text.open (argv[1], ios::in | ios::binary);
 	
 	ifstream in_file(argv[1], ios::binary);
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
 	
 	// Input patterns
 	ifstream is_patterns;
 	is_patterns.open (argv[3], ios::in | ios::binary);
 	
 	// Input output file
 	string output_filename = argv[5]; 
 	
 	// Input index file
 	string index_name = argv[6];
    	ifstream is_index;
 	is_index.open (argv[6], ios::in | ios::binary);

    	std::chrono::steady_clock::time_point  start_bd = std::chrono::steady_clock::now();
 	unordered_set<INT> text_anchors;
	
   	unsigned char c = 0;
  	INT text_size = 0;
	for (INT i = 0; i < text_file_size; i++)
	{	
		is_text.read(reinterpret_cast<char*>(&c), 1);
		
		if( (unsigned char) c == '\n' )
			continue;
		else
		{
			alphabet.insert( (unsigned char) c );
			text_size++;
		}
		
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
	
	
	INT hash = karp_rabin_hashing::init();
	INT power = karp_rabin_hashing::pow_mod_mersenne(hash, k, 61);
	
   	ifstream is_block;
	is_block.open (argv[1], ios::in | ios::binary);
	unsigned char * text_block = ( unsigned char * ) malloc (  ( block + 1 ) * sizeof ( unsigned char ) );
	INT * rank = ( INT * ) malloc( ( block  ) *  sizeof( INT ) );
	unsigned char * suffix_block = ( unsigned char * ) malloc (  ( ell  ) * sizeof ( unsigned char ) );
		
	c = 0;
	INT count = 0;
	INT pos = 0;
	 	
	for (INT i = 0; i < text_size; i++)
	{	
		is_block.read(reinterpret_cast<char*>(&c), 1);
			
		if( (unsigned char) c != '\n' )
		{
			text_block[count] = (unsigned char) c ;
			count++;
			if( count == block || i == text_size - 1 )
			{
				text_block[count] = '\0';
					
				bd_anchors( text_block, pos, ell, k, text_anchors, rank, power );
					
				memcpy( &suffix_block[0], &text_block[ block - ell + 1], ell -1 );
				memcpy( &text_block[0], &suffix_block[0], ell -1 );
					
				pos = pos + ( block - ell + 1 );
				count = ell - 1;	
			}
		}
	}
		
	is_block.close();
	free( text_block );
	free( suffix_block );
	free( rank );	
		
	INT g = text_anchors.size();
	INT n = text_size;
	
	std::chrono::steady_clock::time_point  end_bd = std::chrono::steady_clock::now();
	std::cout <<"bd construction took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_bd - start_bd).count() << " [ms]" << std::endl;
	cout<<"The text is of length "<< n << ", its alphabet size is "<< alphabet.size()<<", and it has "<<g<<" bd-anchors of order "<<ell<<endl;
	cout<<"The density is "<<(double) g / text_size<<endl;
	
	unsigned char * text_string = ( unsigned char * ) malloc (  ( text_size + 1 ) * sizeof ( unsigned char ) );
	ifstream is_full;
 	is_full.open (argv[1], ios::in | ios::binary);
 	
	c = 0;
	for (INT i = 0; i < text_size; i++)
	{	
		is_full.read(reinterpret_cast<char*>(&c), 1);
		text_string[i] = (unsigned char) c;
	}
	is_full.close();
	
	std::chrono::steady_clock::time_point  start_index = std::chrono::steady_clock::now();
	
	/* Constructing right and left compacted tries */
	
	vector<INT> * RSA = new vector<INT>();
	vector<INT> * RLCP = new vector<INT>();
	
	string sa_index_name = index_name + ".RSA";
 	
	ifstream in_RSA(sa_index_name, ios::binary);
	in_RSA.seekg (0, in_RSA.end);
	INT file_size_sa = in_RSA.tellg();
	
	ifstream is_RSA;
 	is_RSA.open (sa_index_name, ios::in | ios::binary);
	
	string lcp_index_name = index_name + ".RLCP";
	
	ifstream is_RLCP;
 	is_RLCP.open (lcp_index_name, ios::in | ios::binary);
 	
	ifstream in_RLCP(lcp_index_name, ios::binary);
	in_RLCP.seekg (0, in_RLCP.end);
	INT file_size_lcp = in_RLCP.tellg();
	
	vector<INT> * anchors_vector = new vector<INT>();
	for(auto &anchor: text_anchors)
		anchors_vector->push_back(anchor);
		
	text_anchors.clear();
	
	if( !(is_RSA) || !(is_RLCP )  )
	{
		ssa(text_string, anchors_vector, sa_index_name, lcp_index_name, RSA, RLCP, hash );
	 
	} 	
	else 
	{
		if( file_size_sa > 0 )
		{
			
		   	string sa = "";
		   	INT sa_int = 0;
			c = 0;
			
			for (INT i = 0; i < file_size_sa; i++)
			{	
				is_RSA.read(reinterpret_cast<char*>(&c), 1);
			
				if( (unsigned char) c == '\n' )
				{
					sa_int = stoi( sa);
					RSA->push_back( sa_int );
					sa = "";
				}
				else sa += (unsigned char) c;
				
			}
			is_RSA.close();
		
		
		}
		
	
		if( file_size_lcp > 0 )
		{	
			
		   	string lcp = "";
		   	INT lcp_int = 0;
			c = 0;
			
			for (INT i = 0; i < file_size_lcp; i++)
			{	
				is_RLCP.read(reinterpret_cast<char*>(&c), 1);
				
				if( (unsigned char) c == '\n' )
				{
					lcp_int = stoi(lcp);
					RLCP->push_back(lcp_int);
					lcp = "";
				}
				else lcp += (unsigned char) c;
				
			}
			is_RLCP.close();
		}
	
	}
	
	cout<<"Right Compacted trie constructed "<<endl;
	
	vector<INT> * LSA = new vector<INT>();
	vector<INT> * LLCP = new vector<INT>();
	
  	/* We reverse the string for the left direction and also overwrite all other DSs */
  	reverse(text_string);
  	
	sa_index_name = index_name + ".LSA";
	ifstream is_LSA;
 	is_LSA.open (sa_index_name, ios::in | ios::binary);
 	
	ifstream in_LSA(sa_index_name, ios::binary);
	in_LSA.seekg (0, in_LSA.end);
	file_size_sa = in_LSA.tellg();
	
	lcp_index_name = index_name +".LLCP" ;
	
	ifstream is_LLCP;
 	is_LLCP.open (lcp_index_name, ios::in | ios::binary);
	
	ifstream in_LLCP(lcp_index_name, ios::binary);
	in_LLCP.seekg (0, in_LLCP.end);
	file_size_lcp = in_LLCP.tellg();


	for(INT i = 0; i<anchors_vector->size(); i++)
	{
		anchors_vector->at(i) = ( text_size - 1 ) - anchors_vector->at(i);
	}
	

	if ( !(is_LSA) || !(is_LLCP) )
	{
		ssa(text_string, anchors_vector, sa_index_name, lcp_index_name, LSA, LLCP, hash );
	}
	
	else
	{
		if( file_size_sa  > 0 )
		{
		   	string sa = "";
		   	INT sa_int = 0;
			c = 0;
			
			for (INT i = 0; i < file_size_sa; i++)
			{	
				is_LSA.read(reinterpret_cast<char*>(&c), 1);
				
				if( (unsigned char) c == '\n')
				{
					sa_int = stoi( sa);
					LSA->push_back( sa_int );
					sa = "";
				}
				else sa += (unsigned char) c;
				
			}
			is_LSA.close();	
		}
	
		if( file_size_lcp  > 0 )
		{
		   	string lcp = "";
		   	INT lcp_int = 0;
			c = 0;
			
			for (INT i = 0; i < file_size_lcp; i++)
			{	
				is_LLCP.read(reinterpret_cast<char*>(&c), 1);
				
				if( (unsigned char) c == '\n' )
				{
					lcp_int = stoi( lcp );
					LLCP->push_back(lcp_int);
					lcp = "";
				}
				else lcp += (unsigned char) c;
				
			}
			is_LLCP.close();	
		}
	}
	
	cout<<"Left Compacted trie constructed "<<endl;

  	delete( anchors_vector );
	/* After constructing the tries these DSs over the whole string are not needed anymore, our data structure must be of size O(g) */

  	/* The following RMQ data structures are used for spelling pattern over the LSA and RSA */
  
	string rmq_left_suffix = index_name +".lrmq";
  		
	ifstream in_rmq_left(rmq_left_suffix, ios::binary);
  	rmq_succinct_sct<> lrmq;
  	
  	
  	if( in_rmq_left )
  	{
  	
  		load_from_file(lrmq, rmq_left_suffix); 
	}
  	else
  	{
	  	int_vector<> llcp( g , 0 ); // create a vector of length n and initialize it with 0s

		
		for ( INT i = 0; i < g; i ++ )
		{
			llcp[i] = LLCP->at(i);
			
		}

		util::assign(lrmq, rmq_succinct_sct<>(&llcp));
		
		util::clear(llcp);

		store_to_file(lrmq, rmq_left_suffix);
	}
	
	cout<<"Left RMQ DS constructed "<<endl;
	string rmq_right_suffix = index_name+ ".rrmq";
	
	ifstream in_rmq_right(rmq_right_suffix, ios::binary);
  	  
  	int_vector<> rlcp( g , 0 ); // create a vector of length n and initialize it with 0s
  	rmq_succinct_sct<> rrmq;
  	
  	if( in_rmq_right )
  	{

  		load_from_file(rrmq, rmq_right_suffix);
  	}
  	else
  	{
		int_vector<> rlcp( g , 0 ); // create a vector of length n and initialize it with 0s

		for ( INT i = 0; i < g; i ++ )
		{
			rlcp[i] = RLCP->at(i);
		}
		
		util::assign(rrmq, rmq_succinct_sct<>(&rlcp));
		
		util::clear(rlcp);
		 
	  	
	  	
	  	store_to_file(rrmq, rmq_right_suffix);
	}	
	 
	cout<<"Right RMQ DS constructed "<<endl;
  	cout<<"The whole index is constructed"<<endl;
	reverse( text_string );
	std::chrono::steady_clock::time_point  end_index = std::chrono::steady_clock::now();
	std::cout <<"Index construction took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_index- start_index + end_bd - start_bd).count() << " [ms]" << std::endl;

	std::chrono::steady_clock::time_point  begin_pt = std::chrono::steady_clock::now();
  	INT num_seqs = 0;           // the total number of patterns considered
	INT max_len_pattern = 0;
	INT ALLOC_SIZE = 180224;
	INT seq_len = 0;
	INT max_alloc_seq_len = 0;
	INT max_alloc_seqs = 0;
	unsigned char ** patterns = NULL;
	
	while ( is_patterns.read(reinterpret_cast<char*>(&c), 1) )
	{
		if( num_seqs >= max_alloc_seqs )
		{
			patterns = ( unsigned char ** ) realloc ( patterns,   ( max_alloc_seqs + ALLOC_SIZE ) * sizeof ( unsigned char* ) );
			patterns[ num_seqs ] = NULL;
			
			max_alloc_seqs += ALLOC_SIZE;
		}
		
		if( seq_len != 0 && c == '\n' )
		{
			patterns[ num_seqs ][ seq_len ] = '\0';
			
			num_seqs++;

			if( seq_len > max_len_pattern)
				max_len_pattern = seq_len;
			
			seq_len = 0;
			max_alloc_seq_len = 0;
			
			patterns[ num_seqs ] = NULL;
		}
		else 
		{
			if ( seq_len >= max_alloc_seq_len )
			{
				patterns[ num_seqs ] = ( unsigned char * ) realloc ( patterns[ num_seqs ],   ( max_alloc_seq_len + ALLOC_SIZE ) * sizeof ( unsigned char ) );
				max_alloc_seq_len += ALLOC_SIZE;
			}
			
			patterns[ num_seqs ][ seq_len ] = (unsigned char) c;	
			seq_len++;	
		}
	} 
	is_patterns.close();
	
	INT *f = new INT[ell<<1];
  	ofstream pattern_output;
	pattern_output.open(output_filename);
	
	unsigned char * left_pattern = ( unsigned char * ) malloc (  ( max_len_pattern + 1 ) * sizeof ( unsigned char ) );
	unsigned char * first_window = ( unsigned char * ) malloc (  ( max_len_pattern + 1 ) * sizeof ( unsigned char ) );
	unsigned char * right_pattern = ( unsigned char * ) malloc (  ( max_len_pattern + 1 ) * sizeof ( unsigned char ) );
			
	INT hits = 0;
	for(INT i = 0; i<num_seqs; i++)
   	{
 		INT pattern_size = strlen( (char*) patterns[i] );
   	
  		if ( pattern_size < ell )
  		{
  			pattern_output<< patterns[i] << " skipped: its length is less than ell!\n";
  			continue;
  		}
		
		memcpy( &first_window[0], &patterns[i][0], ell );
		first_window[ell] = '\0';
		
  		INT j = red_minlexrot( first_window, ell, k, power );
  		
		if ( pattern_size - j >= j ) //if the right part is bigger than the left part, then search the right part to get a smaller interval on RSA (on average)
		{ 
			
			INT right_pattern_size = pattern_size-j;
			memcpy( &right_pattern[0], &patterns[i][j], pattern_size-j );
			right_pattern[pattern_size - j] = '\0';
			
			pair<INT,INT> right_interval = pattern_matching ( right_pattern, text_string, RSA, RLCP, rrmq, g, right_pattern_size, text_size );
  												

			if(right_interval.first > right_interval.second)
			{
  				pattern_output<< patterns[i] << " was not found in the text!\n";
				continue;
			}	
		
			for(INT t = right_interval.first; t <= right_interval.second; t++ ) //this can be a large interval and only one occurrence is valid.
			{
				INT index = RSA->at(t);
				INT jj = j;		//this is the index of the anchor in the pattern
				index--; 	jj--;	//jump the index of the anchor and start looking on the left
				while ( ( jj >= 0 ) && ( index >= 0 ) && ( text_string[index] == patterns[i][jj] ) )
				{
					index--; jj--;
				}
				if ( jj < 0 ) //we have matched the pattern completely
				{
					pattern_output<< patterns[i] <<" found at position "<< index + 1 << " of the text"<<endl;
					hits++;
				}					
				else	pattern_output<< patterns[i] << " was not found in the text!\n";
			}
		}
		else //otherwise, search the left part to get a smaller interval on LSA (on average)
		{ 
			INT s = 0;
			INT left_pattern_size  = j+1;
			for(INT a = j; a>=0; a--)
			{
				left_pattern[s] = patterns[i][a];
				s++;
			}
			left_pattern[j+1] = '\0';
			
			
			pair<INT,INT> left_interval = rev_pattern_matching ( left_pattern, text_string, LSA, LLCP, lrmq, g, left_pattern_size, text_size );
  														
			if(left_interval.first > left_interval.second)	
			{
  				pattern_output<< patterns[i] << " was not found in the text!\n";
				continue;
			}
			for(INT t = left_interval.first; t <= left_interval.second; t++ ) //this can be a large interval and only one occurrence is valid.
			{
				INT index = n-1-LSA->at(t);
				INT jj = j;		//this is the index of the anchor in the pattern
				index++; 	jj++;	//jump the index of the anchor and start looking on the right
				while ( ( jj < pattern_size ) && ( index < n ) && ( text_string[index] == patterns[i][jj] ) )
				{
					index++; jj++;
				}
				if ( jj == pattern_size ) //we have matched the pattern completely
				{ 
					if ( index == n - 1 )	
						pattern_output<< patterns[i] <<" found at position "<< index - pattern_size + 1 << " of the text"<<endl;					
					else			
						
						pattern_output<< patterns[i] <<" found at position "<<  index - pattern_size << " of the text"<<endl;
					hits++;
				}
				else	pattern_output<< patterns[i] << " was not found in the text!\n";
			}
			
		}
	
   	}
 	
 	std::chrono::steady_clock::time_point  end_pt = std::chrono::steady_clock::now();
	std::cout <<"Pattern matching took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_pt - begin_pt).count() << " [ms]" << std::endl;
	std::cout <<"Occurrences: "<< hits <<endl;
  	
	for( INT i = 0; i < num_seqs; i ++ )
        	free (patterns[i]);
        free (patterns);
	free ( RSA );
  	free ( RLCP );
  	free ( LSA );
  	free ( LLCP );
  	free( left_pattern );
  	free( first_window );
  	free( right_pattern );
  	free( text_string );

	return 0;
}

