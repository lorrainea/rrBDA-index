/**
    rrBDA-index_I_ext: Randomized Reduced Bi-directional Anchors (using external memory)
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
#include "grid.h"
#include "rrbda-index_I.h"
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

//unordered_set<INT> draws;
typedef grid_point point;
typedef grid_query query;

double vm, vm0, rss, rss0;

/* Sorting comparison */
bool sort_sa(const pair<INT,INT> &a,const pair<INT,INT> &b)
{
       return a.first<b.first;
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

/* Constructs the right compacted trie given the anchors and the SA of the whole string in O(n) time */
void right_compacted_trie ( unordered_set<INT> &anchors, INT n, INT * RSA, INT * RLCP, INT g, INT ram_use, string sa_fname, string lcp_fname )
{
	stream_reader<uint40>* SA =  new stream_reader <uint40> (sa_fname, ram_use);
	stream_reader<uint40>* LCP =  new stream_reader <uint40> (lcp_fname, ram_use);
	
	INT ii = 0;
	INT minLCP = n;
	INT prevSA = 0;
	INT currSA = 0;
	INT currLCP = 0;

	for( INT i = 0; i <= n; i++ ) // in lex order
	{
		
		/* If the ith lex suffix is an anchor then add it to the compacted trie (encoded in arrays RSA and RLCP) */
		prevSA = currSA;
		currSA = SA->read();
		currLCP = LCP->read();
		
		auto it = anchors.find( currSA );
		if( it != anchors.end() )
		{
			
			RSA[ii] = currSA;		// store this suffix
			
			if ( ii == 0 )	RLCP[ii] = 0; 	// if it is the first time the RLCP = LCP = 0
			else
			{
				if ( prevSA == RSA[ii-1] )     // if the immediately prior suffix was added
					RLCP[ii] = currLCP;	// then the LCP value is the correct one for RLCP
				else
					RLCP[ii] = std::min(minLCP, currLCP);	//otherwise, we should take the minimum in the range
					
				
			}
  			//cout<<"RSA[i]: "<< RSA[ii] <<" RLCP[i]: "<<RLCP[ii]<<"\n"; getchar();
			minLCP = n; // set this to something high to get the _FIRST_ next minimum value, because a new range STARTS
			ii++;
		}
		else /* Do not add it but remember the minLCP seen so far in a range*/
		{
			if ( currLCP < minLCP )
				minLCP = currLCP;
		}
	}
	
	delete(SA);
	delete(LCP);
}

/* Constructs the left compacted trie given the anchors and the SA of the whole string in O(n) time */
void left_compacted_trie ( unordered_set<INT> &anchors, INT n, INT * LSA, INT * LLCP, INT g, INT ram_use, string sa_fname, string lcp_fname )
{

	stream_reader<uint40>* SA =  new stream_reader <uint40> (sa_fname, ram_use);
	stream_reader<uint40>* LCP =  new stream_reader <uint40> (lcp_fname, ram_use);
	
		
	INT ii = 0;
	INT minLCP = n;
	INT prevSA = 0;
	INT currSA = 0;
	INT currLCP = 0;

	for( INT i = 0; i <= n; i++ ) // in lex order
	{
		/* If the ith lex suffix is an anchor then add it to the compacted trie (encoded in arrays RSA and RLCP) */
		prevSA = currSA;
		currSA = SA->read();
		currLCP = LCP->read();
		auto it = anchors.find( ( n - 1 ) - currSA );
		//cout<<currSA<<" "<<currLCP<<endl;
		if( it != anchors.end() )
		{
			
			LSA[ii] = currSA;		// store this suffix

			if ( ii == 0 )	LLCP[ii] = 0; 	// if it is the first time the RLCP = LCP = 0
			else
			{
				if ( prevSA == LSA[ii-1] ) // if the immediately prior suffix was added
					LLCP[ii] = currLCP;	//then the LCP value is the correct one for RLCP
				else
					LLCP[ii] = std::min(minLCP, currLCP);	//otherwise, we should take the minimum in the range
				
			}

			//cout<<"LSA[i]: "<< RSA[ii] <<" LLCP[i]: "<< RLCP[ii]<<"\n"; getchar();
			minLCP = n; //set this to something high to get the FIRST next minimum value
			ii++;
		}
		else /* Do not add it but remember the minLCP seen so far in a range*/
		{
			if ( currLCP < minLCP )
				minLCP = currLCP;
		}
	}
	
	delete(SA);
	delete(LCP);
}

int main(int argc, char **argv)
{
	unordered_set<char> alphabet;

	if( argc < 7 )
 	{
        	cout<<"Wrong arguments!\n";
 		cout<<"./rbda-index_I <text_file> <ell> <pattern_file> <block_size> <ram_use> <output_filename> <index_filename>\n";
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

    	std::chrono::steady_clock::time_point  start_bd = std::chrono::steady_clock::now();
 	unordered_set<INT> text_anchors;

	string bd = index_name + ".bd";
 	
 	ifstream is_bd_anchors;
 	is_bd_anchors.open (bd, ios::in | ios::binary);
 	
	ifstream in_bd_anchors(bd, ios::binary);
 	in_bd_anchors.seekg (0, in_bd_anchors.end);
   	INT file_size = in_bd_anchors.tellg();
   	string bd_anchor = "";
   	INT bd_anchor_int = 0;
   	
   	char c = 0;
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
	
   	if( file_size > 0 )
	{
		// Read in from .bd file
	    	c = 0;
		for (INT i = 0; i < file_size; i++)
		{	
			is_bd_anchors.read(reinterpret_cast<char*>(&c), 1);
		
			if( (unsigned char) c == '\n' )
			{
				bd_anchor_int = stoi( bd_anchor );
				text_anchors.insert( bd_anchor_int );
				bd_anchor = "";
			}
			else bd_anchor += (unsigned char) c;
			
		}
		is_bd_anchors.close();
	}
	else
	{
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
		
		ofstream bd_output;
		bd_output.open(bd);
			
		for (auto &anchor : text_anchors)	
			bd_output<<anchor<<endl;
				
		bd_output.close();
		
		is_block.close();
		free( text_block );
		free( suffix_block );
		free( rank );	
		
	}
			
	INT g = text_anchors.size();
	INT n = text_size;
	
	std::chrono::steady_clock::time_point  end_bd = std::chrono::steady_clock::now();
	std::cout <<"bd construction took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_bd - start_bd).count() << "[ms]" << std::endl;
	cout<<"The text is of length "<< n << ", its alphabet size is "<< alphabet.size()<<", and it has "<<g<<" bd-anchors of order "<<ell<<endl;
	cout<<"The density is "<<(double) g / text_size<<endl;
	
	string text_string = "";
	ifstream is_full;
 	is_full.open (argv[1], ios::in | ios::binary);
 	
	c = 0;
	for (INT i = 0; i < text_size; i++)
	{	
		is_full.read(reinterpret_cast<char*>(&c), 1);
		
		if( (unsigned char) c == '\n'  )
			continue;
			
		else text_string.push_back( (unsigned char) c );
	}
	is_full.close();
	
	std::chrono::steady_clock::time_point  start_index = std::chrono::steady_clock::now();
	
	/* Constructing right and left compacted tries */
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
	
	string rsa = index_name + ".RSA";
	
	ifstream is_RSA;
 	is_RSA.open (rsa, ios::in | ios::binary);
 	
	ifstream in_RSA(rsa, ios::binary);

   	string sa_fname = index_name + "_SA.sa5";
	ifstream in_SA(sa_fname, ios::binary);
	in_RSA.seekg (0, in_RSA.end);
	INT file_size_sa = in_RSA.tellg();
	
	if( !(in_SA)  )
	{

	  	char commandesa[ sa_fname.length() + 1000 ];
	  	char * fullpathstart = dirname(realpath(argv[0], NULL));
	  	char command1[ sa_fname.length() + 1000 ];
	  	strcpy(command1, fullpathstart);
	  	strcat(command1, "/psascan/construct_sa %s -m %ldMi -o %s");
	  	sprintf(commandesa, command1, argv[1], ram_use, sa_fname.c_str());
	  	int outsa=system(commandesa);
	  	
	}
	
	if( file_size_sa > 0 )
	{
	   	string sa = "";
	   	INT sa_int = 0;
		c = 0;
		INT p = 0;
		
		for (INT i = 0; i < file_size_sa; i++)
		{	
			is_RSA.read(reinterpret_cast<char*>(&c), 1);
			
			if( (unsigned char) c == '\n' )
			{
				sa_int = stoi( sa);
				RSA[p] = sa_int;
				sa = "";
				p++;
			}
			else sa += (unsigned char) c;
			
		}
		is_RSA.close();
	
	
	}
	
	string rlcp =  index_name + ".RLCP";
	
	ifstream is_RLCP;
 	is_RLCP.open (rlcp, ios::in | ios::binary);
 	
	ifstream in_RLCP(rlcp, ios::binary);
	in_RLCP.seekg (0, in_RLCP.end);
	file_size = in_RLCP.tellg();
	
   	
   	string lcp_fname = index_name + "_LCP.lcp5";
	ifstream in_LCP(lcp_fname, ios::binary);
	
	if ( !(in_LCP) || file_size <= 0 || file_size_sa <= 0 )
	{
		if( !(in_LCP ) )
		{
			char commande[ sa_fname.length() + lcp_fname.length() + 1000];
			char * fullpathstart = dirname(realpath(argv[0], NULL));
			char command2[ sa_fname.length() + lcp_fname.length()  + 1000];
			strcpy(command2, fullpathstart);
			strcat(command2, "/sparsePhi/src/construct_lcp_parallel -m %ldG -o %s -s %s %s");
			sprintf(commande, command2, ram_use, lcp_fname.c_str(), sa_fname.c_str(), argv[1]);
			int out=system(commande);
		}
		
	    	
	  	right_compacted_trie ( text_anchors, n, RSA, RLCP, g, ram_use, sa_fname, lcp_fname );
	  	
	  	ofstream rsa_output;
		rsa_output.open(rsa);
		
		for(INT i = 0; i<g; i++)	
			rsa_output<<RSA[i]<<endl;
			
		rsa_output.close();
		
		ofstream rlcp_output;
		rlcp_output.open(rlcp);
		
		for(INT i = 0; i<g; i++)	
			rlcp_output<<RLCP[i]<<endl;
			
		rlcp_output.close();
		
		cout<<"Right Compacted trie constructed "<<endl;	
	}
	
	if( file_size > 0 )
	{
	   	string lcp = "";
	   	INT lcp_int = 0;
		c = 0;
		INT p = 0;
		
		for (INT i = 0; i < file_size; i++)
		{	
			is_RLCP.read(reinterpret_cast<char*>(&c), 1);
			
			if( (unsigned char) c == '\n' )
			{
				lcp_int = stoi(lcp);
				RLCP[p] = lcp_int;
				lcp = "";
				p++;
			}
			else lcp += (unsigned char) c;
			
		}
		is_RLCP.close();
	
	
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
	
	char * output_reverse;
	const char * reversed_text = "_reverse";
	
	output_reverse = (char *) malloc(strlen(argv[1])+9);
	strcpy( output_reverse, argv[1]);
	strcat(output_reverse, reversed_text);
		
  	/* We reverse the string for the left direction and also overwrite all other DSs */
  	std::ofstream output_r;
  	output_r.open (output_reverse);
  	reverse(text_string.begin(), text_string.end());
  	output_r << text_string;
    	output_r.close();
 
	string lsa = index_name + ".LSA";
  	
  	ifstream is_LSA;
 	is_LSA.open (lsa, ios::in | ios::binary);
 	
  	ifstream in_LSA(lsa, ios::binary);
  	
 
  	in_LSA.seekg (0, in_LSA.end);
  	file_size_sa = in_LSA.tellg();
  	
   	string sa_fname_reverse = index_name +"reverse_SA.sa5";
	ifstream in_SA_reverse(sa_fname_reverse, ios::binary);
	
	if ( !(in_SA_reverse)  )
	{
	  	char commandesa_reverse[ sa_fname_reverse.length() + 1000 ];
	  	char * fullpathstart_reverse = dirname(realpath(argv[0], NULL));
	  	char command1_reverse[ sa_fname_reverse.length() + 1000 ];
	  	strcpy(command1_reverse, fullpathstart_reverse);
	  	strcat(command1_reverse, "/psascan/construct_sa %s -m %ldMi -o %s");
	  	sprintf(commandesa_reverse, command1_reverse, output_reverse, ram_use, sa_fname_reverse.c_str());
	  	int outsa_reverse=system(commandesa_reverse);
	}
	
	if( file_size_sa  > 0 )
	{
	   	string sa = "";
	   	INT sa_int = 0;
		c = 0;
		INT p = 0;
		
		for (INT i = 0; i < file_size_sa; i++)
		{	
			is_LSA.read(reinterpret_cast<char*>(&c), 1);
			
			if( (unsigned char) c == '\n' )
			{
				sa_int = stoi( sa);
				LSA[p] = sa_int;
				sa = "";
				p++;
			}
			else sa += (unsigned char) c;
			
		}
		is_LSA.close();	
	}
	

	string llcp = index_name + ".LLCP";
	ifstream is_LLCP;
 	is_LLCP.open (llcp, ios::in | ios::binary);
	
	ifstream in_LLCP(llcp, ios::binary);
	in_LLCP.seekg (0, in_LLCP.end);
	file_size = in_LLCP.tellg();
   	
 	string lcp_fname_reverse = index_name + "_reverse_LCP.lcp5";
        ifstream in_LCP_reverse(lcp_fname_reverse, ios::binary);
	
	
        if( !(in_LCP_reverse) || file_size <= 0 || file_size_sa <= 0 )
	{
	
		if( !(in_LCP_reverse ) )
		{
			char commande_reverse[ sa_fname_reverse.length() + lcp_fname_reverse.length() + 1000];
			char * fullpathstart_reverse = dirname(realpath(argv[0], NULL));
			char command2_reverse[ sa_fname_reverse.length() + lcp_fname_reverse.length() + 1000];
			strcpy(command2_reverse, fullpathstart_reverse);
			strcat(command2_reverse, "/sparsePhi/src/construct_lcp_parallel -m %ldG -o %s -s %s %s");
			sprintf(commande_reverse, command2_reverse, ram_use, lcp_fname_reverse.c_str(), sa_fname_reverse.c_str(), output_reverse);
			int out_reverse=system(commande_reverse);
		}
		
		left_compacted_trie ( text_anchors, n, LSA, LLCP, g, ram_use, sa_fname_reverse, lcp_fname_reverse );
  		
  		ofstream lsa_output;
		lsa_output.open(lsa);
		
		for(INT i = 0; i<g; i++)	
			lsa_output<<LSA[i]<<endl;
			
		lsa_output.close();
		
		ofstream llcp_output;
		llcp_output.open(llcp);
		
		for(INT i = 0; i<g; i++)	
			llcp_output<<LLCP[i]<<endl;
			
		llcp_output.close();
	
		cout<<"Left Compacted trie constructed"<<endl;
	}
	
	
	if( file_size > 0 )
	{
	   	string lcp = "";
	   	INT lcp_int = 0;
		c = 0;
		INT p = 0;
		
		for (INT i = 0; i < file_size; i++)
		{	
			is_LLCP.read(reinterpret_cast<char*>(&c), 1);
			
			if( (unsigned char) c == '\n' )
			{
				lcp_int = stoi(lcp);
				LLCP[p] = lcp_int;
				lcp = "";
				p++;

			}
			else lcp += (unsigned char) c;
			
		}
		is_LLCP.close();
	}
	/* After constructing the tries these DSs over the whole string are not needed anymore, our data structure must be of size O(g) */
  	text_anchors.clear();
  	
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
	  	int_vector<> llcp_rmq( g , 0 ); // create a vector of length n and initialize it with 0s

		
		for ( INT i = 0; i < g; i ++ )
		{
			llcp_rmq[i] = LLCP[i];
			
		}

		util::assign(lrmq, rmq_succinct_sct<>(&llcp_rmq));
		
		util::clear(llcp_rmq);

		store_to_file(lrmq, rmq_left_suffix);
	}
	
	cout<<"Left RMQ DS constructed "<<endl;
	string rmq_right_suffix = index_name+ ".rrmq";
	
	ifstream in_rmq_right(rmq_right_suffix, ios::binary);
  	  
  	int_vector<> rlcp_rmq( g , 0 ); // create a vector of length n and initialize it with 0s
  	rmq_succinct_sct<> rrmq;
  	
  	if( in_rmq_right )
  	{

  		load_from_file(rrmq, rmq_right_suffix);
  	}
  	else
  	{
		int_vector<> rlcp_rmq( g , 0 ); // create a vector of length n and initialize it with 0s

		for ( INT i = 0; i < g; i ++ )
		{
			rlcp_rmq[i] = RLCP[i];
		}
		
		util::assign(rrmq, rmq_succinct_sct<>(&rlcp_rmq));
		
		util::clear(rlcp_rmq);
		 
	  	
	  	
	  	store_to_file(rrmq, rmq_right_suffix);
	}	
	 
	cout<<"Right RMQ DS constructed "<<endl;
  	//Forming coordinate points using LSA and RSA

	string grid_suffix =  index_name+".grid";
  		
  	ifstream in_grid(grid_suffix, ios::binary);
  	in_grid.seekg (0, in_grid.end);
  	file_size = in_grid.tellg();
  	
  	ifstream is_grid;
 	is_grid.open (grid_suffix, ios::in | ios::binary);
   	
   	std::vector<point> points;
   	grid construct;
   	
 
   	if( file_size > 0 )
  	{  	
  		load_from_file(construct, grid_suffix);
  	}
  	else
  	{
	  	vector<pair<INT,INT>> l;
		vector<pair<INT,INT>> r;
		for ( INT i = 0; i < g; i++ )
	  	{
	  		l.push_back( make_pair( n-1-LSA[i], i) );
	  		r.push_back( make_pair( RSA[i], i ) );
	 	}
	 	
	 	sort(l.begin(),l.end(),sort_sa);
		sort(r.begin(),r.end(),sort_sa);
		
	  	for ( INT i = 0; i < g; i++ )
	  	{
	 
			point to_insert;
			to_insert.row = l.at(i).second+1;
			to_insert.col = r.at(i).second+1;
			
			to_insert.level = 0;
			to_insert.label = l.at(i).first;
			points.push_back(to_insert); 
	  	}
	  		
		construct.build( points, 0 );
		
		store_to_file( construct, grid_suffix );	
	}
	
	
	cout<<"The grid is constructed"<<endl;  
  	cout<<"The whole index is constructed"<<endl;
  	
	std::chrono::steady_clock::time_point  end_index = std::chrono::steady_clock::now();
	std::cout <<"Index took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_index- start_index).count() << "[ms]" << std::endl;

	vector<vector<unsigned char> > all_patterns;
    	vector<unsigned char> pattern;
    	c = 0;
    	while (is_patterns.read(reinterpret_cast<char*>(&c), 1))
    	{
        	if(c == '\n')
        	{
  			if(pattern.empty())	break;
  			all_patterns.push_back(pattern);
  			pattern.clear();
        	}
        	else	pattern.push_back((unsigned char)c);
    	}
    	is_patterns.close();
    	pattern.clear();

	vector<string> new_all_pat;
	for(auto &it_pat : all_patterns)	new_all_pat.push_back(string(it_pat.begin(), it_pat.end()));
	all_patterns.clear();
  	
	std::chrono::steady_clock::time_point  begin_pt = std::chrono::steady_clock::now();
    	reverse(text_string.begin(), text_string.end()); 				//I re-reverse to take the original string
  
	ofstream pattern_output;
	pattern_output.open(output_filename);
	for(auto &pattern : new_all_pat)
   	{
  		if ( pattern.size() < ell )
  		{
  			pattern_output<<"Pattern skipped: its length is less than ell!\n";
  			continue;
  		}
  		
		string first_window = pattern.substr(0, ell);
		INT j = red_minlexrot( first_window, ell, k, power );
		string left_pattern = pattern.substr(0, j+1);
	  	reverse(left_pattern.begin(), left_pattern.end());
	  		
	  	
		pair<INT,INT> left_interval = rev_pattern_matching ( left_pattern, text_string, LSA, LLCP, lrmq, g );  
	  		
		string right_pattern = pattern.substr(j, pattern.size()-j);
		pair<INT,INT> right_interval = pattern_matching ( right_pattern, text_string, RSA, RLCP, rrmq, g );
	  	if ( left_interval.first <= left_interval.second  && right_interval.first <= right_interval.second )
	  	{
	  		//Finding rectangle containing bd-anchors in grid
	  		grid_query rectangle;
	  		rectangle.row1 = left_interval.first+1;
	  		rectangle.row2 = left_interval.second+1;
	  		rectangle.col1 = right_interval.first+1;
	  		rectangle.col2 = right_interval.second+1;
	  		
				
			vector<long unsigned int> result;
			construct.search_2d(rectangle, result);
			
			for(INT i = 0; i<result.size(); i++)
			{
				pattern_output<<pattern<<" found at position "<<RSA[result.at(i)-1]-j<<" of the text"<<endl;	
			}
		
  		}
	 	else	pattern_output<<"No occurrences found!\n";
			
		
  		
 	}
	pattern_output.close();
 	
	std::chrono::steady_clock::time_point  end_pt = std::chrono::steady_clock::now();
	std::cout <<"Pattern matching took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_pt - begin_pt).count() << "[ms]" << std::endl;
  	delete ( RSA );
  	delete ( RLCP );
  	delete ( LSA );
  	delete ( LLCP );

	return 0;
}