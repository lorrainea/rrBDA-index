/**
    rrBDA-index_int: Randomized Reduced Bi-directional Anchors (using external memory)
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

#include <unordered_set>
#include "utils.h"
#include "stream.h"
#include "uint40.h"
#include <math.h>
#include "rrbda-index_ext.h"
#include "krfp.h"

using namespace std;
using namespace sdsl;

#ifdef _USE_64
#include <divsufsort64.h>                                         // include header for suffix sort
#endif

#ifdef _USE_32
#include <divsufsort.h>                                       	  // include header for suffix sort
#endif

void reverse( unsigned char * &s, INT text_size)
{
	INT length = text_size ;
	
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


INT compute_index( uint64_t hash, string index_name, INT text_size, INT g, unsigned char * text_string, unordered_set<INT> &text_anchors, INT * RSA, INT * RLCP, INT * LSA, INT * LLCP, rmq_succinct_sct<> &lrmq, rmq_succinct_sct<> &rrmq, char * arg0, char * arg1, INT ram_use )
{
	/* Constructing right and left compacted tries */

	string rsa = index_name + ".RSA";
	unsigned char c = 0;
	
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
	  	char * fullpathstart = dirname(realpath(arg0, NULL));
	  	char command1[ sa_fname.length() + 1000 ];
	  	strcpy(command1, fullpathstart);
	  	strcat(command1, "/psascan/construct_sa %s -m %ldMi -o %s");
	  	sprintf(commandesa, command1, arg1, ram_use, sa_fname.c_str());
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
				sa_int = stol( sa);
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
	INT file_size = in_RLCP.tellg();
	
   	
   	string lcp_fname = index_name + "_LCP.lcp5";
	ifstream in_LCP(lcp_fname, ios::binary);
	
	if ( !(in_LCP) || file_size <= 0 || file_size_sa <= 0 )
	{
		if( !(in_LCP ) )
		{
			char commande[ sa_fname.length() + lcp_fname.length() + 1000];
			char * fullpathstart = dirname(realpath(arg0, NULL));
			char command2[ sa_fname.length() + lcp_fname.length()  + 1000];
			strcpy(command2, fullpathstart);
			strcat(command2, "/sparsePhi/src/construct_lcp_parallel -m %ldG -o %s -s %s %s");
			sprintf(commande, command2, ram_use, lcp_fname.c_str(), sa_fname.c_str(), arg1);
			int out=system(commande);
		}
		
	    	
	  	right_compacted_trie ( text_anchors, text_size, RSA, RLCP, g, ram_use, sa_fname, lcp_fname );
	  	
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
				lcp_int = stol(lcp);
				RLCP[p] = lcp_int;
				lcp = "";
				p++;
			}
			else lcp += (unsigned char) c;
			
		}
		is_RLCP.close();
	
	
	}
	
	char * output_reverse;
	const char * reversed_text = "_reverse";
	
	output_reverse = (char *) malloc(strlen(arg1)+9);
	strcpy( output_reverse, arg1);
	strcat(output_reverse, reversed_text);
		
  	/* We reverse the string for the left direction and also overwrite all other DSs */
  	std::ofstream output_r;
  	output_r.open (output_reverse);
  	
  	reverse(text_string, text_size);
  	output_r << text_string;
    	output_r.close();
 
	string lsa = index_name + ".LSA";
  	
  	ifstream is_LSA;
 	is_LSA.open (lsa, ios::in | ios::binary);
 	
  	ifstream in_LSA(lsa, ios::binary);
  	
 
  	in_LSA.seekg (0, in_LSA.end);
  	file_size_sa = in_LSA.tellg();
  	
   	string sa_fname_reverse = index_name +"_reverse_SA.sa5";
	ifstream in_SA_reverse(sa_fname_reverse, ios::binary);
	
	if ( !(in_SA_reverse)  )
	{
	  	char commandesa_reverse[ sa_fname_reverse.length() + 1000 ];
	  	char * fullpathstart_reverse = dirname(realpath(arg0, NULL));
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
				sa_int = stol( sa);
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
			char * fullpathstart_reverse = dirname(realpath(arg0, NULL));
			char command2_reverse[ sa_fname_reverse.length() + lcp_fname_reverse.length() + 1000];
			strcpy(command2_reverse, fullpathstart_reverse);
			strcat(command2_reverse, "/sparsePhi/src/construct_lcp_parallel -m %ldG -o %s -s %s %s");
			sprintf(commande_reverse, command2_reverse, ram_use, lcp_fname_reverse.c_str(), sa_fname_reverse.c_str(), output_reverse);
			int out_reverse=system(commande_reverse);
		}
		
		left_compacted_trie ( text_anchors, text_size, LSA, LLCP, g, ram_use, sa_fname_reverse, lcp_fname_reverse );
  		
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
	
	free( output_reverse );
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
				lcp_int = stol(lcp);
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
  	cout<<"The whole index is constructed"<<endl;
  	reverse(text_string, text_size); 
	
	return 0;
}	
