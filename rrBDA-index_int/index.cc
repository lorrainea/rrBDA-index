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
#include "rrbda-index_int.h"
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
 

INT compute_index( uint64_t hash, string index_name, INT text_size, INT g, unsigned char * text_string, unordered_set<INT> &text_anchors, vector<INT> * RSA, vector<INT> * RLCP, vector<INT> * LSA, vector<INT> * LLCP, rmq_succinct_sct<> &lrmq, rmq_succinct_sct<> &rrmq )
{
	vector<INT> * anchors_vector = new vector<INT>();	

	unsigned char c = 0;
	
	/* Constructing right and left compacted tries */
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
	
	
	for(auto &anchor: text_anchors)
		anchors_vector->push_back(anchor);
		
	text_anchors.clear();
	
	if( !(is_RSA) || !(is_RLCP )  )
	{
		ssa(text_string, text_size, anchors_vector, sa_index_name, lcp_index_name, RSA, RLCP, hash );
	 
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
					sa_int = stol( sa);
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
					lcp_int = stol(lcp);
					RLCP->push_back(lcp_int);
					lcp = "";
				}
				else lcp += (unsigned char) c;
				
			}
			is_RLCP.close();
		}
	
	}
	
	cout<<"Right Compacted trie constructed "<<endl;
	
  	/* We reverse the string for the left direction and also overwrite all other DSs */
  	reverse(text_string, text_size);
  	
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
		ssa(text_string, text_size, anchors_vector, sa_index_name, lcp_index_name, LSA, LLCP, hash );
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
					sa_int = stol( sa);
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
					lcp_int = stol( lcp );
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
  	
  	reverse( text_string, text_size );
  	
  	

	return 0;
	
}	
