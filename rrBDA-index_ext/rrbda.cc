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

/* Computes the bd-anchors of a string of length n in O(n) time */
INT bd_anchors(  unsigned char * seq, INT pos, INT ell, INT k, unordered_set<INT> &anchors, INT * FP, INT power)
{

	INT w = ell;
	INT n = strlen ( (char*) seq );
	INT fp = 0;

        for(INT j = 0; j<k; j++)
                fp =  karp_rabin_hashing::concat( fp, seq[j] , 1 );

        FP[0] = fp; 
        
        deque<pair<INT,utils::FP_>> min_fp = {};
	vector<utils::FP_> minimizers;
	
        // find all fingerprints for all k substrings
        for( INT j = 1; j<=n-k; j++)
        {
                fp = karp_rabin_hashing::concat( fp, seq[j+k-1] , 1 );
                fp = karp_rabin_hashing::subtract_fast( fp, seq[j-1] , power );

                FP[j] = fp; 
               
        }


	/* Compute reduced bd-anchors for every window of size ell */
	for( INT j = 0; j<=n-w; j++ )
	{
	
		if( j == 0 )
		{
			for ( INT l = 0; l <= w-k; l++) 
		   	{
		 		while ( !min_fp.empty() && FP[l] < min_fp.back().first )
		 			min_fp.pop_back();
		 
			       	utils::FP_ potential_bd;
				potential_bd.start_pos = l;
				potential_bd.fp_pos = l;
							
				min_fp.push_back(std::make_pair(FP[l], potential_bd));
		    	}
		
		}
		else
		{
			while( min_fp.front().second.start_pos < j )
				min_fp.pop_front();
				
			while (!min_fp.empty() && min_fp.back().first > FP[j+w-k])
				min_fp.pop_back();
					
			utils::FP_ potential_bd;
			potential_bd.start_pos = j+w-k;
			potential_bd.fp_pos = j+w-k;
				
			min_fp.push_back(std::make_pair(FP[j+w-k], potential_bd));
			
		}
			
		INT min_ = min_fp.at(0).first;
		for(INT i = 0; i<min_fp.size(); i++)
		{
			if( min_fp.at(i).first == min_ )
			{
				minimizers.push_back( min_fp.at(i).second );
			}
			else if( min_fp.at(i).first >  min_ )
				break;
		}
		
		char a = ' ';
		char b = ' ';
		
		INT a_pos = 0;
		INT b_pos = 0;
		bool cont = false;
		
		/* Filter draws if there are more than one minimum fp, otherwise only one potential bd-anchor in window */			
		if( minimizers.size() > 1 )
		{ 	
			INT smallest_fp_pos = minimizers.at(0).start_pos;
		
       		
			for(INT i = 1; i<minimizers.size(); i++ )
			{
				a_pos = smallest_fp_pos+k;
				b_pos = minimizers.at(i).start_pos+k;
				
				cont = false;
			
				for(INT c = k; c<w; c++)
				{
					a = seq[a_pos];
					b = seq[b_pos];

					if( b_pos >= j + w )
					{
						b_pos = j;
						cont = true;
						break;
					}
					
					if( b < a )
					{
						smallest_fp_pos =  minimizers.at(i).start_pos;
						break;
					}
					else if( b > a )
						break;
					
					a_pos++;
					b_pos++;
				}
				
				if( cont == true )
				{
					cont  = false;
					for(INT c = 0; c<w; c++)
					{
		      				a = seq[a_pos];
						b = seq[b_pos];
							
						if( a_pos >= j + w )
						{
							a_pos = j;
							cont = true;
							break;
						}
					
						if( b < a )
						{
							smallest_fp_pos =  minimizers.at(i).start_pos;
							break;
						}
						else if( b > a )
							break;
						
						a_pos++;
						b_pos++;
					}
				}
				
				
				if( cont == true )
				{
					for(INT c = 0; c<w; c++)
					{
		      				a = seq[a_pos];
						b = seq[b_pos];
							
						if( b_pos >= j + w || b < a )
						{
							smallest_fp_pos =  minimizers.at(i).start_pos;
							break;
						}
						else if( b > a )
							break;
						
						a_pos++;
						b_pos++;
					}
					
				}	
			}	
			
			anchors.insert( smallest_fp_pos+pos );
		}	
		else 
		{
			anchors.insert( minimizers.at(0).start_pos+pos );
		}
			
		
		minimizers.clear();
	}
						
	return 0;
}
