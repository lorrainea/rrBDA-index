#include <unordered_set>
#include "utils.h"
#include "stream.h"
#include "uint40.h"
#include <math.h>
#include "bda-index_I.h"
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
INT bd_anchors(  unsigned char * seq, INT pos, INT ell, INT k, unordered_set<INT> &anchors, INT * FP )
{

	cout<<seq<<" k:"<<k<<endl;
	INT w = ell;
	INT n = strlen ( (char*) seq );
	INT fp = 0;

        for(INT j = 0; j<k; j++)
                fp =  karp_rabin_hashing::concat( fp, seq[j] , 1 );

        FP[0] = fp; 
        
        deque<pair<INT,utils::FP_>> min_fp = {};
	vector<utils::FP_> minimizers;
	
        // find all fingerprints for all k substrings
        for( INT j = 1; j<n; j++)
        {
                fp = karp_rabin_hashing::concat( fp, seq[j+k-1] , 1 );
                fp = karp_rabin_hashing::subtract( fp, seq[j-1] , k );

                FP[j] = fp; 
        }

	for( INT j = 0; j<=n-k; j++)
        {
               cout<<" fp "<<FP[j]<<" "<<j<<" "<<seq[j]<<endl;
        } 
        
     
   	for ( INT j = 0; j < w - k - 1; j++) 
   	{
 		while ( !min_fp.empty() && FP[j] < min_fp.back().first )
 			min_fp.pop_back();
 
       		utils::FP_ potential_bd;
		potential_bd.start_pos = j;
		potential_bd.fp_pos = j;
				
		min_fp.push_back(std::make_pair(FP[j], potential_bd));
		
    	}
    	
	/* Compute reduced bd-anchors for every window of size ell */
	
	INT i = w - k - 1;
	for( INT j = 0; j<=n-w; j++ )
	{
		
		while (!min_fp.empty() && min_fp.back().first > FP[i])
			min_fp.pop_back();
					
		utils::FP_ potential_bd;
		potential_bd.start_pos = i;
		potential_bd.fp_pos = i;
				
		min_fp.push_back(std::make_pair(FP[i], potential_bd));
		
	
		while( min_fp.front().second.start_pos <= i - w + k)
		{
			min_fp.pop_front();
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
		
		i++;
	
		/* Filter draws if there are more than one minimum fp, otherwise only one potential bd-anchor in window */			
		if( minimizers.size() > 1 )
		{ 	
			INT minimum = 0;

			for(INT i = 1; i<minimizers.size(); i++)
			{
		
				INT dist_to_end = w;
				
				INT fp_pos = minimizers.at(i).fp_pos;
				INT min_fp_pos = minimizers.at(minimum).fp_pos;
				
				if( ( (j+ w ) - fp_pos ) < dist_to_end )
					dist_to_end = ( (j+ w ) - fp_pos );
				
				if( ( (j+ w ) -  min_fp_pos ) < dist_to_end )
					dist_to_end = ( (j+ w ) -  min_fp_pos );
				
				INT min_inv = min( seq[pos+min_fp_pos], seq[pos+fp_pos]) ;
				
					
				INT max_inv = max( seq[pos+min_fp_pos], seq[pos+fp_pos]) ;
				
				INT lcp1 = 0; 
				
				while ( seq[min_fp_pos+lcp1] == seq[fp_pos+lcp1] )
					lcp1++;
			
				
				if( lcp1 < dist_to_end )
				{
					
					if( seq[ min_fp_pos+lcp1 ] > seq[fp_pos+lcp1 ] )
					{
						minimum = i;
					}
				}
				else
				{
					
					
					min_fp_pos =  min_fp_pos + min(lcp1,dist_to_end) ;
					fp_pos = j;
				
					dist_to_end = w;
					if( ( (j+ w ) - fp_pos ) < dist_to_end )
						dist_to_end = ( (j+ w ) - fp_pos );
				
					if( ( (j+ w ) -  min_fp_pos ) < dist_to_end )
						dist_to_end = ( (j+ w ) -  min_fp_pos );
				
					INT min_inv = min( seq[min_fp_pos], seq[fp_pos]) ;
					
						
					INT lcp2 = 0; 
				
					while ( seq[min_fp_pos+lcp2] == seq[fp_pos+lcp2] )
						lcp2++;
				
					if( lcp2 < dist_to_end )
					{
						
						if( seq[ min_fp_pos+lcp2 ] > seq[ fp_pos+lcp2 ] )
						{
							minimum = i;
						}
					}
					else
					{
						
						
						
					 	min_fp_pos = min_fp_pos + min(lcp2,dist_to_end);
						fp_pos = fp_pos + min(lcp2,dist_to_end);
						dist_to_end = w;
						if( ( (minimizers.at(i).start_pos) - fp_pos ) < dist_to_end )
							dist_to_end = ( (minimizers.at(i).start_pos) - fp_pos );
						
						if( ( (minimizers.at(i).start_pos) -  min_fp_pos ) < dist_to_end )
							dist_to_end = ( (minimizers.at(i).start_pos) -  min_fp_pos );
						
						INT min_inv = min( seq[pos+min_fp_pos], seq[pos+fp_pos]) ;
						INT max_inv = max( seq[pos+min_fp_pos], seq[pos+fp_pos]) ;
						
						INT lcp3 = 0; 
						
						while ( seq[min_fp_pos+lcp3] == seq[fp_pos+lcp3] )
							lcp3++;
						
						if( lcp3 < dist_to_end )
						{
							if( seq[min_fp_pos+lcp3 ] > seq[fp_pos+lcp3 ] )
							{
								minimum = i;
							}
						}
					}
					
					
				
				}
			}
			anchors.insert( minimizers.at(minimum).start_pos+pos );
			cout<<minimizers.at(minimum).start_pos<<" "<<FP[ minimizers.at(minimum).start_pos ]<<endl;
		

		}
		else 
		{
			anchors.insert( minimizers.at(0).start_pos+pos );
			cout<<minimizers.at(0).start_pos<<" "<<FP[ minimizers.at(0).start_pos ]<<endl;
		}
			
		
		minimizers.clear();
	}
						
	return 0;
}
