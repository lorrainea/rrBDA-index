#include <unordered_set>
#include "utils.h"
#include "stream.h"
#include "uint40.h"
#include <math.h>
#include "bda-index_I.h"
#include "krfp.h"
#include "rk_lce.h"
#include "includes.h"
#include <deque>
#include <cstring>
#include "cyclichash.h"

using namespace std;

/* Kasai et al algorithm for O(n)-time LCP construction */
INT LCParray ( unsigned char * text, INT n, INT * SA, INT * ISA, INT * LCP )
{
        INT i=0, j=0;

        LCP[0] = 0;
        for ( i = 0; i < n; i++ ) // compute LCP[ISA[i]]
                if ( ISA[i] != 0 )
                {
                        if ( i == 0) j = 0;
                        else j = (LCP[ISA[i-1]] >= 2) ? LCP[ISA[i-1]]-1 : 0;
                        while ( text[i+j] == text[SA[ISA[i]-1]+j] )
                                j++;
                        LCP[ISA[i]] = j;
                }
        return ( 1 );
}

/* Computes the bd-anchors of a string of length n in O(n) time */
INT bd_anchors(  unsigned char * seq, string filename, INT ell, INT k, unordered_set<INT> &anchors, INT option  )
{

	auto lce = rklce::rk_lce("");
	
	if ( option == 1)
	{
		cout << "Building LCP data structure and allocating memory for SA ... " << endl;
	
		lce = rklce::rk_lce(filename);
	
		cout << " Size of LCE structure (Bytes) " << lce.bit_size()/8 << endl;

	}
	
	karp_rabin_hashing::init();
	
	//CyclicHash<> hf(k, 4);
  
	INT w = ell;
	INT n = strlen ( (char*) seq );
	

	INT fp = 0;
        INT smallest_fp = fp;

        INT * FP = ( INT * ) malloc( ( w  ) *  sizeof( INT ) );

        for(INT j = 0; j<k; j++)
        	//hf.eat(seq[j]);
                fp =  karp_rabin_hashing::concat( fp, seq[j] , 1 );

        FP[0] = fp; //hf.hashvalue;
        INT pos = 1;
        
        deque<pair<INT,utils::FP_>> min_fp = {};
	vector<utils::FP_> minimizers;
	
        // find all fingerprints for all k substrings in first window
        for(INT j = 1; j<=w-k; j++)
        {
        	//hf.update(seq[j+k-1], seq[j-1]);
                fp = karp_rabin_hashing::concat( fp, seq[j+k-1] , 1 );
                fp = karp_rabin_hashing::subtract( fp, seq[j-1] , k );

                FP[pos] = fp; //hf.hashvalue;
                pos++;
        }

        // minimum fp in first window
        for (INT j = 0; j <= w - k ; j++)
        {
                while ( !min_fp.empty() && FP[j] < min_fp.back().first )
                        min_fp.pop_back();

		utils::FP_ potential_bd;
		potential_bd.start_pos = j;
		potential_bd.fp_pos = j;
				
				
                min_fp.push_back(std::make_pair(FP[j], potential_bd));
        }

	/* Compute reduced bd-anchors for every window of size ell */
	INT i = w - k + 1;
	//hf.update(seq[i], seq[i-k]);
	//fp = hf.hashvalue;
	fp = karp_rabin_hashing::concat( fp, seq[i] , 1 );
        fp = karp_rabin_hashing::subtract( fp, seq[i-k] , k );
                
	for( INT j = 1; j<=n-w; j++ )
	{
		while (!min_fp.empty() && min_fp.back().first > fp)
			min_fp.pop_back();

		utils::FP_ potential_bd;
		potential_bd.start_pos = i;
		potential_bd.fp_pos = i;
				
		min_fp.push_back(std::make_pair(fp, potential_bd));
		
		while( min_fp.front().second.start_pos < i - w + k )
			min_fp.pop_front();
			
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
		//hf.update(seq[i], seq[i-k]);
		//fp = hf.hashvalue;
		fp = karp_rabin_hashing::concat( fp, seq[i] , 1 );
                fp = karp_rabin_hashing::subtract( fp, seq[i-k] , k );
                
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
				
				INT lcp1 = 0; 
				
				if( option == 1 )
					lcp1 = lce.LCE( min_fp_pos, fp_pos );	
				else
				{
				
					while ( seq[min_fp_pos+lcp1] == seq[fp_pos+lcp1] && lcp1 <= dist_to_end  )
						lcp1++;
				}

				//cout<<min_fp_pos<<" "<<fp_pos<<" "<<lcp1<<endl;
				if( lcp1 < dist_to_end )
				{
					
					if( seq[ min_fp_pos+lcp1 ] > seq[ fp_pos+lcp1 ] )
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
					
					
					INT lcp2 = 0; 
				
					if( option == 1 )
						lcp2 = lce.LCE( min_fp_pos, fp_pos );	
					else
					{
					
						while ( seq[min_fp_pos+lcp2] == seq[fp_pos+lcp2] && lcp2 <= dist_to_end )
							lcp2++;
					}
					//	cout<<min_fp_pos<<" "<<fp_pos<<" "<<lcp2<<endl;
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
						
						INT lcp3 = 0; 
				
						if( option == 1 )
							lcp3 = lce.LCE( min_fp_pos, fp_pos );	
						else
						{
						
							while ( seq[min_fp_pos+lcp3] == seq[fp_pos+lcp3] && lcp3 <= dist_to_end )
								lcp3++;
						}
						//	cout<<min_fp_pos<<" "<<fp_pos<<" "<<lcp3<<endl;
						if( lcp3 < dist_to_end )
						{
							if( seq[ min_fp_pos+lcp3 ] > seq[ fp_pos+lcp3 ] )
							{
								minimum = i;
							}
						}
					}
					
					
				
				}
			}
			anchors.insert( minimizers.at(minimum).start_pos );
		

		}
		else 
		{
			anchors.insert( minimizers.at(0).start_pos );
		}
			
		
		minimizers.clear();
	}
	
	//cout<<"LCE alphabet size "<<lce.alphabet_size()<<endl;
			
	free(FP);			
	return 0;
}

