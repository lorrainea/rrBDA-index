/**
    rrBDA-index_int: Randomized Reduced Bi-directional Anchors (using intenrnal memory)
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

INT red_minlexrot( unsigned char * X, INT n, uint64_t k, uint64_t power )
{  
  	uint64_t fp = 0;
    	vector<INT> * draws = new vector<INT>();
	
	for(uint64_t j = 0; j<k; j++)
        	fp =  karp_rabin_hashing::concat( fp, X[j] , 1 );

	draws->push_back(0);
	uint64_t smallest_fp = fp;
  	INT smallest_fp_pos = 0;
  	
  	for(INT j = 1; j<=n-k; j++)
        {
                fp = karp_rabin_hashing::concat( fp,  X[j+k-1] , 1 );
                fp = karp_rabin_hashing::subtract_fast( fp, X[j-1] , power );

                if( fp < smallest_fp )
                {
                	draws->clear();
                	
                	smallest_fp = fp;
                	smallest_fp_pos = j;
                	
                	draws->push_back( j );
                
                }
                else if( fp == smallest_fp )
                {	
                	draws->push_back( j );
                }	
                
       
        }

 	if( draws->size() > 1 )
        {
        	smallest_fp = draws->at(0);
        	
       		for(INT i = 1; i<draws->size(); i++ )
		{
			INT a_pos = smallest_fp_pos+k;
			INT b_pos = draws->at(i)+k;
				
			bool cont = false;
			for(INT j = k; j<n; j++)
			{
				unsigned char a = X[a_pos];
				unsigned char b = X[b_pos];
						
				if( b_pos >= n )
				{
					b_pos = 0;
					cont = true;
					break;
				}
					
				if( b < a )
				{
					smallest_fp_pos =  draws->at(i);
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
				for(INT j = 0; j<n; j++)
				{
		      			unsigned char a = X[a_pos];
					unsigned char b = X[b_pos];
							
					if( a_pos >= n )
					{
						a_pos = 0;
						cont = true;
						break;
					}
					
					if(  b < a )
					{
						smallest_fp_pos =  draws->at(i);
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
				for(INT j = 0; j<n; j++)
				{
		      			unsigned char a = X[a_pos];
					unsigned char b = X[b_pos];
							
					if( b_pos >= n || b < a )
					{
						smallest_fp_pos =  draws->at(i);
						break;
					}
					else if ( b > a )
						break;
					a_pos++;
					b_pos++;
				}
			}	
		}	
        
        }
        delete( draws );
   	return smallest_fp_pos;
}


/* Booth's O(n)-time algorithm -- slightly adapted for efficiency */
INT minlexrot( unsigned char * X, INT *f, INT n)
{  
	INT n_d = n<<1;
  	for(INT i = 0; i < n_d; ++i)	f[i] = (INT) -1;

  	INT k = 0;
  	for (INT j = 1; j < n_d; ++j)
  	{
                unsigned char sj = X[j%n];
                INT i = f[j - k - 1];
                while (i != (INT)-1 && sj != X[(k + i + 1)%n])
                {
                        if (sj < X[(k + i + 1)%n])        k = j - i - 1;
                        i = f[i];
                }
				
                if (i == (INT) - 1 && sj != X[(k + i + 1)%n])
                {
                        if (sj < X[(k+i+1)%n])    k = j;
                        f[j - k] = -1;
                }
                else
                        f[j - k] = i + 1;
   	}
   	return k;
}


/* Computes the length of lcp of two suffixes of two strings */
INT lcp ( unsigned char *  x, INT M, unsigned char * y, INT l, INT a_size, INT w_size )
{
	INT xx = a_size;
	if ( M >= xx ) return 0;
	INT yy = w_size;
	if ( l >= yy ) return 0;

	INT i = 0;
	while ( ( M + i < xx ) && ( l + i < yy ) )
	{
		if ( x[M+i] != y[l+i] )	break;
		i++;
	}
	return i;
}

/* Searching a list of strings using LCP from "Algorithms on Strings" by Crochemore et al. Algorithm takes O(m + log n), where n is the list size and m the length of pattern */
pair<INT,INT> pattern_matching ( unsigned char * w, unsigned char * a, vector<INT> * SA, vector<INT> * LCP, rmq_succinct_sct<> &rmq, INT n, INT w_size, INT a_size )
{

	INT m = w_size; //length of pattern
	INT N = a_size; //length of string
	INT d = -1;
	INT ld = 0;
	INT f = n;
	INT lf = 0;

	pair<INT,INT> interval;

	while ( d + 1 < f )
	{
		
		INT i = (d + f)/2;
		
		/* lcp(i,f) */
		INT lcpif;
		
		if( f == n )
			lcpif = 0;
		else lcpif = LCP->at(rmq ( i + 1, f ) );
			
		/* lcp(d,i) */
		INT lcpdi;
		
		if( i == n )
			lcpdi = 0;
		else lcpdi = LCP->at(rmq ( d + 1, i ) );
	
		if ( ( ld <= lcpif ) && ( lcpif < lf ) )
		{
			d = i;
			ld = lcpif;
		}
		else if ( ( ld <= lf ) && ( lf < lcpif ) ) 	f = i;
		else if ( ( lf <= lcpdi ) && ( lcpdi < ld ) )
		{
			f = i;
			lf = lcpdi;
		}
		else if ( ( lf <= ld ) && ( ld < lcpdi ) )	d = i;
		else
		{
			INT l = std::max (ld, lf);
			l = l + lcp ( a, SA->at( i ) + l, w, l, a_size, w_size );
			if ( l == m ) //lower bound is found, let's find the upper bound
		        {
				INT e = i;
				while ( d + 1 < e )
				{
					INT j = (d + e)/2;

					/* lcp(j,e) */
					INT lcpje;
				
					if( e == n )
						lcpje = 0;
					else lcpje = LCP->at( rmq ( j + 1, e ) );
					
					if ( lcpje < m ) 	d = j;
					else 			e = j;
				}

				/* lcp(d,e) */
				INT lcpde;
				
				if( e == n )
					lcpde = 0;
				else lcpde = LCP->at( rmq ( d + 1, e ) );
				
				if ( lcpde >= m )	d = std::max (d-1,( INT ) -1 );

				e = i;
				while ( e + 1 < f )
				{
					INT j = (e + f)/2;

					/* lcp(e,j) */
					INT lcpej;
					
					if( j == n )
						lcpej = 0;
					else lcpej = LCP->at( rmq ( e + 1, j ) );
					
					if ( lcpej < m ) 	f = j;
					else 			e = j;
				}

				/* lcp(e,f) */
				INT lcpef;
				
				if( f == n )
					lcpef = 0;
				else lcpef = LCP->at(rmq ( e + 1, f ) );
				
				if ( lcpef >= m )	f = std::min (f+1,n);

				interval.first = d + 1;
				interval.second = f - 1;
				return interval;


			}
			else if ( ( l == N - SA->at( i ) ) || ( ( SA->at( i ) + l < N ) && ( l != m ) && ( a[SA->at( i )+l] < w[l] ) ) )
			{
				d = i;
				ld = l;
			}
			else
			{
				f = i;
				lf = l;
			}
		}
	}
	

	interval.first = d + 1;
	interval.second = f - 1;
	return interval;
}

/* Computes the length of lcs of two suffixes of two strings */
INT lcs ( unsigned char *  x, INT M, unsigned char *  y, INT l, INT m )
{
	if ( M < 0 ) return 0;
	INT yy = m;
	if ( l >= yy ) return 0;

	INT i = 0;
	while ( ( M - i >= 0 ) && ( l + i < yy ) )
	{
		if ( x[M-i] != y[l+i] )	break;
		i++;
	}
	return i;
}


pair<INT,INT> rev_pattern_matching (unsigned char *  w, unsigned char *  a, vector<INT> * SA, vector<INT> * LCP, rmq_succinct_sct<> &rmq, INT n, INT w_size, INT a_size )
{
	
	
	INT m = w_size; //length of pattern
	INT N = a_size; //length of string
	INT d = -1;
	INT ld = 0;
	INT f = n;
	INT lf = 0;

	pair<INT,INT> interval;

	while ( d + 1 < f )
	{
		INT i = (d + f)/2;
		INT revSA = N - 1 - SA->at(i);
		//std::unordered_map<pair<INT,INT>, INT, boost::hash<pair<INT,INT> >>::iterator it;

		/* lcp(i,f) */
		INT lcpif;
		
		if( f == n )
			lcpif = 0;
		else lcpif = LCP->at( rmq ( i + 1, f ) );
		
		/* lcp(d,i) */
		INT lcpdi;
		//it = rmq.find(make_pair(d+1, i));
		
		if( i == n )
			lcpdi = 0;
		else lcpdi = LCP->at( rmq ( d + 1, i ) );
		
	
		if ( ( ld <= lcpif ) && ( lcpif < lf ) )
		{
			d = i;
			ld = lcpif;
		}
		else if ( ( ld <= lf ) && ( lf < lcpif ) ) 	f = i;
		else if ( ( lf <= lcpdi ) && ( lcpdi < ld ) )
		{
			f = i;
			lf = lcpdi;
		}
		else if ( ( lf <= ld ) && ( ld < lcpdi ) )	d = i;
		else
		{
			INT l = std::max (ld, lf);
			
			// avoid the function call if revSA-1<0 or l>=w.size() by changing lcs?
			l = l + lcs ( a, revSA - l, w, l, m );
			if ( l == m ) //lower bound is found, let's find the upper bound
		    	{
				INT e = i;
				while ( d + 1 < e )
				{
					INT j = (d + e)/2;

					/* lcp(j,e) */
					INT lcpje;
					
					
					if( e == n )
						lcpje = 0;
					else lcpje = LCP->at( rmq ( j + 1, e ) );
					
					if ( lcpje < m ) 	d = j;
					else 			e = j;
				}

				/* lcp(d,e) */
				INT lcpde;
				
				
				if( e == n )
					lcpde = 0;
				else lcpde = LCP->at( rmq ( d + 1, e ) );
				
				if ( lcpde >= m )	d = std::max (d-1,( INT ) -1 );
			
				e = i;
				while ( e + 1 < f )
				{
					INT j = (e + f)/2;

					/* lcp(e,j) */
					INT lcpej;
					
					
					if( j == n )
						lcpej = 0;
					else lcpej = LCP->at( rmq ( e + 1, j ) );
					
					if ( lcpej < m ) 	f = j;
					else 			e = j;
				}

				/* lcp(e,f) */
				INT lcpef;
					
				if( f == n )
					lcpef = 0;
				else lcpef = LCP->at( rmq ( e + 1, f ) );
				
				if ( lcpef >= m )	f = std::min (f+1,n);

				interval.first = d + 1;
				interval.second = f - 1;
				return interval;


			}
			else if ( ( l == N - SA->at(i) ) || ( ( revSA - l >= 0 ) && ( l != m ) && ( a[revSA - l] < w[l] ) ) )
			{
				d = i;
				ld = l;
			}
			else
			{
				f = i;
				lf = l;
			}

		}
	}
	

	interval.first = d + 1;
	interval.second = f - 1;
	
	return interval;
}


INT query(char * arg3, unsigned char * text_string, string output_filename, INT text_size, vector<INT> * LSA, vector<INT> * LLCP, vector<INT> * RSA, vector<INT> * RLCP, rmq_succinct_sct<> &lrmq, rmq_succinct_sct<> &rrmq, INT g, INT ell, INT power, INT k )
{
	INT num_seqs = 0;           // the total number of patterns considered
	INT max_len_pattern = 0;
	INT ALLOC_SIZE = 180224;
	INT seq_len = 0;
	INT max_alloc_seq_len = 0;
	INT max_alloc_seqs = 0;
	unsigned char ** patterns = NULL;
	
	unsigned char c = 0;
	// Input patterns
 	ifstream is_patterns;
 	is_patterns.open (arg3, ios::in | ios::binary);
	
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
				INT index = text_size-1-LSA->at(t);
				INT jj = j;		//this is the index of the anchor in the pattern
				index++; 	jj++;	//jump the index of the anchor and start looking on the right
				while ( ( jj < pattern_size ) && ( index < text_size ) && ( text_string[index] == patterns[i][jj] ) )
				{
					index++; jj++;
				}
				if ( jj == pattern_size ) //we have matched the pattern completely
				{ 
					if ( index == text_size - 1 )	
						pattern_output<< patterns[i] <<" found at position "<< index - pattern_size + 1 << " of the text"<<endl;					
					else			
						
						pattern_output<< patterns[i] <<" found at position "<<  index - pattern_size << " of the text"<<endl;
					hits++;
				}
			}
			
		}
	
   	}
   	 	
	for( INT i = 0; i < num_seqs; i ++ )
        	free (patterns[i]);
        free (patterns);
        
   	free( left_pattern );
  	free( first_window );
  	free( right_pattern );
  	delete [] f ;
  	
   	return hits;
 	
}
