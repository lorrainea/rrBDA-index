#include <unordered_set>
#include "utils.h"
#include "stream.h"
#include "uint40.h"
#include <math.h>
#include "bda-index_I.h"
#include "krfp.h"

using namespace std;
using namespace sdsl;

INT red_minlexrot( string &X, INT n, INT k, INT power )
{  
  	INT fp = 0;
    	vector<INT> * draws = new vector<INT>();
	
	for(INT j = 0; j<k; j++)
        	fp =  karp_rabin_hashing::concat( fp, X[j] , 1 );

	draws->push_back(0);
	INT smallest_fp = fp;
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
       		smallest_fp_pos = draws->at(0);
       		
        	for(INT i = 1; i<draws->size(); i++ )
        	{
        		INT a_pos = smallest_fp_pos+k;
        		INT b_pos = draws->at(i)+k;
        		
        		for(INT j = k; j<n; j++)
        		{
 
        			char a = X[a_pos];
        			char b = X[b_pos];

				if( a_pos >=n )
				{
					a_pos = 0;
					break;
				}
					
				if( b_pos >= n )
				{
					b_pos = 0;
					break;
				}
				
        			if( b < a )
        			{
        				smallest_fp_pos =  draws->at(i);
        				
        				break;
        			}
        			
        			a_pos++;
        			b_pos++;
    
        		}
        		
        		for(INT j = 0; j<n; j++)
        		{
       
        			char a = X[a_pos];
        			char b = X[b_pos];
        			
        			if( a_pos >=n )
				{
					break;
				}
					
				if( b_pos >= n )
				{
					smallest_fp_pos =  draws->at(i);
					
					break;
				}
        			
        			if( b < a )
        			{
        				smallest_fp_pos =  draws->at(i);
        				
        				break;
        			}
        			
        			a_pos++;
        			b_pos++;
    
        		}
        		
        		
        	}
        
        }
        
        delete( draws );
   	return smallest_fp_pos;
}


/* Booth's O(n)-time algorithm -- slightly adapted for efficiency */
INT minlexrot( string &X, INT *f, INT n)
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
INT lcp ( string & x, INT M, string & y, INT l )
{
	INT xx = x.size();
	if ( M >= xx ) return 0;
	INT yy = y.size();
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
pair<INT,INT> pattern_matching ( string & w, string & a, vector<INT> * SA, vector<INT> * LCP, rmq_succinct_sct<> &rmq, INT n )
{

	INT m = w.size(); //length of pattern
	INT N = a.size(); //length of string
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
			l = l + lcp ( a, SA->at(i) + l, w, l );
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
			else if ( ( l == N - SA->at(i) ) || ( ( SA->at(i) + l < N ) && ( l != m ) && ( a[SA->at(i)+l] < w[l] ) ) )
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
INT lcs ( string & x, INT M, string & y, INT l )
{
	if ( M < 0 ) return 0;
	INT yy = y.size();
	if ( l >= yy ) return 0;

	INT i = 0;
	while ( ( M - i >= 0 ) && ( l + i < yy ) )
	{
		if ( x[M-i] != y[l+i] )	break;
		i++;
	}
	return i;
}


/* Searching a list of strings using LCP from "Algorithms on Strings" by Crochemore et al. Algorithm takes O(m + log n), where n is the list size and m the length of pattern */
pair<INT,INT> rev_pattern_matching ( string & w, string & a, vector<INT> * SA, vector<INT> * LCP, rmq_succinct_sct<> &rmq, INT n )
{
	
	
	INT m = w.size(); //length of pattern
	INT N = a.size(); //length of string
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
			l = l + lcs ( a, revSA - l, w, l );
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
