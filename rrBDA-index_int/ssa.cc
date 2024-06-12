/**
    Sparse Suffix and LCP Array: Simple, Direct, Small, and Fast
    Copyright (C) 2023 Lorraine A. K. Ayad, Grigorios Loukides, 
    Solon P. Pissis and Hilde Verbeek

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

#include <iostream>
#include <cstring>
#include <string>
#include <fstream>
#include <stack>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <vector>
#include <cmath>
#include <sys/time.h>
#include <numeric>
#include <sstream>
#include "krfp.h"
#include "rrbda-index_int.h"
#include "unordered_dense.h"

#define IPS4 false

#if IPS4 == true
#include "../include/ips4o.hpp"
#endif

#define DEBUG false
#define THRESHOLD 1500000

using namespace std;

struct SSA
{
	INT lcp;
	vector<INT> L;
	SSA() : lcp(0), L() {}
};

double gettime( void )
{
	struct timeval ttime;
	gettimeofday( &ttime , 0 );
	return ttime.tv_sec + ttime.tv_usec * 0.000001;
};

double prep_total;
double hash_total;
double gr_total;
double order_total;
double sort_total;

/* Comparator */
auto compare(unsigned char * sequence, vector<INT> * A, INT lcp )
{
	return [sequence, A, lcp](INT a, INT b) 
	{
		return sequence[ A->at(a)+lcp ] < sequence[ A->at(b)+lcp ];
	};
}

/* Compute the KR fingerprint of sequence[ssa..ssa+l-1] using the FP table -- Time is O(min(l,n/s)), where s is the size of the FP table */
uint64_t fingerprint( INT ssa, uint64_t * FP, INT fp_len, INT l, unsigned char * sequence, INT text_size, uint64_t power )
{
	uint64_t fp = 0;
	INT ssa_end = (text_size >= ssa+l) ? ssa+l : text_size; //this is the end of the substring we are interested in PLUS 1
	bool subtract_slow = (text_size >= ssa+l) ? false : true;
	
	if( l > fp_len )
	{
		uint64_t fp_short = 0;
		if( ssa == 0 ) // then there is not short fragment before the substring
		{
			; // do nothing
		} 
		else if( ssa % fp_len  !=  0 ) // we are in-between two stored prefixes
		{
			uint64_t prefix = ssa / fp_len; // this is the prefix we are in but we have something stored on the left
			uint64_t start = 0;
			
			if( prefix != 0 )
			{
				fp_short = FP[prefix - 1];
				start = prefix * fp_len;
			}
			
			for(INT i = start; i<ssa; i++)	fp_short = karp_rabin_hashing::concat( fp_short, sequence[i], 1);

		}
		else 	// we have the fp_short stored and we read it from FP
		{	
			uint64_t prefix = ssa / fp_len;
			fp_short = FP[prefix - 1];
		}	

		uint64_t fp_long = 0;

		if( ssa_end % fp_len  != 0 )
		{
                        uint64_t prefix = ssa_end / fp_len; // this is the prefix we are in but we have something stored on the left
                        uint64_t start = 0;

                        if( prefix != 0 )
                        {
                        	fp_long = FP[prefix - 1];
                        	start = prefix * fp_len;
                        }
                        for(INT i = start; i< ssa_end; i++)	fp_long = karp_rabin_hashing::concat( fp_long, sequence[i] , 1 );
                }
                else
                { 
                	uint64_t prefix = ssa_end / fp_len;
                	fp_long = FP[prefix - 1];
                }

		if( subtract_slow == false )
               	 fp = karp_rabin_hashing::subtract_fast(fp_long, fp_short, power);
                else fp = karp_rabin_hashing::subtract(fp_long, fp_short, ssa_end - ssa);

        }
        else 
        {
        	for(INT i=ssa; i< ssa_end; ++i)	fp =  karp_rabin_hashing::concat( fp, sequence[i], 1 );
	}

	return fp;
}

/* Extend the prefixes of grouped suffixes by length l and re-group the computed KR fingerprints -- Time is O(b.min(l,n/s)), where s is the size of the FP table */
INT group( vector<SSA> &B, vector<INT> * A, uint64_t * FP, INT fp_len, INT l, unsigned char * sequence, INT text_size, INT b, INT &m, INT &z, uint64_t hash_variable )
{
    	vector<SSA> * B_prime = new vector<SSA>();
	(*B_prime).reserve(b); 

	uint64_t power = karp_rabin_hashing::pow_mod_mersenne(hash_variable, l, 61);
	
	const auto Bsz = B.size();
	
	for(INT i = 0; i<Bsz; ++i )
	{
		const INT s = (B)[i].L.size();
		INT k = 0;	
		vector<vector<INT>> vec;
		vector<INT> tmp;

		if( s <= (const INT)z )
		{
			double start = gettime();
			vector<pair<uint64_t,INT> > vec_to_sort;
			for(auto it=(B)[i].L.begin();it!=(B)[i].L.end(); ++it)
			{
				if ( (*A)[(*it)]+(B)[i].lcp + l > text_size ) { tmp.push_back((*it)); continue; }
				uint64_t fp = fingerprint( (*A)[(*it)]+(B)[i].lcp, FP, fp_len, l, sequence, text_size, power );
				vec_to_sort.push_back( make_pair(fp,*it) );
			}
			#if IPS4 == true
				ips4o::sort(vec_to_sort.begin(),vec_to_sort.end());
			#else
				sort(vec_to_sort.begin(),vec_to_sort.end());
			#endif 

			const auto vsz=vec_to_sort.size();
			for(INT i=0;i<vsz;++i)
			{
				if (i == 0 || (i>0 && vec_to_sort[i].first != vec_to_sort[i - 1].first)) //if we have new group 
				{
					vector<INT> new_vec;		
					new_vec.push_back(vec_to_sort[i].second);
					vec.push_back(new_vec);	 //adds into vec a new vector
					k++;
				}
				else	vec.back().push_back(vec_to_sort[i].second); //adds the id into the last added vector in vec
			}
			double end = gettime();
			sort_total += end - start;
		}
		else
		{	
			double start = gettime();
    			auto groups = ankerl::unordered_dense::map<uint64_t, INT >();
			for(auto it=(B)[i].L.begin();it!=(B)[i].L.end(); ++it)
			{
				if ( (*A)[(*it)]+(B)[i].lcp + l > text_size ) { tmp.push_back((*it)); continue; }
				uint64_t fp = fingerprint( (*A)[(*it)]+(B)[i].lcp, FP, fp_len, l, sequence, text_size, power );
				auto itx = groups.find(fp);
				if(itx == groups.end())
				{
					groups[fp] = k;
					vector<INT> v;
					v.push_back(*it);
					vec.push_back(v);
					k++;
				}
				else	vec[itx->second].push_back(*it);
			}
			groups.clear();

			double end = gettime();
			hash_total += end - start;
		}

		double start = gettime();
		
	    	vector<INT>().swap(B[i].L);

		for( INT j = 0; j < k; j++ )
		{	
			const auto itsz=vec[j].size();

			if( itsz == s )
			{
				(B)[i].lcp += l;
				for(const auto& value: vec[j])	(B)[i].L.push_back(value);				

			}
			else if( itsz >= 2 )
			{
				m++; 
				SSA new_ssa;

				for(const auto& value: vec[j])	new_ssa.L.push_back(value);				
				
				new_ssa.lcp = (B)[i].lcp + l;
				B_prime->push_back( new_ssa );
				(B)[i].L.push_back(m);

				A->push_back( (*A)[ vec[j][0] ] );											
			}
			else if ( itsz == 1 )
			{
				(B)[i].L.push_back( vec[j][0] );
			}
		}
		for (auto& it : tmp) (B)[i].L.push_back( it );
		vector<INT>().swap(tmp);
		vector<vector<INT>>().swap(vec);

		double end = gettime();
		gr_total += end - start;
	}

	double start = gettime();
	B.insert(std::end(B), std::begin(*B_prime), std::end(*B_prime));
	delete( B_prime);
	double end = gettime();
	gr_total += end - start;

	return 0;
}

/* Sort the final group members and infer the SSA and SLCP array -- Time is O(b log b) */
INT order( vector<INT> * final_ssa, vector<INT> * final_lcp, vector<SSA> &B, vector<INT> * A, unsigned char * sequence, INT text_size, INT b )
{

	const INT Bsz=B.size();
	for(INT i = 0; i<Bsz; i++)
	{	
		#if IPS4 == true
			ips4o::sort((B)[i].L.begin(), (B)[i].L.end(), compare(sequence,A,(B)[i].lcp));
		#else
			sort((B)[i].L.begin(), (B)[i].L.end(), compare(sequence,A,(B)[i].lcp));
		#endif		
	}
	stack<pair<INT,INT>> S; 
	
	S.push( make_pair<INT, INT>((INT)b, 0) );  //b is the correct first index, not b+1
	
	INT l = 0;
	INT mymax=numeric_limits<INT>::max();
	while( !S.empty() )
	{
		INT i = S.top().first;
		INT l_prime = S.top().second;
		S.pop();
		
		if( l_prime < l )
			l = l_prime;

		if(i>=b) //it is not one of the initial groups
		{
			INT lcp = (B)[i-b].lcp;
			auto myl=(B)[i-b].L;
			for(vector<INT>::reverse_iterator it=myl.rbegin();it!=myl.rend();++it)
			{				  					
				S.push(make_pair<INT, INT>( (INT) *it, (INT) lcp ));
			}
		}
		else
		{
			final_ssa->push_back( (*A)[i] );
			final_lcp->push_back( l );
			l = mymax;
		}	
	}

	return 0;
}


INT ssa(unsigned char * sequence, vector<INT> * ssa_list , string sa_index_name, string lcp_index_name, vector<INT> * final_ssa, vector<INT> * final_lcp, uint64_t hash_variable )
{
	INT z = THRESHOLD;
	INT text_size = strlen( (char*) sequence);

	INT b = ssa_list->size();
	cout<<"Number of suffixes b = " << b << endl;
	
	INT s = 2*b;
	if ( s > text_size ) s = text_size;
	INT fp_len = text_size / s;
	if ( (text_size - fp_len * s) > text_size/s ) 
		s = text_size/fp_len; 
	cout<<"Block length = "<<fp_len<<endl;
	cout<<"Size s of FP table = " << s <<endl<<endl;
	
	// computing fingerprints
	uint64_t * FP =  ( uint64_t * ) calloc( s , sizeof( uint64_t ) );
	uint64_t fp = 0;
	
	std::chrono::steady_clock::time_point start_total = std::chrono::steady_clock::now();
	
	prep_total = 0;
	double start = gettime();
	cout<<"Preprocessing starts"<<endl;
	INT i = 0;
	for(uint64_t j = 0; j < s; ++j )
	{
		for (uint64_t k = 0; k < fp_len; k ++ )fp = karp_rabin_hashing::concat(fp, sequence[i+k], 1);
			FP[j]=fp;
		i += fp_len;
		if ( i + fp_len > text_size ) break;
	}
	cout<<"Preprocessing ends"<<endl<<endl;
	double end = gettime();
	prep_total = end - start;

    	vector<SSA> B;

	vector<INT> * A = new vector<INT>();
	
	vector<INT> * A_prime = new vector<INT>();
	vector<INT> * P = new vector<INT>();
	
	vector<INT> L;

	const auto ssa_list_sz=ssa_list->size();
	for(INT i = 0; i<ssa_list_sz; ++i )
	{	
		A->push_back( (*ssa_list)[i] );
		L.push_back(i);
	}
	
	SSA initial;
	initial.lcp = 0;
	
    	initial.L.assign(L.begin(), L.end()); 

	B.push_back( initial );

	A->push_back( (*ssa_list)[0] );
	INT m = ssa_list->size();

	hash_total = 0;
	gr_total = 0;
	
	INT c1 = 1;
	INT initial_l = 1ULL << static_cast<INT>(log2(c1*text_size/b));
	INT next_initial_l = initial_l * 2 - 1;
	//vector<INT> * final_ssa = new vector<INT>();
	//vector<INT> * final_lcp = new vector<INT>();
	
	cout<<"First run starts"<<endl;
	while( initial_l > 0 )
	{
		cout<< "Initial l: " << initial_l <<", nodes: "<< m <<endl;
		group( B, A, FP, fp_len, initial_l, sequence, text_size, b, m, z, hash_variable );
		initial_l=initial_l>>1;	
	}	
		
	order_total = 0;
	start = gettime();
	
	order( final_ssa, final_lcp, B, A, sequence, text_size, b);	
	end = gettime();
	order_total = end - start;
	cout<<"First run ends"<<endl<<endl;
	
	for(INT i = 0; i<b; i++)
	{
		
		if( (*final_lcp)[i] == next_initial_l || ( i < b-1 && (*final_lcp)[i+1] == next_initial_l ) )
		{
			P->push_back(i);
			A_prime->push_back((*final_ssa)[i]);
		}
	}
	
		
	vector<INT> * final_ssa_prime = new vector<INT>();
	vector<INT> * final_lcp_prime = new vector<INT>();
	
	if( P->size() > 0 )
	{
		cout<<"Second run starts"<<endl;
		
		INT l = 1ULL << static_cast<INT>(log2(text_size));
		vector<SSA>().swap(B);

		b = A_prime->size();
		m = A_prime->size();
		
		vector<INT> L(m); 
	    	std::iota(L.begin(), L.end(), 0); 

		SSA initial;
		initial.lcp = 0;
		initial.L.assign(L.begin(), L.end()); 
		
		B.push_back( initial );
		A_prime->push_back( (*A_prime)[0] );
		
		while( l > 0 )
		{	
			cout<< "l: " << l <<", nodes: "<< m <<endl;
			group( B, A_prime, FP, fp_len, l, sequence, text_size, b, m, z, hash_variable );
			l=l>>1;
		}	
		
		start = gettime();
		order( final_ssa_prime, final_lcp_prime, B, A_prime, sequence, text_size, b);
		end = gettime();
		order_total += end - start;
		
		const auto Psz=P->size();	
		for(INT i = 0; i<Psz; ++i)
		{
			(*final_ssa) [(*P)[i] ] = (*final_ssa_prime)[i];
			if( (*final_lcp)[ (*P)[i] ]== next_initial_l )
				(*final_lcp)[ (*P)[i] ] = (*final_lcp_prime)[i];
		}

		cout<<"Second run ends"<<endl;
	}
	
	free( FP );
	
	std::chrono::steady_clock::time_point end_total = std::chrono::steady_clock::now();

	ofstream output_ssa(sa_index_name);
	for(INT i = 0; i<final_ssa->size(); i++)
		output_ssa<<final_ssa->at(i)<<endl;
	
	ofstream output_lcp(lcp_index_name);
	for(INT i = 0; i<final_lcp->size(); i++)
		output_lcp<<final_lcp->at(i)<<endl;
	
	delete( final_lcp_prime );
	delete( final_ssa_prime );
	delete( A );
	delete( A_prime );
	delete( P );
	
	return 0;
}

