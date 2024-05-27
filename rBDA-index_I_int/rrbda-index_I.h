/**
    rrBDA-index_I_int: Randomized Reduced Bi-directional Anchors (using internal memory)
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

#include <cstdio>
#include <cstdlib>
#include <unordered_set>
#include <string>
#include <sstream>
#include <sdsl/bit_vectors.hpp>                                   
#include <sdsl/rmq_support.hpp>

#ifdef _USE_64
typedef int64_t INT;
#endif

#ifdef _USE_32
typedef int32_t INT;
#endif


using namespace sdsl;
using namespace std;
 
INT bd_anchors( unsigned char * seq, INT pos, INT ell, INT k, unordered_set<INT> &anchors, INT * rank, INT power );
INT red_minlexrot( string &X, INT n, INT r, INT power );
INT ssa(string seq_filename, vector<INT> * ssa_list , string sa_index_name,string lcp_index_name, vector<INT> * final_ssa, vector<INT> * final_lcp, INT hash_variable );
pair<INT,INT> rev_pattern_matching ( string & w, string & a, vector<INT> * SA, vector<INT> * LCP, rmq_succinct_sct<> &rmq, INT n );
pair<INT,INT> pattern_matching ( string & w, string & a, vector<INT> * SA, vector<INT> * LCP, rmq_succinct_sct<> &rmq, INT n );

