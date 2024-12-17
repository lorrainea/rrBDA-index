/**
    rrBDA-index_ext: Randomized Reduced Bi-directional Anchors (using external memory)
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
 
INT bd_anchors( unsigned char * seq, INT pos, INT ell, uint64_t k, unordered_set<INT> &anchors, uint64_t * rank, uint64_t power );
INT red_minlexrot( unsigned char * X, INT n, uint64_t r, uint64_t power );
INT compute_anchors(char * arg1, unordered_set<INT> &text_anchors, INT text_size, INT block, INT ell, INT k, uint64_t power);
INT compute_index( uint64_t hash, string index_name, INT text_size, INT g, unsigned char * text_string, unordered_set<INT> &text_anchors, INT * RSA, INT * RLCP, INT * LSA, INT * LLCP, rmq_succinct_sct<> &lrmq, rmq_succinct_sct<> &rrmq, char * arg0, char * arg1, INT ram_use );
INT query(char * arg3, unsigned char * text_string, string output_filename, INT text_size, INT * LSA, INT * LLCP, INT * RSA, INT * RLCP, rmq_succinct_sct<> &lrmq, rmq_succinct_sct<> &rrmq, INT g, INT ell, INT power, INT k );
pair<INT,INT> rev_pattern_matching ( unsigned char * w, unsigned char * a, INT * SA, INT * LCP, rmq_succinct_sct<> &rmq, INT n, INT w_size, INT a_size );
pair<INT,INT> pattern_matching ( unsigned char * w, unsigned char * a, INT * SA, INT * LCP, rmq_succinct_sct<> &rmq, INT n, INT w_size, INT a_size );

