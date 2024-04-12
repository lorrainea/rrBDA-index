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
pair<INT,INT> rev_pattern_matching ( string & w, string & a, INT * SA, INT * LCP, rmq_succinct_sct<> &rmq, INT n );
pair<INT,INT> pattern_matching ( string & w, string & a, INT * SA, INT * LCP, rmq_succinct_sct<> &rmq, INT n );

