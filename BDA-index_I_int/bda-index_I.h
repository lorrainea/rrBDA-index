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


 
struct Minimizer
 {
   INT               pos;
   INT		     min;
   INT               start_window;
   INT		     end_window;
   
 }; 
 
 
//INT ssa_lcp(rklce::rk_lce lce, vector<INT> * SSA, vector<INT> * LCP );
INT bd_anchors( unsigned char * seq, INT pos, INT ell, INT k, unordered_set<INT> &anchors, INT * rank);
INT red_minlexrot( string &X, INT n, INT r);
//INT ssa_rk_lce(rklce::rk_lce lce, vector<INT> * ssa_list,  char * sa_index_name, char * lcp_index_name, vector<INT> * SSA, vector<INT> * LCP );
INT ssa(string seq_filename, vector<INT> * ssa_list , char * sa_index_name, char * lcp_index_name, vector<INT> * final_ssa, vector<INT> * final_lcp );
pair<INT,INT> rev_pattern_matching ( string & w, string & a, vector<INT> * SA, vector<INT> * LCP, rmq_succinct_sct<> &rmq, INT n );
pair<INT,INT> pattern_matching ( string & w, string & a, vector<INT> * SA, vector<INT> * LCP, rmq_succinct_sct<> &rmq, INT n );

