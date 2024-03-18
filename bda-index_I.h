 #include <cstdio>
#include <cstdlib>
#include <unordered_set>
#include <string>
#include <sstream>

#ifdef _USE_64
typedef int64_t INT;
#endif

#ifdef _USE_32
typedef int32_t INT;
#endif

using namespace std;


 
struct Minimizer
 {
   INT               pos;
   INT		     min;
   INT               start_window;
   INT		     end_window;
   
 }; 
 
 
INT ssa_lcp(auto lce, vector<INT> * SSA, vector<INT> * LCP );
INT bd_anchors(  unsigned char * seq, string filename, INT ell, INT k, unordered_set<INT> &anchors, INT option );
INT red_minlexrot( string &X, INT *f, INT n, INT r );
//pair<INT,INT> rev_pattern_matching ( string & w, string & a, INT * SA, INT * LCP, rmq_succinct_sct<> &rmq, INT n );
//pair<INT,INT> pattern_matching ( string & w, string & a, INT * SA, INT * LCP, rmq_succinct_sct<> &rmq, INT n );

