
//#include "utils.h"
#include "stream.h"
#include "uint40.h"
#include <math.h>
#include "grid.h"
#include "rbda-index_I.h"
#include <sdsl/bit_vectors.hpp>                                   // include header for bit vectors
#include <sdsl/rmq_support.hpp>	
#include <sdsl/io.hpp>
#include "krfp.h"

#ifdef _USE_64
#include <divsufsort64.h>                                         // include header for suffix sort
#endif

#ifdef _USE_32
#include <divsufsort.h>                                       	  // include header for suffix sort
#endif

//unordered_set<INT> draws;
typedef grid_point point;
typedef grid_query query;

double vm, vm0, rss, rss0;

/* Sorting comparison */
bool sort_sa(const pair<INT,INT> &a,const pair<INT,INT> &b)
{
       return a.first<b.first;
}

void process_mem_usage(double& vm_usage, double& resident_set)
{
    vm_usage     = 0.0;
    resident_set = 0.0;

    // the two fields we want
    unsigned long vsize;
    long rss;
    {
        std::string ignore;
        std::ifstream ifs("/proc/self/stat", std::ios_base::in);
        ifs >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
                >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
                >> ignore >> ignore >> vsize >> rss;
    }

    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    vm_usage = vsize / 1024.0;
    resident_set = rss * page_size_kb;
}


int main(int argc, char **argv)
{
	unordered_set<char> alphabet;

	if( argc < 7 )
 	{
        	cout<<"Wrong arguments!\n";
 		cout<<"./rbda-index_I <text_file> <ell> <pattern_file> <block_size> <output_filename> <index_filename>\n";
 		exit(-1);
 	}
	
	// Input text file
 	ifstream is_text;
 	is_text.open (argv[1], ios::in | ios::binary);
 	
 	ifstream in_file(argv[1], ios::binary);
   	in_file.seekg(0, ios::end);
   	INT text_file_size = in_file.tellg();

	// Input ell
 	std::string str_ell(argv[2]);
 	
 	INT ell;
 	std::stringstream(str_ell)>>ell;
	
	// Input block size
 	std::string str_block(argv[4]);
 	
 	INT block;
 	std::stringstream(str_block)>>block;
 	
 	// Input patterns
 	ifstream is_patterns;
 	is_patterns.open (argv[3], ios::in | ios::binary);
 	
 	// Input output file
 	string output_filename = argv[5]; 
 	
 	// Input index file
 	string index_name = argv[6];
    	ifstream is_index;
 	is_index.open (argv[6], ios::in | ios::binary);

    	std::chrono::steady_clock::time_point  start_bd = std::chrono::steady_clock::now();
 	unordered_set<INT> text_anchors;
	
	
	string bd = index_name + ".bd";
 	
 	ifstream is_bd_anchors;
 	is_bd_anchors.open (bd, ios::in | ios::binary);
 	
	ifstream in_bd_anchors(bd, ios::binary);
 	in_bd_anchors.seekg (0, in_bd_anchors.end);
   	INT file_size = in_bd_anchors.tellg();
   	string bd_anchor = "";
   	INT bd_anchor_int = 0;
   	
   	char c = 0;
  	INT text_size = 0;
	for (INT i = 0; i < text_file_size; i++)
	{	
		is_text.read(reinterpret_cast<char*>(&c), 1);
		
		if( (unsigned char) c == '\n' )
			continue;
		else
		{
			alphabet.insert( (unsigned char) c );
			text_size++;
		}
		
	}
	is_text.close();
	
	if( text_size < ell )
	{
	
		fprintf( stderr, " Error: Window size (ell) cannot be larger than sequence length!\n");
		return ( 1 );
	}
	
	if( block < ell )
	{
		fprintf( stderr, " Error: Window size (ell) cannot be larger than the block size!\n");
		return ( 1 );
	}
	
	INT k  = ceil(4*log2(ell)/log2(alphabet.size()));
	if( ell - k - 1 < 0 )
		k = 2;
	
	
	INT hash = karp_rabin_hashing::init();
	INT power = karp_rabin_hashing::pow_mod_mersenne(hash, k, 61);
	
   	if( file_size > 0 )
	{
		// Read in from .bd file
	    	c = 0;
		for (INT i = 0; i < file_size; i++)
		{	
			is_bd_anchors.read(reinterpret_cast<char*>(&c), 1);
		
			if( (unsigned char) c == '\n' )
			{
				bd_anchor_int = stoi( bd_anchor );
				text_anchors.insert( bd_anchor_int );
				bd_anchor = "";
			}
			else bd_anchor += (unsigned char) c;
			
		}
		is_bd_anchors.close();
	}
	else
	{
		ifstream is_block;
	 	is_block.open (argv[1], ios::in | ios::binary);
		unsigned char * text_block = ( unsigned char * ) malloc (  ( block + 1 ) * sizeof ( unsigned char ) );
		INT * rank = ( INT * ) malloc( ( block  ) *  sizeof( INT ) );
		unsigned char * suffix_block = ( unsigned char * ) malloc (  ( ell  ) * sizeof ( unsigned char ) );
		
	 	c = 0;
	 	INT count = 0;
	 	INT pos = 0;
	 	
		for (INT i = 0; i < text_size; i++)
		{	
			is_block.read(reinterpret_cast<char*>(&c), 1);
			
			if( (unsigned char) c != '\n' )
			{
				text_block[count] = (unsigned char) c ;
				count++;
				if( count == block || i == text_size - 1 )
				{
					text_block[count] = '\0';
					
					bd_anchors( text_block, pos, ell, k, text_anchors, rank, power );
					
					memcpy( &suffix_block[0], &text_block[ block - ell + 1], ell -1 );
					memcpy( &text_block[0], &suffix_block[0], ell -1 );
					
					pos = pos + ( block - ell + 1 );
					count = ell - 1;
					
					
				}
			}
			
		}
		
		ofstream bd_output;
		bd_output.open(bd);
			
		for (auto &anchor : text_anchors)	
			bd_output<<anchor<<endl;
				
		bd_output.close();
		
		is_block.close();
		free( text_block );
		free( suffix_block );
		free( rank );	
		
	}
			
	INT g = text_anchors.size();
	INT n = text_size;
	
	std::chrono::steady_clock::time_point  end_bd = std::chrono::steady_clock::now();
	std::cout <<"bd construction took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_bd - start_bd).count() << "[ms]" << std::endl;
	cout<<"The text is of length "<< n << ", its alphabet size is "<< alphabet.size()<<", and it has "<<g<<" bd-anchors of order "<<ell<<endl;
	cout<<"The density is "<<(double) g / text_size<<endl;
	
	string text_string = "";
	ifstream is_full;
 	is_full.open (argv[1], ios::in | ios::binary);
 	
	c = 0;
	for (INT i = 0; i < text_size; i++)
	{	
		is_full.read(reinterpret_cast<char*>(&c), 1);
		
		if( (unsigned char) c == '\n'  )
			continue;
			
		else text_string.push_back( (unsigned char) c );
	}
	is_full.close();
	
	std::chrono::steady_clock::time_point  start_index = std::chrono::steady_clock::now();
	
	/* Constructing right and left compacted tries */
	
	vector<INT> * RSA = new vector<INT>();
	vector<INT> * RLCP = new vector<INT>();
	
	string sa_index_name = index_name + ".RSA";
 	
	ifstream in_RSA(sa_index_name, ios::binary);
	in_RSA.seekg (0, in_RSA.end);
	INT file_size_sa = in_RSA.tellg();
	
	ifstream is_RSA;
 	is_RSA.open (sa_index_name, ios::in | ios::binary);
	
	string lcp_index_name = index_name + ".RLCP";
	
	ifstream is_RLCP;
 	is_RLCP.open (lcp_index_name, ios::in | ios::binary);
 	
	ifstream in_RLCP(lcp_index_name, ios::binary);
	in_RLCP.seekg (0, in_RLCP.end);
	INT file_size_lcp = in_RLCP.tellg();
	
	vector<INT> * anchors_vector = new vector<INT>();
	for(auto &anchor: text_anchors)
		anchors_vector->push_back(anchor);
		
	text_anchors.clear();
	
	if( !(is_RSA) || !(is_RLCP )  )
	{
		ssa(argv[1], anchors_vector, sa_index_name, lcp_index_name, RSA, RLCP, hash );
	 
	} 	
	else 
	{
		if( file_size_sa > 0 )
		{
			
		   	string sa = "";
		   	INT sa_int = 0;
			c = 0;
			
			for (INT i = 0; i < file_size_sa; i++)
			{	
				is_RSA.read(reinterpret_cast<char*>(&c), 1);
			
				if( (unsigned char) c == '\n' )
				{
					sa_int = stoi( sa);
					RSA->push_back( sa_int );
					sa = "";
				}
				else sa += (unsigned char) c;
				
			}
			is_RSA.close();
		
		
		}
		
	
		if( file_size_lcp > 0 )
		{	
			
		   	string lcp = "";
		   	INT lcp_int = 0;
			c = 0;
			
			for (INT i = 0; i < file_size_lcp; i++)
			{	
				is_RLCP.read(reinterpret_cast<char*>(&c), 1);
				
				if( (unsigned char) c == '\n' )
				{
					lcp_int = stoi(lcp);
					RLCP->push_back(lcp_int);
					lcp = "";
				}
				else lcp += (unsigned char) c;
				
			}
			is_RLCP.close();
		}
	
	}
	
	cout<<"Right Compacted trie constructed "<<endl;	
	
	vector<INT> * LSA = new vector<INT>();
	vector<INT> * LLCP = new vector<INT>();
	
	string text_name = argv[1];
	string output_reverse = text_name + "_reverse";
		
  	/* We reverse the string for the left direction and also overwrite all other DSs */
  	std::ofstream output_r;
  	output_r.open (output_reverse);
  	reverse(text_string.begin(), text_string.end());
  	output_r << text_string;
    	output_r.close();
 
	sa_index_name = index_name + ".LSA";
	
	ifstream is_LSA;
 	is_LSA.open (sa_index_name, ios::in | ios::binary);
 	
	ifstream in_LSA(sa_index_name, ios::binary);
	in_LSA.seekg (0, in_LSA.end);
	file_size_sa = in_LSA.tellg();
	
	lcp_index_name = index_name +".LLCP" ;
	
	ifstream is_LLCP;
 	is_LLCP.open (lcp_index_name, ios::in | ios::binary);
	
	ifstream in_LLCP(lcp_index_name, ios::binary);
	in_LLCP.seekg (0, in_LLCP.end);
	file_size_lcp = in_LLCP.tellg();


	for(INT i = 0; i<anchors_vector->size(); i++)
	{
		anchors_vector->at(i) = ( text_size - 1 ) - anchors_vector->at(i);
	}
	

	if ( !(is_LSA) || !(is_LLCP) )
	{
		ssa(output_reverse, anchors_vector, sa_index_name, lcp_index_name, LSA, LLCP, hash );
	}
	
	else
	{
		if( file_size_sa  > 0 )
		{
		   	string sa = "";
		   	INT sa_int = 0;
			c = 0;
			
			for (INT i = 0; i < file_size_sa; i++)
			{	
				is_LSA.read(reinterpret_cast<char*>(&c), 1);
				
				if( (unsigned char) c == '\n')
				{
					sa_int = stoi( sa);
					LSA->push_back( sa_int );
					sa = "";
				}
				else sa += (unsigned char) c;
				
			}
			is_LSA.close();	
		}
	
		if( file_size_lcp  > 0 )
		{
		   	string lcp = "";
		   	INT lcp_int = 0;
			c = 0;
			
			for (INT i = 0; i < file_size_lcp; i++)
			{	
				is_LLCP.read(reinterpret_cast<char*>(&c), 1);
				
				if( (unsigned char) c == '\n' )
				{
					lcp_int = stoi( lcp );
					LLCP->push_back(lcp_int);
					lcp = "";
				}
				else lcp += (unsigned char) c;
				
			}
			is_LLCP.close();	
		}
	}
	
	cout<<"Left Compacted trie constructed "<<endl;
	
	delete( anchors_vector );
	/* After constructing the tries these DSs over the whole string are not needed anymore, our data structure must be of size O(g) */

  	/* The following RMQ data structures are used for spelling pattern over the LSA and RSA */
  
	string rmq_left_suffix = index_name +".lrmq";
  		
	ifstream in_rmq_left(rmq_left_suffix, ios::binary);
  	rmq_succinct_sct<> lrmq;
  	
  	
  	if( in_rmq_left )
  	{
  	
  		load_from_file(lrmq, rmq_left_suffix); 
	}
  	else
  	{
	  	int_vector<> llcp( g , 0 ); // create a vector of length n and initialize it with 0s

		
		for ( INT i = 0; i < g; i ++ )
		{
			llcp[i] = LLCP->at(i);
			
		}

		util::assign(lrmq, rmq_succinct_sct<>(&llcp));
		
		util::clear(llcp);

		store_to_file(lrmq, rmq_left_suffix);
	}
	
	cout<<"Left RMQ DS constructed "<<endl;
	string rmq_right_suffix = index_name+ ".rrmq";
	
	ifstream in_rmq_right(rmq_right_suffix, ios::binary);
  	  
  	int_vector<> rlcp( g , 0 ); // create a vector of length n and initialize it with 0s
  	rmq_succinct_sct<> rrmq;
  	
  	if( in_rmq_right )
  	{

  		load_from_file(rrmq, rmq_right_suffix);
  	}
  	else
  	{
		int_vector<> rlcp( g , 0 ); // create a vector of length n and initialize it with 0s

		for ( INT i = 0; i < g; i ++ )
		{
			rlcp[i] = RLCP->at(i);
		}
		
		util::assign(rrmq, rmq_succinct_sct<>(&rlcp));
		
		util::clear(rlcp);
		 
	  	
	  	
	  	store_to_file(rrmq, rmq_right_suffix);
	}	
	 
	cout<<"Right RMQ DS constructed "<<endl;
  	//Forming coordinate points using LSA and RSA

	string grid_suffix =  index_name+".grid";
  		
  	ifstream in_grid(grid_suffix, ios::binary);
  	in_grid.seekg (0, in_grid.end);
  	file_size = in_grid.tellg();
  	
  	ifstream is_grid;
 	is_grid.open (grid_suffix, ios::in | ios::binary);
   	
   	std::vector<point> points;
   	grid construct;
   	
 
   	if( file_size > 0 )
  	{  	
  		load_from_file(construct, grid_suffix);
  	}
  	else
  	{
	  	vector<pair<INT,INT>> l;
		vector<pair<INT,INT>> r;
		for ( INT i = 0; i < g; i++ )
	  	{
	  		l.push_back( make_pair( n-1-LSA->at(i), i) );
	  		r.push_back( make_pair( RSA->at(i), i ) );
	 	}
	 	
	 	sort(l.begin(),l.end(),sort_sa);
		sort(r.begin(),r.end(),sort_sa);
		
	  	for ( INT i = 0; i < g; i++ )
	  	{
	 
			point to_insert;
			to_insert.row = l.at(i).second+1;
			to_insert.col = r.at(i).second+1;
			
			to_insert.level = 0;
			to_insert.label = l.at(i).first;
			points.push_back(to_insert); 
	  	}
	  		
		construct.build( points, 0 );
		
		store_to_file( construct, grid_suffix );	
	}
	
	
	cout<<"The grid is constructed"<<endl;  
  	cout<<"The whole index is constructed"<<endl;
  	
	std::chrono::steady_clock::time_point  end_index = std::chrono::steady_clock::now();
	std::cout <<"Index took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_index- start_index).count() << "[ms]" << std::endl;


	std::chrono::steady_clock::time_point  begin_pt = std::chrono::steady_clock::now();
    	reverse(text_string.begin(), text_string.end()); 				//I re-reverse to take the original string
  

	vector<vector<unsigned char> > all_patterns;
    	vector<unsigned char> pattern;
    	c = 0;
    	while (is_patterns.read(reinterpret_cast<char*>(&c), 1))
    	{
        	if(c == '\n')
        	{
  			if(pattern.empty())	break;
  			all_patterns.push_back(pattern);
  			pattern.clear();
        	}
        	else	pattern.push_back((unsigned char)c);
    	}
    	is_patterns.close();
    	pattern.clear();

	vector<string> new_all_pat;
	for(auto &it_pat : all_patterns)	new_all_pat.push_back(string(it_pat.begin(), it_pat.end()));
	all_patterns.clear();
  	
	std::chrono::steady_clock::time_point  begin = std::chrono::steady_clock::now();
	
	ofstream pattern_output;
	pattern_output.open(output_filename);
	for(auto &pattern : new_all_pat)
   	{
  		if ( pattern.size() < ell )
  		{
  			pattern_output<<"Pattern skipped: its length is less than ell!\n";
  			continue;
  		}
  		
		string first_window = pattern.substr(0, ell);
		INT j = red_minlexrot( first_window, ell, k, power );
		string left_pattern = pattern.substr(0, j+1);
	  	reverse(left_pattern.begin(), left_pattern.end());
	  		
	  	
		pair<INT,INT> left_interval = rev_pattern_matching ( left_pattern, text_string, LSA, LLCP, lrmq, g );  
	  		
		string right_pattern = pattern.substr(j, pattern.size()-j);
		pair<INT,INT> right_interval = pattern_matching ( right_pattern, text_string, RSA, RLCP, rrmq, g );
	  	if ( left_interval.first <= left_interval.second  && right_interval.first <= right_interval.second )
	  	{
	  		//Finding rectangle containing bd-anchors in grid
	  		grid_query rectangle;
	  		rectangle.row1 = left_interval.first+1;
	  		rectangle.row2 = left_interval.second+1;
	  		rectangle.col1 = right_interval.first+1;
	  		rectangle.col2 = right_interval.second+1;
	  		
				
			vector<long unsigned int> result;
			construct.search_2d(rectangle, result);
			
			for(INT i = 0; i<result.size(); i++)
			{
				pattern_output<<pattern<<" found at position "<<RSA->at(result.at(i)-1)-j<<" of the text"<<endl;	
			}
		
  		}
	 	else	pattern_output<<"No occurrences found!\n";
			
		
  		
 	}
	pattern_output.close();
 	
	std::chrono::steady_clock::time_point  end_pt = std::chrono::steady_clock::now();
	std::cout <<"Pattern matching took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_pt - begin_pt).count() << "[ms]" << std::endl;
  	delete ( RSA );
  	delete ( RLCP );
  	delete ( LSA );
  	delete ( LLCP );

	return 0;
}
