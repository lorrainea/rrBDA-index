/*
 *  This file is part of rk-lce.
 *  Copyright (c) by
 *  Nicola Prezza <nicolapr@gmail.com>
 *
 *   rk-lce is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.

 *   rk-lce is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details (<http://www.gnu.org/licenses/>).
 */

/*
 * Modification: This code is modified from the original to read the input 
 * from text files and to compute and output the sparse LCP array from the 
 * LCE data structure. These modifications were done by the authors of SSA.
 */


#include "includes.hpp"
#include "rk_lce.hpp"
#include "bda-index_I.h"

using namespace std;
using namespace rklce;

string output_file = string();
string input_text = string();
string input_pos = string();

rk_lce lce_structure(unsigned char * seq)
{
	string input_text = reinterpret_cast<char*>(seq);

	cout << "Building LCP data structure and allocating memory for SA ... " << endl;
	
	auto lce = rk_lce(input_text);
	
	cout << " Size of LCE structure (Bytes) " << lce.bit_size()/8 << endl;

return lce;
}

INT ssa_lcp(auto lce, vector<INT> * SSA, vector<INT> * LCP )
{
	cout << "Computing suffix array ... " << flush;

	//SA-RK algorithm
	
	std::sort(SSA->begin(), SSA->end(), lce.lex_less_than());
	
	LCP->push_back(0);
	for(INT i = 1; i<SSA->size(); i++)
	{
    		INT lcp = lce.LCE( SSA[i-1], SSA[i] );
    		LCP->push_back(lcp);
	}
    
return 0;
}
