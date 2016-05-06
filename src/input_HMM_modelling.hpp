//
//  input.hpp
//  
//
//  Created by Ginny Li on 19/10/15.
//
//

#ifndef input_HMM_modelling_hpp
#define input_HMM_modelling_hpp

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <unordered_set>
#include <algorithm>
#include <limits>
#include <tuple>
#include <cmath>
#include <cctype>

//http://www.boost.org/doc/libs/master/libs/algorithm/doc/html/algorithm/Searching.html
// #include "boyer_moore_horspool.hpp"
//#include <boost/algorithm/searching/boyer_moore_horspool.hpp>
//#include <boost/foreach.hpp>
//#include <boost/functional/hash.hpp>



using namespace std;

/*
struct Para_values
{
    int block_length;
    int block_overlaps;
    double stop_criterion;
};
*/

//future version



struct Para_values
{
    int block_length;
    int block_overlaps;
    double stop_criterion;
    
    string fasta_name;
    //string peptide_name;
    string pepsite_name;
    
    string pred_sites;
    
    double prob_threshold;
    
    string hydroemi_name;
    
    string pka1emi_name;
    
    string helixproemi_name;
    
    string stericparemi_name;
    
    string polaremi_name;
    
    string volumeemi_name;
    
    string sheetproemi_name;
    
    string aafreqemi_name;
    
    
    string hydro_choose;
    
    string pka1_choose;
    
    string helixpro_choose;
    
    string stericpar_choose;
    
    string polar_choose;
    
    string volume_choose;
    
    string sheetpro_choose;
    
    string aafreq_choose;
    
    string back_states;

    double t_signi;
};



/*
template <class CharT, class Traits>
void check_init_istream(::basic_istream<CharT, Traits>& is)
{
    is.exceptions(::ios_base::failbit);
    is.exceptions(::ios_base::badbit);
}
*/


vector<string> read_input(const string &fn)
{
    
    ifstream is(fn);
    check_init_istream(is);
    vector<string>myinput;
    
    string lines;
    
    for(;;)
    {
        getline(is,lines);
        
        if(!is) break;
        myinput.emplace_back(lines);
        
        
    }
    
    return myinput;
}


vector<string>get_parameterlines(vector<string>myinput)
{
    vector<string> parastring;
    for (size_t i=0; i<myinput.size(); i++)
    {
        if (myinput[i][0]=='>')
        {
            parastring.push_back(myinput[i]);
        }
        
        
    }
    
    return parastring;
    
}

//make changes later


Para_values  get_paravalues(vector<string>parastring)
{
    Para_values paravalues;
    
    
        for (size_t p=0; p<parastring[0].size(); p++)
        {
            if (parastring[0][p]=='=')
            {
                
                paravalues.block_length=stoi(parastring[0].substr(p+1,(parastring[0].size()-p-1)));
                
                
            }
            
        }
    for (size_t p=0; p<parastring[1].size(); p++)
    {
        if (parastring[1][p]=='=')
        {
            
            paravalues.block_overlaps=stoi(parastring[1].substr(p+1,(parastring[1].size()-p-1)));
            
            
        }
        
    }

    for (size_t p=0; p<parastring[2].size(); p++)
    {
        if (parastring[2][p]=='=')
        {
            
            paravalues.stop_criterion=stod(parastring[2].substr(p+1,(parastring[2].size()-p-1)));
            
            
        }
        
    }

    
    for (size_t p=0; p<parastring[3].size(); p++)
    {
        if (parastring[3][p]=='=')
        {
            
            paravalues.fasta_name=parastring[3].substr(p+1,(parastring[3].size()-p-1));
            
            
        }
        
    }

    
//    for (size_t p=0; p<parastring[4].size(); p++)
//    {
//        if (parastring[4][p]=='=')
//        {
//            
//            paravalues.peptide_name=parastring[4].substr(p+1,(parastring[4].size()-p-1));
//            
//            
//        }
//        
//    }
    
    for (size_t p=0; p<parastring[4].size(); p++)
    {
        if (parastring[4][p]=='=')
        {
            
            
            paravalues.pepsite_name=parastring[4].substr(p+1,(parastring[4].size()-p-1));
            
            
        }
        
    }

    for (size_t p=0; p<parastring[5].size(); p++)
    {
        if (parastring[5][p]=='=')
        {
            //just one site now
            
            paravalues.pred_sites=parastring[5].substr(p+1,(parastring[5].size()-p-1));
            
            
        }
        
    }

    
    
    
    for (size_t p=0; p<parastring[6].size(); p++)
    {
        if (parastring[6][p]=='=')
        {
            
            paravalues.prob_threshold=stod(parastring[6].substr(p+1,(parastring[6].size()-p-1)));
            
            
        }
        
    }
    
    
    
    for (size_t p=0; p<parastring[7].size(); p++)
    {
        if (parastring[7][p]=='=')
        {
           
            paravalues.hydroemi_name=parastring[7].substr(p+1,(parastring[7].size()-p-1));
            
            
        }
        
    }
    
    for (size_t p=0; p<parastring[8].size(); p++)
    {
        if (parastring[8][p]=='=')
        {
            
            paravalues.pka1emi_name=parastring[8].substr(p+1,(parastring[8].size()-p-1));
            
            
        }
        
    }

    for (size_t p=0; p<parastring[9].size(); p++)
    {
        if (parastring[9][p]=='=')
        {
            
            paravalues.helixproemi_name=parastring[9].substr(p+1,(parastring[9].size()-p-1));
            
            
        }
        
    }

    
    for (size_t p=0; p<parastring[10].size(); p++)
    {
        if (parastring[10][p]=='=')
        {
            
            paravalues.stericparemi_name=parastring[10].substr(p+1,(parastring[10].size()-p-1));
            
            
        }
        
    }

    for (size_t p=0; p<parastring[11].size(); p++)
    {
        if (parastring[11][p]=='=')
        {
            
            paravalues.polaremi_name=parastring[11].substr(p+1,(parastring[11].size()-p-1));
            
            
        }
        
    }

    for (size_t p=0; p<parastring[12].size(); p++)
    {
        if (parastring[12][p]=='=')
        {
            
            paravalues.volumeemi_name=parastring[12].substr(p+1,(parastring[12].size()-p-1));
            
            
        }
        
    }

    for (size_t p=0; p<parastring[13].size(); p++)
    {
        if (parastring[13][p]=='=')
        {
            
            paravalues.sheetproemi_name=parastring[13].substr(p+1,(parastring[13].size()-p-1));
            
            
        }
        
    }

    
    for (size_t p=0; p<parastring[14].size(); p++)
    {
        if (parastring[14][p]=='=')
        {
            
            paravalues.aafreqemi_name=parastring[14].substr(p+1,(parastring[14].size()-p-1));
            
            
        }
        
    }

    
    for (size_t p=0; p<parastring[15].size(); p++)
    {
        if (parastring[15][p]=='=')
        {
            
            paravalues.hydro_choose=parastring[15].substr(p+1,(parastring[15].size()-p-1));
            
            
        }
        
    }
    

    for (size_t p=0; p<parastring[16].size(); p++)
    {
        if (parastring[16][p]=='=')
        {
            
            paravalues.pka1_choose=parastring[16].substr(p+1,(parastring[16].size()-p-1));
            
            
        }
        
    }
    
    for (size_t p=0; p<parastring[17].size(); p++)
    {
        if (parastring[17][p]=='=')
        {
            
            paravalues.helixpro_choose=parastring[17].substr(p+1,(parastring[17].size()-p-1));
            
            
        }
        
    }
    
    for (size_t p=0; p<parastring[18].size(); p++)
    {
        if (parastring[18][p]=='=')
        {
            
            paravalues.stericpar_choose=parastring[18].substr(p+1,(parastring[18].size()-p-1));
            
            
        }
        
    }
    
    for (size_t p=0; p<parastring[19].size(); p++)
    {
        if (parastring[19][p]=='=')
        {
            
            paravalues.polar_choose=parastring[19].substr(p+1,(parastring[19].size()-p-1));
            
            
        }
        
    }
    
    for (size_t p=0; p<parastring[20].size(); p++)
    {
        if (parastring[20][p]=='=')
        {
            
            paravalues.volume_choose=parastring[20].substr(p+1,(parastring[20].size()-p-1));
            
            
        }
        
    }
    
    
    for (size_t p=0; p<parastring[21].size(); p++)
    {
        if (parastring[21][p]=='=')
        {
            
            paravalues.sheetpro_choose=parastring[21].substr(p+1,(parastring[21].size()-p-1));
            
            
        }
        
    }
    
    
    for (size_t p=0; p<parastring[22].size(); p++)
    {
        if (parastring[22][p]=='=')
        {
            
            paravalues.aafreq_choose=parastring[22].substr(p+1,(parastring[22].size()-p-1));
            
            
        }
        
    }
    
    

    for (size_t p=0; p<parastring[23].size(); p++)
    {
        if (parastring[23][p]=='=')
        {
            
            paravalues.back_states=parastring[23].substr(p+1,(parastring[23].size()-p-1));
            
            
        }
        
    }
    
    
    
    
    

    
    for (size_t p=0; p<parastring[24].size(); p++)
    {
        if (parastring[24][p]=='=')
        {
            
            paravalues.t_signi=stod(parastring[24].substr(p+1,(parastring[24].size()-p-1)));
            
            
        }
        
    }
    
    
    
    

    
    
    
    
       return paravalues;
}


















#endif
