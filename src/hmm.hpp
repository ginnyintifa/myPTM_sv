//
//  hmm.hpp
//  humansites
//
//  Created by Ginny Li on 15-8-25.
//  Copyright (c) 2015年 Ginny. All rights reserved.
//

#ifndef humansites_hmm_hpp
#define humansites_hmm_hpp

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
//#include <unordered_set>
#include <algorithm>
#include <limits>
//#include <tuple>
#include <cstdio>
#include <cctype>
//http://www.boost.org/doc/libs/master/libs/algorithm/doc/html/algorithm/Searching.html
//#include <boost/algorithm/searching/boyer_moore_horspool.hpp>
//#include <boost/functional/hash.hpp>
#include <cmath>

#include "math.hpp"

using namespace std;

vector<vector<double>> eachstate_20aa_emi{2};//if I don't set the length of it,how do I assign value?

//vector<vector<pair<double,double>>>hydro_emission;
//vector<vector<pair<double,double>>>pka1_emission;
//vector<vector<pair<double,double>>>helixpro_emission;
//vector<vector<pair<double,double>>>stericpar_emission;
//vector<vector<pair<double,double>>>polar_emission;
//vector<vector<pair<double,double>>>volume_emission;
//vector<vector<pair<double,double>>>sheetpro_emission;
//
//vector<vector<vector<double>>>hydro_allstates_emission;
//vector<vector<vector<double>>>pka1_allstates_emission;
//vector<vector<vector<double>>>helixpro_allstates_emission;
//vector<vector<vector<double>>>stericpar_allstates_emission;
//vector<vector<vector<double>>>polar_allstates_emission;
//vector<vector<vector<double>>>volume_allstates_emission;
//vector<vector<vector<double>>>sheetpro_allstates_emission;
//

double piprobs[2];

double transmatrix[2][2];

double miu[2];

double sigma2[2];

double tiny;


vector<double> epi(2);

vector<vector<double>> etransmatrix(2,vector<double>(2));

vector<double> emiu(2);



vector<double> esigma2(2);

vector<vector<double>> sumtheta(2,vector<double>(2));

int c=-1;

vector<double> sumgamma(2);




double loglikelihood;

//modi: now I will first get rid of the seq variables.

vector<vector<double>> init_forward(vector<vector<int>>all_aas,const size_t lenn,size_t seqth, /*vector<double>hydro_seq,vector<double>pka1_seq,vector<double>helixpro_seq,*/ const vector<vector<vector<double>>>& hydro_emission,const vector<vector<vector<double>>>& pka1_emission,const vector<vector<vector<double>>>& helixpro_emission,const vector<vector<vector<double>>>& stericpar_emission,const vector<vector<vector<double>>>& polar_emission,vector<vector<vector<double>>>& volume_emission,const vector<vector<vector<double>>>& sheetpro_emission)
{
    vector<vector<double>>forward_my(lenn,vector<double>(2));
    for (int i=0; i<2; i++)
    {
        forward_my[0][i]=log(piprobs[i])+get_logmultinomial(all_aas[0], eachstate_20aa_emi[i])+hydro_emission[seqth][0][i]+pka1_emission[seqth][0][i]+helixpro_emission[seqth][0][i]+stericpar_emission[seqth][0][i]+polar_emission[seqth][0][i]+volume_emission[seqth][0][i]+sheetpro_emission[seqth][0][i];
;
    }
    
    
    for (int t=1; t<lenn; t++)
    {
        for (int i=0; i<2; i++)
        {
            double a=forward_my[t-1][0]+log(transmatrix[0][i]);
            double b=forward_my[t-1][1]+log(transmatrix[1][i]);
            //double c=forward_aa[t-1][2]+log(transmatrix[2][i]);
            //double d=forward_aa[t-1][3]+log(transmatrix[3][i]);
            
            double pro=logsum(a, b);
            
          
            forward_my[t][i]=pro+get_logmultinomial(all_aas[t], eachstate_20aa_emi[i])+hydro_emission[seqth][t][i]+pka1_emission[seqth][t][i]+helixpro_emission[seqth][t][i]+stericpar_emission[seqth][t][i]+polar_emission[seqth][t][i]+volume_emission[seqth][t][i]+sheetpro_emission[seqth][t][i];
            
        
        }
        
    }
    
    
    
    return forward_my;
    
}





////////////////////////////////////////today forward finished////////////////////////////////////

////////////////////////////////////////tuesday backward and beyond///////////////////////////////

vector<vector<double>>init_backward(vector<vector<int>>all_aas, const size_t lenn,size_t seqth,/* vector<double>hydro_seq, vector<double>pka1_seq,vector<double>helixpro_seq,*/ const vector<vector<vector<double>>>&hydro_emission,const vector<vector<vector<double>>>&pka1_emission,const vector<vector<vector<double>>>&helixpro_emission,const vector<vector<vector<double>>>& stericpar_emission,const vector<vector<vector<double>>>& polar_emission,vector<vector<vector<double>>>& volume_emission,const vector<vector<vector<double>>>& sheetpro_emission)
{
    vector<vector<double>>backward_my(lenn,vector<double>(2));
    
    for(int i=0; i<2;i++)
    {

        backward_my[lenn-1][i]=0;
    }
    
    for(int t=(lenn-2);t>=0;t--)
    {
        
            for(int i=0;i<2;i++)
            {
                double aa =log(transmatrix[i][0])+get_logmultinomial(all_aas[t+1], eachstate_20aa_emi[0])+hydro_emission[seqth][t+1][0]+pka1_emission[seqth][t+1][0]+helixpro_emission[seqth][t+1][0]+stericpar_emission[seqth][t+1][0]+polar_emission[seqth][t+1][0]+volume_emission[seqth][t+1][0]+sheetpro_emission[seqth][t+1][0]+backward_my[t+1][0];
                
                double bb =log(transmatrix[i][1])+get_logmultinomial(all_aas[t+1], eachstate_20aa_emi[1])+hydro_emission[seqth][t+1][1]+pka1_emission[seqth][t+1][1]+helixpro_emission[seqth][t+1][1]+stericpar_emission[seqth][t+1][1]+polar_emission[seqth][t+1][1]+volume_emission[seqth][t+1][1]+sheetpro_emission[seqth][t+1][1]+backward_my[t+1][1];
              
                double pro=logsum(aa,bb);
                
                backward_my[t][i]=pro;
                
              
            }
        
      
    }
    
    
    return backward_my;
    
}


vector<vector<double>>assign_for(vector<vector<double>>forward, const size_t lenn)
{
    vector<vector<double>> eforward(lenn,vector<double>(2));
    
    for (int t=0;t<lenn;t++)
    {
        for (int i=0;i<2;i++)
        {
            eforward[t][i]=forward[t][i];
        }
    }
    return eforward;
}





vector<vector<double>>assign_back(vector<vector<double>>backward, const size_t lenn)
{
    vector<vector<double>> ebackward(lenn,vector<double>(2));
    
    for (int t=0;t<lenn;t++)
    {
        for (int i=0;i<2;i++)
        {
            ebackward[t][i]=backward[t][i];
        }
    }
    
    return ebackward;
}





vector<vector<double>> get_gamma(const size_t lenn)
{
    vector<vector<double>>ggamma(lenn,vector<double>(2));
    return ggamma;
}



vector<vector<vector<double>>> get_theta(const size_t lenn)
{
    vector<vector<vector<double>>>theta(lenn,vector<vector<double>>(2,vector<double>(2)));
    return theta;
}


vector<int> get_state(const size_t lenn)
{
    vector<int> histate(lenn);
    return histate;
}


//what does obs_ do?
//modi:set the obs_hydro etc ignored now

vector<vector<vector<double>>> bw_theta(size_t seqth,vector<vector<int>>all_aas,/*vector<double>obs_hydro,vector<double>obs_pka1,vector<double>obs_helixpro,*/const vector<vector<vector<double>>>&hydro_emission,const vector<vector<vector<double>>>&pka1_emission,const vector<vector<vector<double>>>&helixpro_emission,const vector<vector<vector<double>>>& stericpar_emission,const vector<vector<vector<double>>>& polar_emission,vector<vector<vector<double>>>& volume_emission,const vector<vector<vector<double>>>& sheetpro_emission,size_t lenn, vector<vector<double>> &eforward,vector<vector<double>>& ebackward,vector<vector<vector<double>>>& theta)

{
    
    
    for (int t=0;t<lenn-1;t++)
    {
        
        for(int i=0;i<2;i++)
        {
            for (int j=0;j<2;j++)
            {
            
                theta[t][i][j]=eforward[t][i]+etransmatrix[i][j]+get_logmultinomial(all_aas[t+1], eachstate_20aa_emi[j])+hydro_emission[seqth][t+1][j]+pka1_emission[seqth][t+1][j]+helixpro_emission[seqth][t+1][j]+stericpar_emission[seqth][t+1][j]+polar_emission[seqth][t+1][j]+volume_emission[seqth][t+1][j]+sheetpro_emission[seqth][t+1][j]+ebackward[t+1][j];
                
            }
            
            
        }
    
        //to see if the inputs are necessary
        
        //cout<<polar_emission[seqth][t+1][1];
        
        
        double thetaarray[4]={theta[t][0][0],theta[t][0][1],theta[t][1][0],theta[t][1][1]};
        
        
        double dvnormtheta=logsum4(thetaarray);
        
        for (int i=0;i<2; i++)
        {
            for (int j=0;j<2;j++)
            {
                theta[t][i][j]=theta[t][i][j]-dvnormtheta;
                
                
            }
            
        }
        
    }
    
    return theta;
}









vector<vector<double>> bw_gamma (size_t lenn,vector<vector<double>> &eforward,vector<vector<double>>& ggamma,vector<vector<vector<double>>>& theta)

{
    for(int t=0;t<lenn-1;t++)
    {
        for (int i=0;i<2;i++)
        {
            ggamma[t][i]=logsum(theta[t][i][0], theta[t][i][1]);
            
        }
        
        double dnormgamma=logsum(ggamma[t][0],ggamma[t][1]);
        
        for (int i=0;i<2;i++)
        {
            ggamma[t][i]=ggamma[t][i]-dnormgamma;
        }
        
    }
    
    
    
    //get gamma last one
    
    
    for (int i=0;i<2;i++)
    {
        ggamma[lenn-1][i]=eforward[lenn-1][i];
    }
    // normalize ggamma;
    
    double dnormgamma=logsum(ggamma[lenn-1][0],ggamma[lenn-1][1]);

    
    for (int i=0;i<2;i++)
    {
        ggamma[lenn-1][i]=ggamma[lenn-1][i]-dnormgamma;
    }
    
    
    
    
    
    return ggamma;
}



//marginal

//ok I want to try this insert a threshold to it

vector<int>get_histates(size_t lenn,vector<vector<double>>& ggamma, vector<int>& histate,double mythreshold)
{
    for(int t=0; t<lenn;t++)
    {
        
        
        if (exp(ggamma[t][1])>mythreshold)
        {
            histate[t]=1;
        }
        else
            histate[t]=0;
        
       // histate[t]=maxarg(ggamma[t][0],ggamma[t][1]);
        
        
    }
    
    return histate;
    
}



vector<vector<double>>sum_theta(size_t lenn,vector<vector<vector<double>>>& theta)

{
    for (int i=0;i<2;i++)
    {
        for (int j=0;j<2;j++)
        {
            
            
            vector<double> thetaarray(lenn-1);
            
            for(int t=0;t<lenn-1;t++)
            {
                thetaarray[t]=theta[t][i][j];
            }
            
            
            sumtheta[i][j]=logsuma(thetaarray);
            
        }
    }
    
    return sumtheta;
}


vector<double>sum_gamma(size_t lenn, vector<vector<double>>& ggamma)
{
    for (int i=0;i<2;i++)
    {
        
        
        
        vector<double>ggammaarray(lenn-1);
        
        for(int t=0;t<lenn-1;t++)
        {
            ggammaarray[t]=ggamma[t][i];
        }
        sumgamma[i]=logsuma(ggammaarray);
        
    }
    return sumgamma;
}





vector<double> update_epi(vector<vector<vector<double>>>& all_gamma)
{
    
    
    for (int i=0;i<2;i++)
    {
        
        
        vector<double>array_gamma0;
        for (size_t p=0;p<all_gamma.size();p++)
        {
            
            array_gamma0.push_back(all_gamma[p][0][i]);
            
        }
        
        
        epi[i]=logsuma(array_gamma0);
        epi[i]=epi[i]-log(all_gamma.size());
        
        if (epi[i]<-230)
        {
            epi[i]=-230;//log is e based
        }
        
        
    }
    
    double epinorm=logsum(epi[0],epi[1]);
    
    for (int u=0; u<2; u++)
    {
        epi[u]-=epinorm;
    }
    
    
    
    return epi;
}





vector<vector<double>>update_etransmatirx(vector<vector<vector<double>>>& all_sumtheta, vector<vector<double>>& all_sumgamma)
{
    for (int i=0;i<2;i++)
    {
        
        for (int j=0;j<2;j++)
        {
            vector<double>array_trans(all_sumgamma.size());
            
            for(size_t p=0;p<all_sumtheta.size();p++)
            {
                array_trans[p]=all_sumtheta[p][i][j]-all_sumgamma[p][i];
                
            }
            
            etransmatrix[i][j]=logsuma(array_trans)-log(all_sumgamma.size());//to get the average value;
            
            if (etransmatrix[i][j]<-230)
            {
                etransmatrix[i][j]=-230;
            }
            
        }
        
        
        double tonorm=logsum(etransmatrix[i][0],etransmatrix[i][1]);
        
        
        for (int j=0;j<2;j++)
        {
            etransmatrix[i][j]-=tonorm;
        }
    
        cout<<"trans "<<exp(etransmatrix[i][0])<<" "<<exp(etransmatrix[i][1])<<endl;
        
        
        cout<<"sum "<<exp(etransmatrix[i][0])+exp(etransmatrix[i][1])<<endl;
    }
    
    
    return etransmatrix;
}


//since I use prefixed emiu and esigma2 "update_emiu" and "update_esigma2" are not used.

vector<double> update_emiu(vector<vector<double>> myobs, vector<size_t>all_lenn, vector<vector<vector<double>>>& all_gamma)
{
    for(int i=0;i<3;i++)
    {
        double sumgo=0;
        double sumga=0;
        
        for(size_t p=0;p<myobs.size();p++)
        {
            for(size_t t=0;t<all_lenn[p];t++)
            {
               sumgo=sumgo+myobs[p][t]*exp(all_gamma[p][t][i]);
               sumga=sumga+exp(all_gamma[p][t][i]);
            }
        }
        
        emiu[i]=sumgo/(sumga+tiny);
    }
    
    return emiu;
}



vector<double> update_esigma2(vector<vector<double>> myobs,vector<size_t> all_lenn, vector<vector<vector<double>>>& all_gamma)
{
    
    
    for (int i=0;i<3;i++)
    {
        
        double sumga=0;
        double sumsig=0;
        for(size_t p=0; p<myobs.size();p++)
        {
            for (int t=0;t<all_lenn[p];t++)
            {
            sumsig=sumsig+exp(all_gamma[p][t][i])*(myobs[p][t]-emiu[i])*(myobs[p][t]-emiu[i]);
            sumga=sumga+exp(all_gamma[p][t][i]);
            }
        }
       
        esigma2[i]=sumsig/(sumga+tiny);
        
    }
    
    return esigma2;
}


//modi:ignore seq

//note: actually the _emission are not necessary

//note: next stage get rid of them

vector<vector<double>> update_eforward(vector<vector<int>>all_aas,size_t seqth, /*vector<double>hydro_seq,vector<double>pka1_seq,vector<double>helixpro_seq,*/ const size_t lenn,const vector<vector<vector<double>>>& hydro_emission,const vector<vector<vector<double>>>& pka1_emission,vector<vector<vector<double>>>& helixpro_emission,const vector<vector<vector<double>>>& stericpar_emission,const vector<vector<vector<double>>>& polar_emission,vector<vector<vector<double>>>& volume_emission,const vector<vector<vector<double>>>& sheetpro_emission,vector<vector<double>>&eforward)
{
    
    for (int i=0; i<2; i++)
    {
        eforward[0][i]=epi[i]+get_logmultinomial(all_aas[0], eachstate_20aa_emi[i])+hydro_emission[seqth][0][i]+pka1_emission[seqth][0][i]+helixpro_emission[seqth][0][i]+stericpar_emission[seqth][0][i]+polar_emission[seqth][0][i]+volume_emission[seqth][0][i]+sheetpro_emission[seqth][0][i];
    }
    
    
    for (int t=1; t<lenn; t++)
    {
        for (int i=0; i<2; i++)
        {
            double a=eforward[t-1][0]+log(transmatrix[0][i]);
            double b=eforward[t-1][1]+log(transmatrix[1][i]);
            //double c=forward_aa[t-1][2]+log(transmatrix[2][i]);
            //double d=forward_aa[t-1][3]+log(transmatrix[3][i]);
            
            double pro=logsum(a, b);
            
            
            eforward[t][i]=pro+get_logmultinomial(all_aas[t], eachstate_20aa_emi[i])+hydro_emission[seqth][t][i]+pka1_emission[seqth][t][i]+helixpro_emission[seqth][t][i]+stericpar_emission[seqth][t][i]+polar_emission[seqth][t][i]+volume_emission[seqth][t][i]+sheetpro_emission[seqth][t][i];
            
            
        }
        
    }
    
    
    
    return eforward;
    
}




vector<vector<double>>update_ebackward(vector<vector<int>>all_aas, size_t seqth,/* vector<double>hydro_seq, vector<double>pka1_seq, vector<double>helixpro_seq, */const size_t lenn,const vector<vector<vector<double>>>&hydro_emission,const vector<vector<vector<double>>>&pka1_emission, const vector<vector<vector<double>>>& helixpro_emission,const vector<vector<vector<double>>>& stericpar_emission,const vector<vector<vector<double>>>& polar_emission,vector<vector<vector<double>>>& volume_emission,const vector<vector<vector<double>>>& sheetpro_emission,vector<vector<double>>&ebackward)
{
    
    
    for(int i=0; i<2;i++)
    {
        
        ebackward[lenn-1][i]=0;
    }
    
    for(int t=(lenn-2);t>=0;t--)
    {
        
        for(int i=0;i<2;i++)
        {
            double aa =log(transmatrix[i][0])+get_logmultinomial(all_aas[t+1], eachstate_20aa_emi[0])+hydro_emission[seqth][t+1][0]+pka1_emission[seqth][t+1][0]+helixpro_emission[seqth][t+1][0]+stericpar_emission[seqth][t+1][0]+polar_emission[seqth][t+1][0]+volume_emission[seqth][t+1][0]+sheetpro_emission[seqth][t+1][0]+ebackward[t+1][0];
            double bb =log(transmatrix[i][1])+get_logmultinomial(all_aas[t+1], eachstate_20aa_emi[1])+hydro_emission[seqth][t+1][1]+pka1_emission[seqth][t+1][1]+helixpro_emission[seqth][t+1][1]+stericpar_emission[seqth][t+1][1]+polar_emission[seqth][t+1][1]+volume_emission[seqth][t+1][1]+sheetpro_emission[seqth][t+1][1]+ebackward[t+1][1];
            
            double pro=logsum(aa,bb);
            
            ebackward[t][i]=pro;
            
            
        }
        
        
    }
    
    
    return ebackward;
    
}





double get_loglikelihood(size_t lenn, const vector<vector<double>>&eforward)
{
    vector<double> eforaarray;
    
    for (size_t u=0; u<2; u++)
    {
        eforaarray.push_back(eforward[lenn-1][u]);
        
        
    }

    
    return logsuma(eforaarray);
    
}


//viterbi algorithm

//best score along a single path, at time t, which accounts for the first t observations and ends in state Si

//we need to keep track of the arguments which maximized

//now I will make use of the parameters got from the bw algorithm


// viterbi is done for every sequence.


//vector<int> viterbi_foreach(vector<vector<double>> etransmatrix,vector<double>epi,vector<vector<int>> all_aas, const size_t lenn)




#endif
