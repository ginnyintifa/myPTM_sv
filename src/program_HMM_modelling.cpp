//
//  main.cpp
//  humansites
//
//  Created by Ginny Li on 15-7-28.
//  Copyright (c) 2015年 Ginny. All rights reserved.
//

//a modification
//in R want to label ones that are filtered to be  green,



// a modification of my program

//change to only hydro and pka1



//I mistakenly made some wrong changes


#include "site.hpp"
#include "aaprop.hpp"
#include "emisson.hpp"
#include "hmm.hpp"
#include "math.hpp"
#include "input_HMM_modelling.hpp"
#include "ttest.hpp"
#include "chitest.hpp"


//now, next stage of modification, to make the emission in input

//think about how to make it as the vector of input and make this part



int main(int argc, char *argv[])

{
    
    using namespace std;
    
    //there is one modification on mypeptides_site, change it to myfirst_sites for convinience
    
    vector<string>myinput=read_input(argv[1]);
    vector<string>mypara=get_parameterlines(myinput);
    
    Para_values myparavalue=get_paravalues(mypara);
    
    
    
    
    
    vector<int> myaawhich=get_aawhich(myparavalue.pred_sites);
    
    
    vector<prot_type> myproteins=read_fasta(myparavalue.fasta_name);//"uniprot-SGD.fasta");
    vector<string> mypeptides_site=read_tsv(myparavalue.pepsite_name);//"yeastsitepep.tsv");
    
    //vector<string> myfirst_sites;
    
    vector<aa_properties> aminoacid_properties(24);
    
    vector<double>each_hydro;
    vector<double>each_pka1;
    vector<double>each_helixpro;
    vector<double>each_stericpar;
    vector<double>each_polar;
    vector<double>each_volume;
    vector<double>each_sheetpro;
    
    assign_properties(aminoacid_properties);

    assign_each_property(aminoacid_properties,each_hydro,each_pka1,each_helixpro,each_stericpar,each_polar,each_volume,each_sheetpro);
    
    
    
    
    for (size_t i=0; i<mypeptides_site.size(); i++)
    {
        for (size_t p=0; p<mypeptides_site[i].size(); p++)
        {
            int thesign=0;
            
            for (size_t q=0; q<myparavalue.pred_sites.size(); q++)
            {
                if (mypeptides_site[i][p]==myparavalue.pred_sites[q])
                {
                    thesign=1;
                }
                
                
            }
            
            if (thesign==0)
            {
                mypeptides_site[i][p]=toupper(mypeptides_site[i][p]);
            }
        }
    }
    
    
    
    
    
    vector<string>mypeptides=mypeptides_site;
    
    for (size_t i=0; i<mypeptides_site.size(); i++)
    {
        for (size_t p=0; p<mypeptides_site[i].size(); p++)
        {
            mypeptides[i][p]=toupper(mypeptides[i][p]);
        }
    }
    
    
    myproteins=read_concise_pair(myproteins);
    
    mypeptides=read_concise(mypeptides);
    
    
    vector<Pep_in_prot> mylinks=link_peptides_with_sites(mypeptides, mypeptides_site,myparavalue.pred_sites);
    
    pair<vector<prot_type>,vector<vector<Pep_in_prot>>>ppmap=get_mapps(myproteins, mylinks);
    vector<pair<prot_type,vector<site_type>>> myprotsitesall=get_all_seq_sites(ppmap);
    vector<pair<prot_type,vector<site_type>>> myneededprotsiteall=get_needed_myprot_sites_all(myprotsitesall);
    vector<pair<prot_type,vector<site_type>>> myuniquesites=get_unique_site(myneededprotsiteall);
    
    sort(myuniquesites.begin(), myuniquesites.end(), sortpair);
    
    vector<vector<string>> window_for_all=get_windows_for_all(myuniquesites,myparavalue.block_length,myparavalue.block_overlaps);
    
    vector<vector<vector<int>>> aa_compo_forall=get_allseqs__aavector(window_for_all, aminoacid_properties);
    
    vector<vector<int>>aa_vectorforall;
    
    for (size_t i=0; i<myuniquesites.size(); i++)
    {
        vector<int> aa_vectorforone;
        
        aa_vectorforone=grab_aa_vector(aminoacid_properties, myuniquesites[i].first.second);
        
        aa_vectorforall.push_back(aa_vectorforone);
    }
    
    
    ofstream outfile;
    
       
    
    
    
    
    vector<double>seqs_means_hydro=get_means(each_hydro,aa_vectorforall);
    
    vector<double>seqs_means_pka1=get_means(each_pka1,aa_vectorforall);
    
    vector<double>seqs_means_helixpro=get_means(each_helixpro,aa_vectorforall);
    
    vector<double>seqs_means_stericpar=get_means(each_stericpar,aa_vectorforall);
    
    vector<double>seqs_means_polar=get_means(each_polar,aa_vectorforall);
    
    vector<double>seqs_means_volume=get_means(each_volume,aa_vectorforall);
    
    vector<double>seqs_means_sheetpro=get_means(each_sheetpro,aa_vectorforall);
    

    
    
    vector<pair<prot_type,vector<vector<site_type>>>>window_uniquesites_forall=get_unique_site_inblock(myuniquesites,myparavalue.block_length,myparavalue.block_overlaps);
    
    
    //bug happened here
    
    vector<vector<string>>marked_window_forall=marked_windows(window_uniquesites_forall, window_for_all,myparavalue.block_length,myparavalue.block_overlaps);
    vector<vector<int>>windows_states=get_states(window_uniquesites_forall,myparavalue.pred_sites);
    
    cout<<"get states...\n";
    
    
    
    
    vector<vector<double>> hydro_forall=get_allseqs_property(aa_compo_forall,each_hydro,seqs_means_hydro);//normalized
    
    vector<vector<double>> pka1_forall=get_allseqs_property(aa_compo_forall,each_pka1,seqs_means_pka1);
    
    vector<vector<double>> helixpro_forall=get_allseqs_property(aa_compo_forall,each_helixpro,seqs_means_helixpro);
    vector<vector<double>> stericpar_forall=get_allseqs_property(aa_compo_forall,each_stericpar,seqs_means_hydro);//normalized
    vector<vector<double>> polar_forall=get_allseqs_property(aa_compo_forall,each_polar,seqs_means_pka1);
    vector<vector<double>> volume_forall=get_allseqs_property(aa_compo_forall,each_volume,seqs_means_helixpro);
    vector<vector<double>> sheetpro_forall=get_allseqs_property(aa_compo_forall,each_sheetpro,seqs_means_helixpro);
    

    
    
    
    
    
        vector<vector<pair<int, int>>>states_windowrates=get_statesrates(aa_compo_forall, windows_states);
    vector<pair<prot_type,vector<int>>>allseqs_blockwithsites=get_blocks_withsites_all(window_uniquesites_forall);
    
    //this is what I want I think
    
  
    
    
    
    
    
    
    vector<vector<double>>hydro_2states=get_properties_2states(hydro_forall,windows_states);
    
    vector<vector<double>>pka1_2states=get_properties_2states(pka1_forall,windows_states);
    
    vector<vector<double>>helixpro_2states=get_properties_2states(helixpro_forall,windows_states);

    vector<vector<double>>stericpar_2states=get_properties_2states(stericpar_forall,windows_states);
    
    vector<vector<double>>polar_2states=get_properties_2states(polar_forall,windows_states);
    //Xchange
    
    vector<vector<double>>volume_2states=get_properties_2states(volume_forall,windows_states);
    vector<vector<double>>sheetpro_2states=get_properties_2states(sheetpro_forall,windows_states);
    

    
    
    cout<<"@@@@@@@@@@@@@@@@@ "<<hydro_2states.size()<<endl;
   
      outfile.open("hydro_dens_0.tsv");
    
      for (size_t i=0; i<hydro_2states[0].size(); i++)
        {
            outfile<<hydro_2states[0][i]<<endl;
        }
  
    
    outfile.close();
    
    outfile.open("hydro_dens_1.tsv");
    
    for (size_t i=0; i<hydro_2states[1].size(); i++)
    {
        outfile<<hydro_2states[1][i]<<endl;
    }
    
    
    outfile.close();
    
    
    outfile.open("pka1_dens_0.tsv");
    
            for (size_t i=0; i<pka1_2states[0].size(); i++)
        {
            outfile<<pka1_2states[0][i]<<endl;
        }
   
    
    outfile.close();
    
    outfile.open("pka1_dens_1.tsv");
    
    for (size_t i=0; i<pka1_2states[1].size(); i++)
    {
        outfile<<pka1_2states[1][i]<<endl;
    }
    
    
    outfile.close();
    
    
    outfile.open("helixpro_dens_0.tsv");
    
    for (size_t i=0; i<helixpro_2states[0].size(); i++)
    {
        outfile<<helixpro_2states[0][i]<<endl;
    }
    
    
    outfile.close();
    
    outfile.open("helixpro_dens_1.tsv");
    
    for (size_t i=0; i<helixpro_2states[1].size(); i++)
    {
        outfile<<helixpro_2states[1][i]<<endl;
    }
    
    
    outfile.close();
    
    
    outfile.open("stericpar_dens_0.tsv");
    
    for (size_t i=0; i<stericpar_2states[0].size(); i++)
    {
        outfile<<stericpar_2states[0][i]<<endl;
    }
    
    
    outfile.close();
    
    outfile.open("stericpar_dens_1.tsv");
    
    for (size_t i=0; i<stericpar_2states[1].size(); i++)
    {
        outfile<<stericpar_2states[1][i]<<endl;
    }
    
    
    outfile.close();
    

    
    outfile.open("polar_dens_0.tsv");
    
    for (size_t i=0; i<polar_2states[0].size(); i++)
    {
        outfile<<polar_2states[0][i]<<endl;
    }
    
    
    outfile.close();
    
    outfile.open("polar_dens_1.tsv");
    
    for (size_t i=0; i<polar_2states[1].size(); i++)
    {
        outfile<<polar_2states[1][i]<<endl;
    }
    
    
    outfile.close();
    
    outfile.open("volume_dens_0.tsv");
    
    for (size_t i=0; i<volume_2states[0].size(); i++)
    {
        outfile<<volume_2states[0][i]<<endl;
    }
    
    
    outfile.close();
    
    outfile.open("volume_dens_1.tsv");
    
    for (size_t i=0; i<volume_2states[1].size(); i++)
    {
        outfile<<volume_2states[1][i]<<endl;
    }
    
    
    outfile.close();
    

    
    outfile.open("sheetpro_dens_0.tsv");
    
    for (size_t i=0; i<sheetpro_2states[0].size(); i++)
    {
        outfile<<sheetpro_2states[0][i]<<endl;
    }
    
    
    outfile.close();
    
    outfile.open("sheetpro_dens_1.tsv");
    
    for (size_t i=0; i<sheetpro_2states[1].size(); i++)
    {
        outfile<<sheetpro_2states[1][i]<<endl;
    }
    
    
    outfile.close();
    
  
    
    




    
    cout<<"get my peptides properties...\n";
    
    
    
    //from this onward, not only do with aacompoforall
    
    
    
    vector<vector<pair<double,double>>>hydro_emission;
    vector<vector<pair<double,double>>>pka1_emission;
    vector<vector<pair<double,double>>>helixpro_emission;
    vector<vector<pair<double,double>>>stericpar_emission;
    vector<vector<pair<double,double>>>polar_emission;
    vector<vector<pair<double,double>>>volume_emission;
    vector<vector<pair<double,double>>>sheetpro_emission;
    
      if (myparavalue.hydroemi_name.size()>1)//see if user wants to use back infomation
    {
        
        
        
        cout<<myparavalue.hydroemi_name<<endl;
        cout<<myparavalue.pka1emi_name<<endl;
    //Xchange
        cout<<myparavalue.helixproemi_name<<endl;
        cout<<myparavalue.aafreqemi_name<<endl;
        
        
        vector<double>hydrotsv=read_tsv_double(myparavalue.hydroemi_name);
        vector<double>pka1tsv=read_tsv_double(myparavalue.pka1emi_name);
        vector<double>helixprotsv=read_tsv_double(myparavalue.helixproemi_name);
        vector<double>stericpartsv=read_tsv_double(myparavalue.stericparemi_name);
        vector<double>polartsv=read_tsv_double(myparavalue.polaremi_name);
        vector<double>volumetsv=read_tsv_double(myparavalue.volumeemi_name);
        vector<double>sheetprotsv=read_tsv_double(myparavalue.sheetproemi_name);

        vector<double>aafreqtsv=read_tsv_double(myparavalue.aafreqemi_name);
        
        
        cerr<<"finish reading...\n";
        
        hydro_emission=get_emi_structure(hydrotsv);
        pka1_emission=get_emi_structure(pka1tsv);
       helixpro_emission=get_emi_structure(helixprotsv);
        stericpar_emission=get_emi_structure(stericpartsv);
    polar_emission=get_emi_structure(polartsv);
    volume_emission=get_emi_structure(volumetsv);
    sheetpro_emission=get_emi_structure(sheetprotsv);


        
        eachstate_20aa_emi=get_aafreq_structure(aafreqtsv);
        cerr<<"get the emission ready...\n";
        
        
        
        cout<<"eachstate_20aa_emi.. "<<eachstate_20aa_emi.size()<<" "<<eachstate_20aa_emi[0].size()<<" "<<eachstate_20aa_emi[1].size()<<endl;

        
    }
    
    
    
    else
    {
    assign_aaemissions(states_windowrates,eachstate_20aa_emi);
    
    
    hydro_2states=choose_data(hydro_2states);
    pka1_2states=choose_data(pka1_2states);
    helixpro_2states=choose_data(helixpro_2states);
        stericpar_2states=choose_data(stericpar_2states);
        polar_2states=choose_data(polar_2states);
        volume_2states=choose_data(volume_2states);
        sheetpro_2states=choose_data(sheetpro_2states);


    
   hydro_emission=get_emission(hydro_2states);
    pka1_emission=get_emission(pka1_2states);
    helixpro_emission=get_emission(helixpro_2states);
        stericpar_emission=get_emission(stericpar_2states);
        polar_emission=get_emission(polar_2states);
        volume_emission=get_emission(volume_2states);
   sheetpro_emission=get_emission(sheetpro_2states);

    cerr<<"prob calculated...\n";
        
        //change the interface
    }
    
    
    
    
    
    
    vector<vector<vector<double>>>hydro_allstates_emission=get_allstates_emission(hydro_forall, hydro_emission);
    vector<vector<vector<double>>>pka1_allstates_emission=get_allstates_emission(pka1_forall, pka1_emission);
    vector<vector<vector<double>>>helixpro_allstates_emission=get_allstates_emission(helixpro_forall, helixpro_emission);
    vector<vector<vector<double>>>stericpar_allstates_emission=get_allstates_emission(stericpar_forall, stericpar_emission);
    vector<vector<vector<double>>>polar_allstates_emission=get_allstates_emission(polar_forall, polar_emission);
    vector<vector<vector<double>>>volume_allstates_emission=get_allstates_emission(volume_forall, volume_emission);
    vector<vector<vector<double>>>sheetpro_allstates_emission=get_allstates_emission(sheetpro_forall, sheetpro_emission);

    
    cerr<<"emission got...\n";

    
    
    
    cout<<"eachstate_20aa_emi.. "<<eachstate_20aa_emi.size()<<" "<<eachstate_20aa_emi[0].size()<<" "<<eachstate_20aa_emi[1].size()<<endl;
    
    
    outfile.open("aafreq_output.tsv");
    
    outfile<<"AA"<<"\t"<<"State_0 freq"<<"\t"<<"State_1 freq"<<"\t"<<"Ratio"<<endl;
    
    for (size_t i=0; i<20; i++)
    {
        cout<<i<<endl;
        outfile<<aminoacid_properties[i].aaname<<"\t"<<eachstate_20aa_emi[0][i]<<"\t"<<eachstate_20aa_emi[1][i]<<"\t"<<double(eachstate_20aa_emi[1][i])/double(eachstate_20aa_emi[0][i])<<endl;
    }
    
    outfile.close();
    
    
    
    
    
    
    
    if (myparavalue.hydro_choose=="NO")
    {
        for (size_t i=0; i<hydro_allstates_emission.size(); i++)
        {
            for (size_t p=0; p<hydro_allstates_emission[i].size(); p++)
            {
                for (size_t j=0; j<hydro_allstates_emission[i][p].size(); j++)
                {
                    hydro_allstates_emission[i][p][j]=0;
                }
            }
        }
        
    }
    
    if (myparavalue.pka1_choose=="NO")
    {
        for (size_t i=0; i<pka1_allstates_emission.size(); i++)
        {
            for (size_t p=0; p<pka1_allstates_emission[i].size(); p++)
            {
                for (size_t j=0; j<pka1_allstates_emission[i][p].size(); j++)
                {
                    pka1_allstates_emission[i][p][j]=0;
                }
            }
        }
        
    }
    
    if (myparavalue.helixpro_choose=="NO")
    {
        for (size_t i=0; i<helixpro_allstates_emission.size(); i++)
        {
            for (size_t p=0; p<helixpro_allstates_emission[i].size(); p++)
            {
                for (size_t j=0; j<helixpro_allstates_emission[i][p].size(); j++)
                {
                    helixpro_allstates_emission[i][p][j]=0;
                }
            }
        }
        
    }
    
    if (myparavalue.stericpar_choose=="NO")
    {
        for (size_t i=0; i<stericpar_allstates_emission.size(); i++)
        {
            for (size_t p=0; p<stericpar_allstates_emission[i].size(); p++)
            {
                for (size_t j=0; j<stericpar_allstates_emission[i][p].size(); j++)
                {
                    stericpar_allstates_emission[i][p][j]=0;
                }
            }
        }
        
    }
    
    
    
    if (myparavalue.polar_choose=="NO")
    {
        for (size_t i=0; i<polar_allstates_emission.size(); i++)
        {
            for (size_t p=0; p<polar_allstates_emission[i].size(); p++)
            {
                for (size_t j=0; j<polar_allstates_emission[i][p].size(); j++)
                {
                    polar_allstates_emission[i][p][j]=0;
                }
            }
        }
        
    }
    

    if (myparavalue.volume_choose=="NO")
    {
        for (size_t i=0; i<volume_allstates_emission.size(); i++)
        {
            for (size_t p=0; p<volume_allstates_emission[i].size(); p++)
            {
                for (size_t j=0; j<volume_allstates_emission[i][p].size(); j++)
                {
                    volume_allstates_emission[i][p][j]=0;
                }
            }
        }
        
    }
    

    if (myparavalue.sheetpro_choose=="NO")
    {
        for (size_t i=0; i<sheetpro_allstates_emission.size(); i++)
        {
            for (size_t p=0; p<sheetpro_allstates_emission[i].size(); p++)
            {
                for (size_t j=0; j<sheetpro_allstates_emission[i][p].size(); j++)
                {
                    sheetpro_allstates_emission[i][p][j]=0;
                }
            }
        }
        
    }
    
    
    
    
    piprobs[0]=0.5;
    piprobs[1]=0.5;
    
    
    transmatrix[0][0]=0.5;
    transmatrix[0][1]=0.5;
    
    transmatrix[1][0]=0.5;
    transmatrix[1][1]=0.5;
    
    epi[0]=0.5;
    epi[1]=0.5;
    
    etransmatrix[0][0]=log(0.5);
    etransmatrix[0][1]=log(0.5);
    
    etransmatrix[1][0]=log(0.5);
    etransmatrix[1][1]=log(0.5);
    
    
   // int K=2;
    
    sumtheta[0][0]=0.0;
    sumtheta[0][1]=0.0;
   
    sumtheta[1][0]=0.0;
    sumtheta[1][1]=0.0;
    
    sumgamma[0]=0.0;
    sumgamma[1]=0.0;
    
   
    vector<vector<vector<double>>>all_eforward{aa_compo_forall.size()};
    
    vector<vector<vector<double>>>all_ebackward{aa_compo_forall.size()};
    
    vector<size_t>all_lenn(aa_compo_forall.size());
    
    vector<vector<vector<double>>>all_gamma{aa_compo_forall.size()};
    vector<vector<vector<vector<double>>>>all_theta{aa_compo_forall.size()};

    vector<vector<int>>all_histate{aa_compo_forall.size()};
    
    vector<vector<vector<double>>>all_sumtheta{aa_compo_forall.size()};
    vector<vector<double>>all_sumgamma{aa_compo_forall.size()};
    
   
    vector<vector<double>> all_loglikelihood(aa_compo_forall.size(),vector<double>(1000));

    
      for(size_t i=0;i<aa_compo_forall.size();i++)
    {
        
        
        size_t lenn=aa_compo_forall[i].size();//get length of each observations
        
        vector<vector<double>>forward_my=init_forward(aa_compo_forall[i], lenn,i, /*hydro_forall[i],pka1_forall[i],helixpro_forall[i],*/ hydro_allstates_emission,pka1_allstates_emission, helixpro_allstates_emission,stericpar_allstates_emission,polar_allstates_emission,volume_allstates_emission, sheetpro_allstates_emission);
        
        
        loglikelihood=get_loglikelihood(lenn, forward_my);
        
        
        all_loglikelihood[i][0]=loglikelihood;
        
        vector<vector<double>>backward_my=init_backward(aa_compo_forall[i], lenn,i, /*hydro_forall[i],pka1_forall[i],helixpro_forall[i],*/ hydro_allstates_emission,pka1_allstates_emission,helixpro_allstates_emission,stericpar_allstates_emission,polar_allstates_emission,volume_allstates_emission, sheetpro_allstates_emission);
        
        vector<vector<double>> eforward=assign_for(forward_my, lenn);
        
        vector<vector<double>>ebackward=assign_back(backward_my, lenn);
        vector<vector<double>> ggamma=get_gamma(lenn);
        
        vector<vector<vector<double>>>theta=get_theta(lenn);
        
        vector<int>histate=get_state(lenn);
        
        
        all_eforward[i]=eforward;
        
        all_ebackward[i]=ebackward;
        
        
        all_gamma[i]=ggamma;
        
        all_theta[i]=theta;
        all_lenn[i]=lenn;
        all_histate[i]=histate;
    }
    vector<int>signs(aa_compo_forall.size());
    
    cerr<<"HMM starts...\n";
    
    
    
    
    
    do
    {
        cout<<c<<endl;
        
        c=c+1;
        
        for(size_t i=0;i<aa_compo_forall.size();i++)
            
        {
            
           
            
            
            all_theta[i]=bw_theta(i,aa_compo_forall[i],/*hydro_forall[i],pka1_forall[i],helixpro_forall[i], */hydro_allstates_emission,pka1_allstates_emission,helixpro_allstates_emission,stericpar_allstates_emission,polar_allstates_emission,volume_allstates_emission, sheetpro_allstates_emission,all_lenn[i],all_eforward[i],all_ebackward[i],all_theta[i]);
            
            
            all_gamma[i]=bw_gamma(all_lenn[i], all_eforward[i], all_gamma[i], all_theta[i]);
            
            
            
            all_histate[i]=get_histates(all_lenn[i], all_gamma[i], all_histate[i],myparavalue.prob_threshold);
            all_sumtheta[i]=sum_theta(all_lenn[i],all_theta[i]);
            all_sumgamma[i]=sum_gamma(all_lenn[i],all_gamma[i]);
            
            
        }
       
        epi=update_epi(all_gamma);
        etransmatrix=update_etransmatirx(all_sumtheta,all_sumgamma);
        
        for(size_t i=0;i<aa_compo_forall.size();i++)
        {
          
            
            all_eforward[i]=update_eforward(aa_compo_forall[i],i,/*hydro_forall[i],pka1_forall[i],  helixpro_forall[i],*/all_lenn[i], hydro_allstates_emission,pka1_allstates_emission,helixpro_allstates_emission,stericpar_allstates_emission,polar_allstates_emission,volume_allstates_emission, sheetpro_allstates_emission, all_eforward[i]);
            
            all_ebackward[i]=update_ebackward(aa_compo_forall[i],i,/*hydro_forall[i],pka1_forall[i],helixpro_forall[i],*/ all_lenn[i], hydro_allstates_emission, pka1_allstates_emission,helixpro_allstates_emission ,stericpar_allstates_emission,polar_allstates_emission,volume_allstates_emission, sheetpro_allstates_emission, all_ebackward[i]);
            
            all_loglikelihood[i][c+1]=get_loglikelihood(all_lenn[i], all_eforward[i]);
            
            
        }
        
        for(size_t i=0;i<aa_compo_forall.size();i++)
        {
            if(abs(all_loglikelihood[i][c+1]-all_loglikelihood[i][c])<myparavalue.stop_criterion)
                signs[i]=0;
            else
                signs[i]=1;
            
        }
        
        
       // cout<<"see if it converges after converging "<<signs[3]<<endl;
        
    }
    
    
    
    while(sumofintegers(signs)>0&&c<20);
    
    cout<<"HMM ends "<<endl;
    
    
    
    
    
    
    
    vector<vector<int>>  psp_windows_states=read_tsv_int(myparavalue.back_states);
    
    
    
   // cout<<"VERY IMPORTANT INFO"<<endl;
    
    cout<<psp_windows_states.size()<<" "<<psp_windows_states[1].size()<<endl;
    
    
    
    
    
    
    vector<vector<int>>modi_histate=all_histate;
    
  
    int truepositive=0;
    int falsepositive=0;
    int truenegative=0;
    int falsenegative=0;

    
    
    
    for (size_t l=0; l<all_histate.size(); l++)
    {
        for (size_t j=0; j<all_histate[l].size(); j++)
        {
           
            int sig=0;
            
            for (size_t f=0; f<myaawhich.size(); f++)
            {
                if (aa_compo_forall[l][j][myaawhich[f]]!=0)
                {
                    sig=1;
                }
            }
            
            if (sig!=1)
            {
                modi_histate[l][j]=0;
            }
            
        }
    }
    
    
    
    
    for (size_t k=0; k<modi_histate.size(); k++)
    {
        for (size_t p=0; p<modi_histate[k].size(); p++)
        {
            if (psp_windows_states[k][p]==1&modi_histate[k][p]==1)
            {
                truepositive+=1;
            }
            else if (psp_windows_states[k][p]==0&modi_histate[k][p]==1)
            {
                falsepositive+=1;
            }
            
            else if (psp_windows_states[k][p]==0&modi_histate[k][p]==0)
            {
                truenegative+=1;
            }
            else if (psp_windows_states[k][p]==1&modi_histate[k][p]==0)
            {
                falsenegative+=1;
            }
            
            
        }
    }
    
    
    
    //this is to program the t test part
    
    
    double tcomp_cutoff=myparavalue.prob_threshold;
    
    
    double alpha=myparavalue.t_signi;
    
    vector<double>each_prot_hydro_sig(aa_compo_forall.size());
    vector<double>each_prot_pka1_sig(aa_compo_forall.size());
    vector<double>each_prot_helix_prob_sig(aa_compo_forall.size());
    vector<double>each_prot_steric_para_sig(aa_compo_forall.size());
    vector<double>each_prot_polar_sig(aa_compo_forall.size());
    vector<double>each_prot_volume_sig(aa_compo_forall.size());
    vector<double>each_prot_sheet_prob_sig(aa_compo_forall.size());
    
    
    
    
    vector<vector<double>>aa_freq_test_all(aa_compo_forall.size());
    
    
    
    
    //Im going to make the indicator for 20 amino acids, hehe
    
    
   // cout<<"proceed to test"<<" "<<tcomp_cutoff<<" "<<alpha<<endl;
    
    for (size_t i=0; i<aa_compo_forall.size(); i++)
        
    {
        
        
        
        
        vector<int>ups;
        vector<int>downs;
        
        
        for (size_t p=0; p<aa_compo_forall[i].size(); p++)
        {
            if (exp(all_gamma[i][p][1])>tcomp_cutoff)
            {
                ups.push_back(p);
            }
            else
            {
                downs.push_back(p);
            }
            
        }
        
        
    
        
        
      //  cout<<i<<" size of ups and size of downs"<<ups.size()<<" "<<downs.size()<<endl;
        
        
        
        if(ups.size()>3&&downs.size()>3)
            
        {
        
        vector<double>hydro_up;
        vector<double>hydro_down;
        vector<double>pka_up;
        vector<double>pka_down;
        vector<double>helix_up;
        vector<double>helix_down;
        vector<double>steric_up;
        vector<double>steric_down;
        vector<double>polar_up;
        vector<double>polar_down;
        vector<double>volume_up;
        vector<double>volume_down;
        vector<double>sheet_up;
        vector<double>sheet_down;

        
        for (size_t u=0; u<ups.size(); u++)
        {
            hydro_up.push_back(hydro_forall[i][ups[u]]);
            pka_up.push_back(pka1_forall[i][ups[u]]);
            helix_up.push_back(helixpro_forall[i][ups[u]]);
            steric_up.push_back(stericpar_forall[i][ups[u]]);
            polar_up.push_back(polar_forall[i][ups[u]]);
            volume_up.push_back(volume_forall[i][ups[u]]);
            sheet_up.push_back(sheetpro_forall[i][ups[u]]);
            
        }
        
        for (size_t d=0; d<downs.size(); d++)
        {
            hydro_down.push_back(hydro_forall[i][downs[d]]);
            pka_down.push_back(pka1_forall[i][downs[d]]);
            helix_down.push_back(helixpro_forall[i][downs[d]]);
            steric_down.push_back(stericpar_forall[i][downs[d]]);
            polar_down.push_back(polar_forall[i][downs[d]]);
            volume_down.push_back(volume_forall[i][downs[d]]);
            sheet_down.push_back(sheetpro_forall[i][downs[d]]);
            
        }

        
        
        
        double hydro_up_mean=getmean(hydro_up);
        double hydro_up_sd=sqrt(getvariance(hydro_up,hydro_up_mean));
        
        double hydro_down_mean=getmean(hydro_down);
        double hydro_down_sd=sqrt(getvariance(hydro_down,hydro_down_mean));
        
        double pka_up_mean=getmean(pka_up);
        double pka_up_sd=sqrt(getvariance(pka_up,pka_up_mean));
        
        double pka_down_mean=getmean(pka_down);
        double pka_down_sd=sqrt(getvariance(pka_down,pka_down_mean));
        
        double helix_up_mean=getmean(helix_up);
        double helix_up_sd=sqrt(getvariance(helix_up,helix_up_mean));
        
        double helix_down_mean=getmean(helix_down);
        double helix_down_sd=sqrt(getvariance(helix_down,helix_down_mean));
        
        double steric_up_mean=getmean(steric_up);
        double steric_up_sd=sqrt(getvariance(steric_up,steric_up_mean));
        
        double steric_down_mean=getmean(steric_down);
        double steric_down_sd=sqrt(getvariance(steric_down,steric_down_mean));
        

        double polar_up_mean=getmean(polar_up);
        double polar_up_sd=sqrt(getvariance(polar_up,polar_up_mean));
        
        double polar_down_mean=getmean(polar_down);
        double polar_down_sd=sqrt(getvariance(polar_down,polar_down_mean));
        
        double volume_up_mean=getmean(volume_up);
        double volume_up_sd=sqrt(getvariance(volume_up,volume_up_mean));
        
        double volume_down_mean=getmean(volume_down);
        double volume_down_sd=sqrt(getvariance(volume_down,volume_down_mean));
        

        double sheet_up_mean=getmean(sheet_up);
        double sheet_up_sd=sqrt(getvariance(sheet_up,sheet_up_mean));
        
        double sheet_down_mean=getmean(sheet_down);
        double sheet_down_sd=sqrt(getvariance(sheet_down,sheet_down_mean));
        
            
          //  cout<<"begin to do ttest "<<endl;
  each_prot_hydro_sig[i]=two_samples_t_test_equal_sd(hydro_up_mean,hydro_up_sd,ups.size(),hydro_down_mean,hydro_down_sd,downs.size(),alpha);
            
            
                      //  cout<<"get hydro..."<<endl;
        each_prot_pka1_sig[i]=two_samples_t_test_equal_sd(pka_up_mean,pka_up_sd,ups.size(),pka_down_mean,pka_down_sd,downs.size(),alpha);
//        
      each_prot_helix_prob_sig[i]=two_samples_t_test_equal_sd(helix_up_mean,helix_up_sd,ups.size(),helix_down_mean,helix_down_sd,downs.size(),alpha);
//
//        
      each_prot_steric_para_sig[i]=two_samples_t_test_equal_sd(steric_up_mean,steric_up_sd,ups.size(),steric_down_mean,steric_down_sd,downs.size(),alpha);
//
//        
      each_prot_polar_sig[i]=two_samples_t_test_equal_sd(polar_up_mean,polar_up_sd,ups.size(),polar_down_mean,polar_down_sd,downs.size(),alpha);
//
      each_prot_volume_sig[i]=two_samples_t_test_equal_sd(volume_up_mean,volume_up_sd,ups.size(),volume_down_mean,volume_down_sd,downs.size(),alpha);
//
      each_prot_sheet_prob_sig[i]=two_samples_t_test_equal_sd(sheet_up_mean,sheet_up_sd,ups.size(),sheet_down_mean,sheet_down_sd,downs.size(),alpha);
//
       
           // cout<<"start to see chi..."<<endl;
            
            vector<vector<int>> each_obs=cal_obs(aa_compo_forall[i],ups,downs);
            
           // cout<<"get each observation for 20 aa "<<endl;
            
            for (size_t l=0; l<20; l++)
            {
                
                
                double chi_stat=cal_chi_stats(each_obs[l][0],each_obs[l][1],each_obs[l][2],each_obs[l][3]);
                
                
                if (chi_stat>0)
                {
                    aa_freq_test_all[i].push_back(cal_pvalue(chi_stat));

                }
                else
                    aa_freq_test_all[i].push_back(8);

                
                
                
                //exam and test and add to protein level summarya
                
                
            }
            
        
            
            
        }
        
        
        else
        {
            
            
            each_prot_hydro_sig[i]=8;
            each_prot_pka1_sig[i]=8;
            each_prot_helix_prob_sig[i]=8;
            each_prot_steric_para_sig[i]=8;
            each_prot_polar_sig[i]=8;
            each_prot_volume_sig[i]=8;
            each_prot_sheet_prob_sig[i]=8;
            
            
            
            
            
            for (size_t l=0; l<20; l++)
            {
                aa_freq_test_all[i].push_back(8);
            }
        }
        
        
       // cout<<"check each freq's size "<<aa_freq_test_all[i].size()<<endl;
        
        
    }

    cout<<"after the testing part..."<<endl;
    
    
    
    
    
    
    
    
    outfile.open("full_output.tsv");
    
    outfile<<"Prot"<<"\t"<<"ID"<<"\t"<<"Seq"<<"\t"<<"Window"<<"\t"<<"Start_position"<<"\t"<<"End_position"<<"\t"<<"State_0"<<"\t"<<"State_1"<<"\t"<<"Optimal_state"<<"\t"<<"Filtered"<<"\t"<<"Observed_state"<<"\t"<<"PSP_state"<<"\t"<<"Hydrophobicity"<<"\t"<<"Pka1"<<"\t"<<"Helix_probability"<<"\t"<<"Steric_parameter"<<"\t"<<"Polar"<<"\t"<<"Volume"<<"\t"<<"Sheet_probability"<<"\t";
    
    for (size_t k=0; k<20; k++)
    {
        outfile<<aminoacid_properties[k].aaname<<"\t";
    }
    
    outfile<<endl;
    
    if (myparavalue.block_overlaps==0)
    {
        for (size_t r=0; r<aa_compo_forall.size(); r++)
        {
            for (size_t g=0; g<aa_compo_forall[r].size();g++ )
            {
                outfile<<r+1<<"\t"<<myuniquesites[r].first.first<<"\t"<<g+1<<"\t"<<marked_window_forall[r][g]<<"\t"<<g*myparavalue.block_length+1<<"\t"<<g*myparavalue.block_length+window_for_all[r][g].size()<<"\t"<<exp(all_gamma[r][g][0])<<"\t"<<exp(all_gamma[r][g][1])<<"\t"<<all_histate[r][g]<<"\t"<<modi_histate[r][g]<<"\t"<<windows_states[r][g]<<"\t"<<psp_windows_states[r][g]<<"\t"<<hydro_forall[r][g]<<"\t"<<pka1_forall[r][g]<<"\t"<<helixpro_forall[r][g]<<"\t"<<stericpar_forall[r][g]<<"\t"<<polar_forall[r][g]<<"\t"<<volume_forall[r][g]<<"\t"<<sheetpro_forall[r][g]<<"\t";
                
                for (size_t l=0; l<20; l++)
                {
                    outfile<<aa_compo_forall[r][g][l]<<"\t";
                }
                
                outfile<<endl;
            }
            
        
            
        }

    }
    else
    {
        for (size_t r=0; r<aa_compo_forall.size(); r++)
        {
            for (size_t g=0; g<aa_compo_forall[r].size()-1;g++ )
            {
                outfile<<r+1<<"\t"<<myuniquesites[r].first.first<<"\t"<<g+1<<"\t"<<marked_window_forall[r][g]<<"\t"<<g*(myparavalue.block_length-myparavalue.block_overlaps)+1<<"\t"<<g*(myparavalue.block_length-myparavalue.block_overlaps)+window_for_all[r][g].size()<<"\t"<<exp(all_gamma[r][g][0])<<"\t"<<exp(all_gamma[r][g][1])<<"\t"<<all_histate[r][g]<<"\t"<<modi_histate[r][g]<<"\t"<<windows_states[r][g]<<"\t"<<psp_windows_states[r][g]<<"\t"<<hydro_forall[r][g]<<"\t"<<pka1_forall[r][g]<<"\t"<<helixpro_forall[r][g]<<"\t"<<stericpar_forall[r][g]<<"\t"<<polar_forall[r][g]<<"\t"<<volume_forall[r][g]<<"\t"<<sheetpro_forall[r][g]<<"\t";
                
                for (size_t l=0; l<20; l++)
                {
                    outfile<<aa_compo_forall[r][g][l]<<"\t";
                }
                
                outfile<<endl;
            }
            
            if (marked_window_forall[r][aa_compo_forall[r].size()-1].size()>=myparavalue.block_overlaps)
            {
                int q=aa_compo_forall[r].size()-1;
                
                outfile<<r+1<<"\t"<<myuniquesites[r].first.first<<"\t"<<q+1<<"\t"<<marked_window_forall[r][q]<<"\t"<<q*(myparavalue.block_length-myparavalue.block_overlaps)+1<<"\t"<<q*(myparavalue.block_length-myparavalue.block_overlaps)+window_for_all[r][q].size()<<"\t"<<exp(all_gamma[r][q][0])<<"\t"<<exp(all_gamma[r][q][1])<<"\t"<<all_histate[r][q]<<"\t"<<modi_histate[r][q]<<"\t"<<windows_states[r][q]<<"\t"<<psp_windows_states[r][q]<<"\t"<<hydro_forall[r][q]<<"\t"<<pka1_forall[r][q]<<"\t"<<helixpro_forall[r][q]<<"\t"<<stericpar_forall[r][q]<<"\t"<<polar_forall[r][q]<<"\t"<<volume_forall[r][q]<<"\t"<<sheetpro_forall[r][q]<<"\t";
                
                for (size_t l=0; l<20; l++)
                {
                    outfile<<aa_compo_forall[r][q][l]<<"\t";
                }
                
                outfile<<endl;

            }
            
            
        }
        
        
    }
    
    
    
    
    outfile.close();
    
    vector<int>each_prot_possiblesites(aa_compo_forall.size(),0);
  
    for (size_t i=0; i<myuniquesites.size(); i++)
    {
        for (size_t f=0; f<myparavalue.pred_sites.size(); f++)
        {
            for (size_t p=0; p<myuniquesites[i].first.second.size(); p++)
            {
                if (myuniquesites[i].first.second[p]==toupper(myparavalue.pred_sites[f]))
                {
                    each_prot_possiblesites[i]+=1;
                }

            }
        }
    }
    
    
    vector<vector<size_t>>sites_positions_inprot(aa_compo_forall.size());
    
    for (size_t i=0; i<aa_compo_forall.size(); i++)
    {
        for (size_t p=0; p<myuniquesites[i].second.size(); p++)
        {
            sites_positions_inprot[i].push_back(myuniquesites[i].second[p].position);
            
           // cout<<i<<" "<<p<<" position of observation "<<myuniquesites[i].second[p].position<<endl;
        }
    }
    
    for (size_t i=0; i<sites_positions_inprot.size(); i++)
    {
        sites_positions_inprot[i]=read_concise_position(sites_positions_inprot[i]);
    }
    
    
    
    
    
    vector<int>each_prot_predictedblocks(all_histate.size(),0);
    for (size_t i=0; i<all_histate.size(); i++)
    {
        for (size_t p=0; p<all_histate[i].size(); p++)
        {
            each_prot_predictedblocks[i]+=all_histate[i][p];

        }
        
    }
    
    vector<int>each_prot_observedblocks(windows_states.size(),0);
    for (size_t i=0; i<windows_states.size(); i++)
    {
        for (size_t p=0; p<windows_states[i].size(); p++)
        {
            each_prot_observedblocks[i]+=windows_states[i][p];
        }
    }
    
    
    vector<int>each_prot_pspblocks(psp_windows_states.size(),0);
    for (size_t i=0; i<psp_windows_states.size(); i++)
    {
        if (psp_windows_states[i][0]!=8)
        {
            for (size_t p=0; p<psp_windows_states[i].size(); p++)
            {
                
                each_prot_pspblocks[i]+=psp_windows_states[i][p];
            }

        }
        else
            each_prot_pspblocks[i]=-1;
        
    }
    
    
    cout<<"generate results..."<<endl;
    
    
    outfile.open("proteinlevel_summary.tsv");
    
        outfile<<"Prot"<<"\t"<<"ID"<<"\t"<<"Protein_length"<<"\t"<<"NO.protein_blocks"<<"\t"<<"NO.observed_blocks"<<"\t"<<"NO.psp_blocks"<<"\t"<<"NO.predicted_blocks"<<"\t"<<"NO.possible_sites"<<"\t"<<"NO.observed_sites"<<"\t"<<"Hydrophobicity"<<"\t"<<"Pka1"<<"\t"<<"Helix_probability"<<"\t"<<"Steric_parameter"<<"\t"<<"Polar"<<"\t"<<"Volume"<<"\t"<<"Sheet_probability"<<"\t";
    
    for (size_t k=0; k<20; k++)
           {
                outfile<<aminoacid_properties[k].aaname<<"\t";
           }
    
        outfile<<endl;
    
    

        for (size_t i=0; i<aa_compo_forall.size(); i++)
        {
            outfile<<i+1<<"\t"<<myuniquesites[i].first.first<<"\t"<<myuniquesites[i].first.second.size()<<"\t"<<aa_compo_forall[i].size()<<"\t"<<each_prot_observedblocks[i]<<"\t"<<each_prot_pspblocks[i]<<"\t"<<each_prot_predictedblocks[i]<<"\t"<<each_prot_possiblesites[i]<<"\t"<<sites_positions_inprot[i].size()<<"\t"<<each_prot_hydro_sig[i]<<"\t"<<each_prot_pka1_sig[i]<<"\t"<<each_prot_helix_prob_sig[i]<<"\t"<<each_prot_steric_para_sig[i]<<"\t"<<each_prot_polar_sig[i]<<"\t"<<each_prot_volume_sig[i]<<"\t"<<each_prot_sheet_prob_sig[i]<<"\t";
        
    
    
    
                for (size_t l=0; l<20; l++)
                {
                       outfile<<aa_freq_test_all[i][l]<<"\t";
                    }

            outfile<<endl;

        }
    
    
    outfile.close();
    
    
    cout<<double(truepositive)/double(truepositive+falsenegative)<<endl;
    
    
    cout<<double(truenegative)/double(truenegative+falsepositive)<<endl;
    
 
    
    cout<<double(truepositive)/double(truepositive+falsepositive)<<endl;
    
    
    cout<<double(truenegative)/double(truenegative+falsenegative)<<endl;
    
    
    cout<<double(truenegative+truepositive)/double(truenegative+truepositive+falsepositive+falsenegative)<<endl;
    
    
    cout<<truepositive<<" "<<falsepositive<<" "<<truenegative<<" "<<falsenegative<<endl;
    
    
    outfile.open("reduced_output.tsv");
    
    outfile<<"Prot"<<"\t"<<"ID"<<"\t"<<"Seq"<<"\t"<<"State_1"<<"\t"<<"Observed_state"<<"\t"<<"Psp_state"<<"\t"<<"Filtered_state"<<"\t"<<"Hydrophobicity"<<"\t"<<"Pka1"<<"\t"<<"Helix_probability"<<"\t"<<"Steric_parameter"<<"\t"<<"Polarizability"<<"\t"<<"Volume"<<"\t"<<"Sheet_probability"<<endl;
    
    
    
    
    for (size_t r=0; r<aa_compo_forall.size(); r++)
    {
        outfile<<r+1<<"\t"<<myuniquesites[r].first.first<<"\t"<<1<<"\t"<<exp(all_gamma[r][0][1])<<"\t"<<windows_states[r][0]<<"\t"<<psp_windows_states[r][0]<<"\t"<<modi_histate[r][0]<<"\t"<<hydro_forall[r][0]<<"\t"<<pka1_forall[r][0]<<"\t"<<helixpro_forall[r][0]<<"\t"<<stericpar_forall[r][0]<<"\t"<<polar_forall[r][0]<<"\t"<<volume_forall[r][0]<<"\t"<<sheetpro_forall[r][0]<<endl;
        

        
        
        for (size_t g=1; g<aa_compo_forall[r].size();g++ )
        {
            outfile<<r+1<<"\t"<<"OM"<<"\t"<<g+1<<"\t"<<exp(all_gamma[r][g][1])<<"\t"<<windows_states[r][g]<<"\t"<<psp_windows_states[r][g]<<"\t"<<modi_histate[r][g]<<"\t"<<hydro_forall[r][g]<<"\t"<<pka1_forall[r][g]<<"\t"<<helixpro_forall[r][g]<<"\t"<<stericpar_forall[r][g]<<"\t"<<polar_forall[r][g]<<"\t"<<volume_forall[r][g]<<"\t"<<sheetpro_forall[r][g]<<endl;
            
            
           
        }
        
       }
    outfile.close();
    
    
    outfile.open("s_ID.tsv");
    for (size_t i=0; i<aa_compo_forall.size(); i++)
    {
        outfile<<i+1<<"\t"<<myuniquesites[i].first.first<<endl;
    }
    outfile.close();
    
    
    
    
    //this below are what I added.
    
    
    vector<vector<pair<int,int>>>predstates_aafreq{2};
    
    
    predstates_aafreq=get_statesrates(aa_compo_forall,modi_histate);
    
    
    
   // cout<<"get predstates_aafreq "<<predstates_aafreq.size()<<" "<<predstates_aafreq[1].size()<<" "<<predstates_aafreq[0][1].first<<" "<<predstates_aafreq[0][1].second<<endl;
    
    
    vector<vector<double>> eachstate_pred_aafreq{2};
    
    
    
    assign_aaemissions(predstates_aafreq,eachstate_pred_aafreq);
    
    
    
   // cout<<"see how it is assigned "<<eachstate_pred_aafreq.size()<<" "<<eachstate_pred_aafreq[0].size()<<" "<<eachstate_pred_aafreq[0][1]<<endl;
    
    
    
    outfile.open("aa_pred_states_output.tsv");
    
    outfile<<"AA"<<"\t"<<"State_0 freq"<<"\t"<<"State_1 freq"<<"\t"<<"Ratio"<<endl;
    
    for (size_t i=0; i<20; i++)
    {
        cout<<i<<endl;
        outfile<<aminoacid_properties[i].aaname<<"\t"<<eachstate_pred_aafreq[0][i]<<"\t"<<eachstate_pred_aafreq[1][i]<<"\t"<<double(eachstate_pred_aafreq[1][i])/double(eachstate_pred_aafreq[0][i])<<endl;
    }
    
    outfile.close();
    
    
    
    
    
    
}

