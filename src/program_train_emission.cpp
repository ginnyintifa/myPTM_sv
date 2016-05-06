//
//  main.cpp
//  humansites
//
//  Created by Ginny Li on 15-7-28.
//  Copyright (c) 2015年 Ginny. All rights reserved.
//

//a modification
//in R want to label ones that are filtered to be 1 as green,



// a modification of my program

//change to only hydro and pka1



//I mistakenly made some wrong changes


#include "site.hpp"
#include "aaprop.hpp"
#include "emisson.hpp"

#include "math.hpp"
#include "input_train_emission.hpp"

int main(int argc, char *argv[])

{
    
    using namespace std;
    
    
    
    vector<string>myinput=read_input(argv[1]);
    vector<string>mypara=get_parameterlines(myinput);
    
    Para_values myparavalue=get_paravalues(mypara);
    vector<int> myaawhich=get_aawhich(myparavalue.pred_sites);
    
    
    vector<prot_type> myproteins=read_fasta(myparavalue.fasta_name);//function using here is dependent on the format of file.
    
    
    cout<<"see how many proteins i have "<<myproteins.size()<<endl;
    
    vector<aa_properties> aminoacid_properties(24);
    
    assign_properties(aminoacid_properties);
    
    vector<double>each_hydro;
    vector<double>each_pka1;
    vector<double>each_helixpro;
    vector<double>each_stericpar;
    vector<double>each_polar;
    vector<double>each_volume;
    vector<double>each_sheetpro;
    
    
    assign_each_property(aminoacid_properties,each_hydro,each_pka1,each_helixpro,each_stericpar,each_polar,each_volume,each_sheetpro);

    
  

        cout<<myparavalue.backinfo<<endl;
        vector<string> mybackinfo=read_tsv(myparavalue.backinfo);
        
        cout<<mybackinfo.size()<<endl;
        for (size_t i=0; i<mybackinfo.size(); i++)
        {
            for (size_t p=0; p<mybackinfo[i].size(); p++)
            {
                int thesign=0;
                
                for (size_t q=0; q<myparavalue.pred_sites.size(); q++)
                {
                    if (mybackinfo[i][p]==myparavalue.pred_sites[q])
                    {
                        thesign=1;
                    }
                    
                    
                }
                
                if (thesign==0)
                {
                    mybackinfo[i][p]=toupper(mybackinfo[i][p]);
                }
            }
        }
        
        
        cout<<"get signs...\n";
        
        
        vector<string>mybackinfopeptides=mybackinfo;
        
        for (size_t i=0; i<mybackinfo.size(); i++)
        {
            for (size_t p=0; p<mybackinfo[i].size(); p++)
            {
                mybackinfopeptides[i][p]=toupper(mybackinfopeptides[i][p]);
            }
        }
        
        cout<<"get the peptides...\n";
        
        myproteins=read_concise_pair(myproteins);
        
        mybackinfopeptides=read_concise(mybackinfopeptides);
        
        
        vector<Pep_in_prot> mybacklinks=link_peptides_with_sites(mybackinfopeptides, mybackinfo,myparavalue.pred_sites);
        cout<<"get links...\n";
        
        pair<vector<prot_type>,vector<vector<Pep_in_prot>>>backppmap=get_mapps(myproteins, mybacklinks);
        
        cout<<"get maps...\n";
        vector<pair<prot_type,vector<site_type>>> mybackprotsitesall=get_all_seq_sites(backppmap);
        vector<pair<prot_type,vector<site_type>>> mybackneededprotsiteall=get_needed_myprot_sites_all(mybackprotsitesall);
        vector<pair<prot_type,vector<site_type>>> mybackuniquesites=get_unique_site(mybackneededprotsiteall);
        
        sort(mybackuniquesites.begin(), mybackuniquesites.end(), sortpair);
        
        vector<vector<string>> backwindow_for_all=get_windows_for_all(mybackuniquesites,myparavalue.block_length,myparavalue.block_overlaps);
        
        
        vector<vector<vector<int>>> backaa_compo_forall=get_allseqs__aavector(backwindow_for_all, aminoacid_properties);
        
        vector<vector<int>>backaa_vectorforall;//for the whole sequence
        
        for (size_t i=0; i<mybackuniquesites.size(); i++)
        {
            vector<int> backaa_vectorforone;
            
            backaa_vectorforone=grab_aa_vector(aminoacid_properties, mybackuniquesites[i].first.second);
            
            backaa_vectorforall.push_back(backaa_vectorforone);
        }
    
       vector<double>backseqs_means_hydro=get_means(each_hydro,backaa_vectorforall);
    
        vector<double>backseqs_means_pka1=get_means(each_pka1,backaa_vectorforall);
    
        vector<double>backseqs_means_helixpro=get_means(each_helixpro,backaa_vectorforall);
    
        vector<double>backseqs_means_stericpar=get_means(each_stericpar,backaa_vectorforall);
    
        vector<double>backseqs_means_polar=get_means(each_polar,backaa_vectorforall);
    
        vector<double>backseqs_means_volume=get_means(each_volume,backaa_vectorforall);
    
        vector<double>backseqs_means_sheetpro=get_means(each_sheetpro,backaa_vectorforall);
    
    
    

    
    
        vector<pair<prot_type,vector<vector<site_type>>>>backwindow_uniquesites_forall=get_unique_site_inblock(mybackuniquesites,myparavalue.block_length,myparavalue.block_overlaps);
        
        vector<vector<string>>backmarked_window_forall=marked_windows(backwindow_uniquesites_forall, backwindow_for_all,myparavalue.block_length,myparavalue.block_overlaps);
        vector<vector<int>>backwindows_states=get_states(backwindow_uniquesites_forall,myparavalue.pred_sites);
        
        
        
        
        vector<vector<double>> backhydro_forall=get_allseqs_property(backaa_compo_forall,each_hydro,backseqs_means_hydro);//normalized
    
        vector<vector<double>> backpka1_forall=get_allseqs_property(backaa_compo_forall,each_pka1,backseqs_means_pka1);
    
        vector<vector<double>> backhelixpro_forall=get_allseqs_property(backaa_compo_forall,each_helixpro,backseqs_means_helixpro);
    vector<vector<double>> backstericpar_forall=get_allseqs_property(backaa_compo_forall,each_stericpar,backseqs_means_hydro);//normalized
    vector<vector<double>> backpolar_forall=get_allseqs_property(backaa_compo_forall,each_polar,backseqs_means_pka1);
    vector<vector<double>> backvolume_forall=get_allseqs_property(backaa_compo_forall,each_volume,backseqs_means_helixpro);
     vector<vector<double>> backsheetpro_forall=get_allseqs_property(backaa_compo_forall,each_sheetpro,backseqs_means_helixpro);

    
        vector<vector<pair<int, int>>>backstates_windowrates=get_statesrates(backaa_compo_forall, backwindows_states);
        vector<pair<prot_type,vector<int>>>backallseqs_blockwithsites=get_blocks_withsites_all(backwindow_uniquesites_forall);
    
    
    
    
    
    
    vector<vector<double>>backhydro_2states=get_properties_2states(backhydro_forall,backwindows_states);
        
        vector<vector<double>>backpka1_2states=get_properties_2states(backpka1_forall,backwindows_states);
    //Xchange
    
        vector<vector<double>>backhelixpro_2states=get_properties_2states(backhelixpro_forall,backwindows_states);
    
    vector<vector<double>>backstericpar_2states=get_properties_2states(backstericpar_forall,backwindows_states);
    
    vector<vector<double>>backpolar_2states=get_properties_2states(backpolar_forall,backwindows_states);
    //Xchange
    
    vector<vector<double>>backvolume_2states=get_properties_2states(backvolume_forall,backwindows_states);
    vector<vector<double>>backsheetpro_2states=get_properties_2states(backsheetpro_forall,backwindows_states);
    
    
    
    
    
    vector<vector<double>> eachstate_20aa_emi{2};
    
    vector<vector<pair<double,double>>>hydro_emission;
    vector<vector<pair<double,double>>>pka1_emission;
    vector<vector<pair<double,double>>>helixpro_emission;
    vector<vector<pair<double,double>>>stericpar_emission;
    vector<vector<pair<double,double>>>polar_emission;
    vector<vector<pair<double,double>>>volume_emission;
    vector<vector<pair<double,double>>>sheetpro_emission;

    
    
        assign_aaemissions(backstates_windowrates,eachstate_20aa_emi);//note this!!!!
        
        
        backhydro_2states=choose_data(backhydro_2states);
        backpka1_2states=choose_data(backpka1_2states);
    //Xchange
        backhelixpro_2states=choose_data(backhelixpro_2states);
    backstericpar_2states=choose_data(backstericpar_2states);
    backpolar_2states=choose_data(backpolar_2states);
    backvolume_2states=choose_data(backvolume_2states);
    backsheetpro_2states=choose_data(backsheetpro_2states);
    
    
    
    
    
 
        
        
        hydro_emission=get_emission(backhydro_2states);//note this!!!didn't change name of emi
        pka1_emission=get_emission(backpka1_2states);//note this!!!
    //Xchange
        helixpro_emission=get_emission(backhelixpro_2states);
    stericpar_emission=get_emission(backstericpar_2states);
    polar_emission=get_emission(backpolar_2states);
    volume_emission=get_emission(backvolume_2states);
    sheetpro_emission=get_emission(backsheetpro_2states);
        
        //this thing to be noticed
        
        cerr<<"prob calculated...\n";
        
        //something wrong here, I should use the right peptide
        
        
    
    
    ofstream outfile;
    
    
    //think about this
    outfile.open("hydro_emi.tsv");
    
        for (size_t i=0 ; i<hydro_emission[1].size(); i++)
        {
            outfile<<hydro_emission[0][i].first<<endl;
        }
    
    for (size_t i=0 ; i<hydro_emission[1].size(); i++)
    {
        outfile<<hydro_emission[0][i].second<<endl;
    }
    
    for (size_t i=0 ; i<hydro_emission[1].size(); i++)
    {
        outfile<<hydro_emission[1][i].second<<endl;
    }
  
    outfile.close();
    
   
    
    
    outfile.open("pka1_emi.tsv");
    for (size_t i=0 ; i<pka1_emission[1].size(); i++)
    {
        outfile<<pka1_emission[0][i].first<<endl;
    }
    
    for (size_t i=0 ; i<pka1_emission[1].size(); i++)
    {
        outfile<<pka1_emission[0][i].second<<endl;
    }
    
    for (size_t i=0 ; i<pka1_emission[1].size(); i++)
    {
        outfile<<pka1_emission[1][i].second<<endl;
    }
    
    outfile.close();
    
    
    
    outfile.open("helixpro_emi.tsv");
    for (size_t i=0 ; i<helixpro_emission[1].size(); i++)
    {
        outfile<<helixpro_emission[0][i].first<<endl;
    }
    
    for (size_t i=0 ; i<helixpro_emission[1].size(); i++)
    {
        outfile<<helixpro_emission[0][i].second<<endl;
    }
    
    for (size_t i=0 ; i<helixpro_emission[1].size(); i++)
    {
        outfile<<helixpro_emission[1][i].second<<endl;
    }
    
    outfile.close();
    
    
    outfile.open("stericpar_emi.tsv");
    for (size_t i=0 ; i<stericpar_emission[1].size(); i++)
    {
        outfile<<stericpar_emission[0][i].first<<endl;
    }
    
    for (size_t i=0 ; i<stericpar_emission[1].size(); i++)
    {
        outfile<<stericpar_emission[0][i].second<<endl;
    }
    
    for (size_t i=0 ; i<stericpar_emission[1].size(); i++)
    {
        outfile<<stericpar_emission[1][i].second<<endl;
    }
    
    outfile.close();
    
    outfile.open("polar_emi.tsv");
    for (size_t i=0 ; i<polar_emission[1].size(); i++)
    {
        outfile<<polar_emission[0][i].first<<endl;
    }
    
    for (size_t i=0 ; i<polar_emission[1].size(); i++)
    {
        outfile<<polar_emission[0][i].second<<endl;
    }
    
    for (size_t i=0 ; i<polar_emission[1].size(); i++)
    {
        outfile<<polar_emission[1][i].second<<endl;
    }
    
    outfile.close();
    
    
    outfile.open("volume_emi.tsv");
    for (size_t i=0 ; i<volume_emission[1].size(); i++)
    {
        outfile<<volume_emission[0][i].first<<endl;
    }
    
    for (size_t i=0 ; i<volume_emission[1].size(); i++)
    {
        outfile<<volume_emission[0][i].second<<endl;
    }
    
    for (size_t i=0 ; i<volume_emission[1].size(); i++)
    {
        outfile<<volume_emission[1][i].second<<endl;
    }
    
    outfile.close();
    

    outfile.open("sheetpro_emi.tsv");
    for (size_t i=0 ; i<sheetpro_emission[1].size(); i++)
    {
        outfile<<sheetpro_emission[0][i].first<<endl;
    }
    
    for (size_t i=0 ; i<sheetpro_emission[1].size(); i++)
    {
        outfile<<sheetpro_emission[0][i].second<<endl;
    }
    
    for (size_t i=0 ; i<sheetpro_emission[1].size(); i++)
    {
        outfile<<sheetpro_emission[1][i].second<<endl;
    }
    
    outfile.close();
    
       outfile.open("aafreq_emi.tsv");
    
        
    for (size_t k=0; k<20; k++)
    {
        outfile<<eachstate_20aa_emi[0][k]<<endl;
       
    }
    
    
    
    for (size_t k=0; k<20; k++)
    {
        outfile<<eachstate_20aa_emi[1][k]<<endl;
        
    }
    
        outfile.close();
    
    
    cerr<<"emission got...\n";

    
    
    
}

