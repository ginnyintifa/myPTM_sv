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
#include "input_external_information.hpp"

int main(int argc, char *argv[])

{
    
    using namespace std;
    
    //there is one modification on mypeptides_site, change it to myfirst_sites for convinience
    
    vector<string>myinput=read_input(argv[1]);
    vector<string>mypara=get_parameterlines(myinput);
    
    Para_values myparavalue=get_paravalues(mypara);
    vector<int> myaawhich=get_aawhich(myparavalue.pred_sites);
    
    
    vector<prot_type> myproteins=read_fasta(myparavalue.fasta_name);//"uniprot-SGD.fasta");
    
    
    myproteins=read_concise_pair(myproteins);

    
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
    vector<string> myexpinfo=read_tsv(myparavalue.expinfo);
    
    
    
    
    cout<<myexpinfo.size()<<endl;
    for (size_t i=0; i<myexpinfo.size(); i++)
    {
        for (size_t p=0; p<myexpinfo[i].size(); p++)
        {
            int thesign=0;
            
            for (size_t q=0; q<myparavalue.pred_sites.size(); q++)
            {
                if (myexpinfo[i][p]==myparavalue.pred_sites[q])
                {
                    thesign=1;
                }
                
                
            }
            
            if (thesign==0)
            {
                myexpinfo[i][p]=toupper(myexpinfo[i][p]);
            }
        }
    }
    
    
    cout<<"get signs...\n";
    
    
    vector<string>myexpinfopeptides=myexpinfo;
    
    
    for (size_t i=0; i<myexpinfo.size(); i++)
    {
        for (size_t p=0; p<myexpinfo[i].size(); p++)
        {
            myexpinfopeptides[i][p]=toupper(myexpinfopeptides[i][p]);
        }
    }
    
    cout<<"get the peptides...\n";
    
    
    myexpinfopeptides=read_concise(myexpinfopeptides);//this is the phospho data
    
    
    vector<Pep_in_prot> myexplinks=link_peptides_with_sites(myexpinfopeptides, myexpinfo,myparavalue.pred_sites);
    cout<<"get links...\n";
    
    
    //do the mapping thing;
    
    pair<vector<prot_type>,vector<vector<Pep_in_prot>>>expppmap=get_mapps(myproteins, myexplinks);
    
    cout<<"get maps...\n";
    vector<pair<prot_type,vector<site_type>>> myexpprotsitesall=get_all_seq_sites(expppmap);
    vector<pair<prot_type,vector<site_type>>> myexpneededprotsiteall=get_needed_myprot_sites_all(myexpprotsitesall);
    
    
    vector<pair<prot_type,vector<site_type>>> myexpuniquesites=get_unique_site(myexpneededprotsiteall);
    
    
    

    
    
    
    sort(myexpuniquesites.begin(), myexpuniquesites.end(), sortpair);
    
    
    cout<<myexpuniquesites.size()<<endl;
    

    
    vector<vector<string>> expwindow_for_all=get_windows_for_all(myexpuniquesites,myparavalue.block_length,myparavalue.block_overlaps);

    
    vector<prot_type>proteinsnow;
    
    for (size_t i=0; i<myexpuniquesites.size(); i++)
    {
        proteinsnow.push_back(myexpuniquesites[i].first);
    }
    
    cout<<proteinsnow.size()<<endl;
    
    
    
    
    
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
        
        mybackinfopeptides=read_concise(mybackinfopeptides);//this is the phospho data
        
        
        vector<Pep_in_prot> mybacklinks=link_peptides_with_sites(mybackinfopeptides, mybackinfo,myparavalue.pred_sites);
        cout<<"get links...\n";
    
    
    //do the mapping thing;
        
        pair<vector<prot_type>,vector<vector<Pep_in_prot>>>backppmap=get_mapps(proteinsnow, mybacklinks);
        
        cout<<"get maps...\n";
    
    
        vector<pair<prot_type,vector<site_type>>> mybackprotsitesall=get_all_seq_sites(backppmap);
        vector<pair<prot_type,vector<site_type>>> mybackneededprotsiteall=get_needed_myprot_sites_all(mybackprotsitesall);
    
    
        vector<pair<prot_type,vector<site_type>>> mybackuniquesites=get_unique_site(mybackneededprotsiteall);
    
    
    //check what they do on not mapped sequence
    
    for (size_t i=0; i<mybackuniquesites.size(); i++)
    {
        cout<<i<<" each mapp's size "<<mybackuniquesites[i].second.size()<<endl;
    }
    
    
    //double check;
    cout<<mybackuniquesites[0].first.first<<endl;
        
       // sort(mybackuniquesites.begin(), mybackuniquesites.end(), sortpair);
    
    
    
    
    
        vector<vector<string>> backwindow_for_all=get_windows_for_all(mybackuniquesites,myparavalue.block_length,myparavalue.block_overlaps);
        
        
        //vector<vector<vector<int>>> backaa_compo_forall=get_allseqs__aavector(backwindow_for_all, aminoacid_properties);
        
       // vector<vector<int>>backaa_vectorforall;//for the whole sequence
        
        //for (size_t i=0; i<mybackuniquesites.size(); i++)
        //{
          //  vector<int> backaa_vectorforone;
            
            //backaa_vectorforone=grab_aa_vector(aminoacid_properties, mybackuniquesites[i].first.second);
            
            //backaa_vectorforall.push_back(backaa_vectorforone);
        //}
    
    
    
    
        vector<pair<prot_type,vector<vector<site_type>>>>backwindow_uniquesites_forall=get_unique_site_inblock(mybackuniquesites,myparavalue.block_length,myparavalue.block_overlaps);
        
       // vector<vector<string>>backmarked_window_forall=marked_windows(backwindow_uniquesites_forall, backwindow_for_all,myparavalue.block_length,myparavalue.block_overlaps);
        vector<vector<int>>backwindows_states=get_states(backwindow_uniquesites_forall,myparavalue.pred_sites);
        
        

    
       // vector<vector<pair<int, int>>>backstates_windowrates=get_statesrates(backaa_compo_forall, backwindows_states);
       // vector<pair<prot_type,vector<int>>>backallseqs_blockwithsites=get_blocks_withsites_all(backwindow_uniquesites_forall);
    
    
    
    
    ofstream outfile;
    
    
    //actually I want to output it in the order of the exp sequence , if no just set all as NAN yes maybe -1
    
    vector<vector<int>>mapped_backwindows_states{expwindow_for_all.size()};

    
    for (size_t i=0; i<expwindow_for_all.size(); i++)
    {
        //cout<<expwindow_for_all[i].size();
        
        for (size_t p=0; p<expwindow_for_all[i].size(); p++)
        {
            
            
            mapped_backwindows_states[i].push_back(8);
        }
    }
    
    
   // cout<<mapped_backwindows_states.size()<<"############ "<<backwindows_states.size()<<"######## "<<proteinsnow.size()<<endl;
    
    
    for (size_t i=0;i<mapped_backwindows_states.size() ;i++)
    {
        for (size_t h=0; h<backwindow_uniquesites_forall.size(); h++)
        {
            
            
            if (backwindow_uniquesites_forall[h].first.first==proteinsnow[i].first)
            {
                
                cout<<h<<" "<<backwindow_uniquesites_forall[h].first.first<<endl;
                
                if (mapped_backwindows_states[i].size()!=backwindows_states[h].size())
                {
                    cout<<mapped_backwindows_states[i].size()<<" "<<backwindows_states[h].size()<<endl;
                }
                
                for (size_t u=0; u<mapped_backwindows_states[i].size(); u++)
                {
                    cout<<backwindows_states[h][u]<<" ";
                    
                    mapped_backwindows_states[i][u]=backwindows_states[h][u];

                }
                
                cout<<endl;
                //can add break here;
                
            }
            
            
        }
    }
    
    
    
    
    
    
    
    
    
    outfile.open("back_states.tsv");
    
    for (size_t i=0; i<mapped_backwindows_states.size(); i++)
    {
        for (size_t p=0; p<mapped_backwindows_states[i].size(); p++)
        {
            outfile<<mapped_backwindows_states[i][p]<<"\t";
        }
        
        outfile<<endl;
    }
    
    outfile.close();
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
}

