#include <vector>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TStopwatch.h"
#include "TVector3.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TROOT.h"
#include "TDirectory.h"
#include "TCutG.h"
#include "TRandom.h"

UInt_t getIndex(vector<UInt_t> v, UInt_t K) 
{ 
    auto it = find(v.begin(), v.end(), K); 
  
    // If element was found 
    if (it != v.end())  
    { 
      
        // calculating the index 
        // of K 
        UInt_t index = it - v.begin(); 
        return index;
    } 
    else { 
        // If the element is not 
        // present in the vector 
        return -1;
    } 
} 

void TripleCoinc(const char* inputFileName)
{
    std::vector<UInt_t> *fGeTetra_E_coinc = 0;
    std::vector<UInt_t> *fGeTetra_Index = 0;
    std::vector<UInt_t> *fGeBeta_Index = 0;

    TFile* data_file = TFile::Open(inputFileName, "UPDATE");

    TTree* data_tree = (TTree*)data_file->Get("tcoinc");
    
    data_tree->SetBranchAddress("GeTetra_E_coinc", &fGeTetra_E_coinc);
    data_tree->SetBranchAddress("GeTetra_Index", &fGeTetra_Index);
    data_tree->SetBranchAddress("GeBeta_Index", &fGeBeta_Index);

    TH1I *GammaBetaNeutronCoinc = new TH1I("GammaBetaNeutronCoinc", "GammaBetaNeutronCoinc", 7000, 0, 7000);

    UInt_t entries = data_tree->GetEntries();

    std::vector<vector<UInt_t>> listTetra;
    std::vector<UInt_t> listBeta;

    listTetra.push_back(vector<UInt_t>(0));
    listTetra.push_back(vector<UInt_t>(0));

    for(UInt_t i = 0; i < entries; i++)
    {   
        data_tree->GetEntry(i);

        UInt_t sizeTetra = fGeTetra_Index->size();
        UInt_t sizeBeta = fGeBeta_Index->size();

        for(UInt_t j = 0; j < sizeTetra; ++j)
        {
            listTetra[0].push_back(fGeTetra_Index->at(j));
            listTetra[1].push_back(fGeTetra_E_coinc->at(j));
        }

        for(UInt_t k = 0; k < sizeBeta; ++k)
        {
            listBeta.push_back(fGeBeta_Index->at(k));
        }
    }

    std::reverse(listTetra[0].begin(), listTetra[0].end());
    std::reverse(listTetra[1].begin(), listTetra[1].end());

    std::vector<UInt_t> true_inter;

    vector<UInt_t> intersection(listTetra[0].size() + listBeta.size()); 
    vector<UInt_t>::iterator it, st; 
  
    it = set_intersection(listTetra[0].begin(), 
                          listTetra[0].end(), 
                          listBeta.begin(), 
                          listBeta.end(), 
                          intersection.begin()); 

    for(st = intersection.begin(); st != it; ++st)
    {
        true_inter.push_back(*st);
    }

    std::vector<UInt_t> indexes;
    std::vector<UInt_t> TripleCoincEnergies;

    for(UInt_t i = 0; i < true_inter.size(); ++i)
    {
        indexes.push_back(getIndex(listTetra[0], true_inter[i]));
    }

    for(UInt_t i = 0; i < indexes.size(); ++i)
    {
        TripleCoincEnergies.push_back(listTetra[1][indexes[i]]);
        GammaBetaNeutronCoinc->Fill(listTetra[1][indexes[i]]);
    }

GammaBetaNeutronCoinc->Write();

}