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

Double_t getIndex(vector<Double_t> v, Double_t K) 
{ 
    auto it = find(v.begin(), v.end(), K); 
  
    // If element was found 
    if (it != v.end())  
    { 
      
        // calculating the index 
        // of K 
        Double_t index = it - v.begin(); 
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
    std::vector<Double_t> *fGeTetra_E_coinc = 0;
    std::vector<Double_t> *fGeTetra_t_coinc = 0;
    std::vector<Double_t> *fGeTetra_Index = 0;
    std::vector<Double_t> *fGeBeta_Index = 0;

    std::vector<Double_t> fTripleCoinc_ECond;
    std::vector<Double_t> fTripleCoinc_tCond;

    TFile* data_file = TFile::Open(inputFileName, "UPDATE");

    TTree* data_tree = (TTree*)data_file->Get("tcoinc");

    TTree* triple_coinc = new TTree("triple_coinc", "triple_coinc");

    TBranch* branch_E = triple_coinc->Branch("TripleCoinc_ECond", &fTripleCoinc_ECond);
    TBranch* branch_t = triple_coinc->Branch("TripleCoinc_tCond", &fTripleCoinc_tCond);
    
    data_tree->SetBranchAddress("GeTetra_ECond", &fGeTetra_E_coinc);
    data_tree->SetBranchAddress("GeTetra_tCond", &fGeTetra_t_coinc);
    data_tree->SetBranchAddress("GeTetra_Index", &fGeTetra_Index);
    data_tree->SetBranchAddress("GeBeta_Index", &fGeBeta_Index);

    TH1D *GammaBetaNeutronECoinc = new TH1D("GammaBetaNeutronECoinc", "GammaBetaNeutronECoinc", 7000, 0, 7000);
    TH1D *GammaBetaNeutrontCoinc = new TH1D("GammaBetaNeutrontCoinc", "GammaBetaNeutrontCoinc", 10000, 0, 10000);

    GammaBetaNeutronECoinc->Reset();
    GammaBetaNeutrontCoinc->Reset();

    Double_t entries = data_tree->GetEntries();

    std::vector<vector<Double_t>> listTetra;
    std::vector<Double_t> listBeta;

    listTetra.push_back(vector<Double_t>(0));
    listTetra.push_back(vector<Double_t>(0));
    listTetra.push_back(vector<Double_t>(0));

    for(Double_t i = 0; i < entries; i++)
    {   
        data_tree->GetEntry(i);

        Double_t sizeTetra = fGeTetra_Index->size();
        Double_t sizeBeta = fGeBeta_Index->size();

        for(Double_t j = 0; j < sizeTetra; ++j)
        {
            listTetra[0].push_back(fGeTetra_Index->at(j));
            listTetra[1].push_back(fGeTetra_E_coinc->at(j));
            listTetra[2].push_back(fGeTetra_t_coinc->at(j));
        }

        for(Double_t k = 0; k < sizeBeta; ++k)
        {
            listBeta.push_back(fGeBeta_Index->at(k));
        }
    }

    std::reverse(listTetra[0].begin(), listTetra[0].end());
    std::reverse(listTetra[1].begin(), listTetra[1].end());
    std::reverse(listTetra[2].begin(), listTetra[2].end());

    std::vector<Double_t> true_inter;

    vector<Double_t> intersection(listTetra[0].size() + listBeta.size()); 
    vector<Double_t>::iterator it, st; 
  
    it = set_intersection(listTetra[0].begin(), 
                          listTetra[0].end(), 
                          listBeta.begin(), 
                          listBeta.end(), 
                          intersection.begin()); 

    for(st = intersection.begin(); st != it; ++st)
    {
        true_inter.push_back(*st);
    }

    std::vector<Double_t> indexes;

    for(Double_t i = 0; i < true_inter.size(); ++i)
    {
        indexes.push_back(getIndex(listTetra[0], true_inter[i]));
    }

    for(Double_t i = 0; i < indexes.size(); ++i)
    {
        fTripleCoinc_ECond.clear();
        fTripleCoinc_tCond.clear();

        fTripleCoinc_ECond.push_back(listTetra[1][indexes[i]]);
        fTripleCoinc_tCond.push_back(listTetra[2][indexes[i]]);
        GammaBetaNeutronECoinc->Fill(listTetra[1][indexes[i]]);
        GammaBetaNeutrontCoinc->Fill(listTetra[2][indexes[i]]);

        triple_coinc->Fill();
    }

    data_file->Write();
}