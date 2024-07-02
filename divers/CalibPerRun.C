#include <vector>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>

#include <root/RtypesCore.h>
#include <root/TFile.h>
#include <root/TTree.h>
#include <root/TROOT.h>
#include <root/TH1.h>
#include <root/TH2.h>
#include <root/TRandom.h>
#include <root/TF1.h>
#include <root/TCanvas.h>

void CalibPerRun(const char* InputFileName, const char* OutputFileName, Double_t a, Double_t b)
{   
    TFile* data_file = TFile::Open(InputFileName);

    TTree* tree_triple_coinc = (TTree*)data_file->Get("triple_coinc");
    TTree* tree_double_coinc = (TTree*)data_file->Get("tcoinc");

    std::vector<Double_t> *ffTripleCoinc_ECond = 0;
    std::vector<Double_t> *ffGeBeta_ECond = 0;

    tree_triple_coinc->SetBranchAddress("TripleCoinc_ECond", &ffTripleCoinc_ECond);
    tree_double_coinc->SetBranchAddress("GeBeta_ECond", &ffGeBeta_ECond);

    TH1D* hist_triple = new TH1D("hist_triple", "hist_triple", 7000, 0, 7000);
    TH1D* hist_double = new TH1D("hist_double", "hist_double", 7000, 0, 7000);

    for(Long64_t i = 0; i < tree_triple_coinc->GetEntries(); ++i)
    {
        tree_triple_coinc->GetEntry(i);
    
        for(ULong_t j = 0; j < ffTripleCoinc_ECond->size(); ++j)
        {
            hist_triple->Fill((a * ffTripleCoinc_ECond->at(j)) + b);
        }
    }

    for(Long64_t i = 0; i < tree_double_coinc->GetEntries(); ++i)
    {
        tree_double_coinc->GetEntry(i);
    
        for(ULong_t j = 0; j < ffGeBeta_ECond->size(); ++j)
        {
            hist_double->Fill((a * ffGeBeta_ECond->at(j)) + b);
        }

        cout << (100.*i/tree_double_coinc->GetEntries()) << endl;
    }

    TFile* output_file = new TFile(OutputFileName, "recreate");

    hist_triple->Write();
    hist_double->Write();

}