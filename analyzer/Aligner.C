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

void Aligner(const char* inputFileName, Double_t maxTime)
{   
    cout << "Removing missaligned cycles" << endl;

    std::vector<Double_t> *fBeta_Time = 0;
    std::vector<Double_t> *fBeta_Cycle = 0;

    std::vector<Double_t> *fTetra_Time = 0;
    std::vector<Double_t> *fTetra_Cycle = 0;

    std::vector<Double_t> fBeta_Bad_Cycle;
    std::vector<Double_t> fTetra_Bad_Cycle;
    
    std::vector<UInt_t> fAlignement;

    TFile* data_file = TFile::Open(inputFileName, "UPDATE");

    TTree* data_tree = (TTree*)data_file->Get("tsingle");

    TBranch* branch = data_tree->Branch("Alignement", &fAlignement);

    TH1D* AlignedBeta_Time_single = new TH1D("AlignedBeta_Time_single","AlignedBeta_Time_single",30000,0,30000);
    TH1D* AlignedTetra_Time_single = new TH1D("AlignedTetra_Time_single","AlignedTetra_Time_single",30000,0,30000);

    data_tree->SetBranchAddress("Beta_Time_single", &fBeta_Time);
    data_tree->SetBranchAddress("Beta_Cycle", &fBeta_Cycle);

    data_tree->SetBranchAddress("Tetra_Time_single", &fTetra_Time);
    data_tree->SetBranchAddress("Tetra_Cycle", &fTetra_Cycle);

    Long64_t entries = data_tree->GetEntries();

    for(Long64_t i = 0; i < entries; i++)
    {   
        data_tree->GetEntry(i);

        for(ULong_t j = 0; j < fBeta_Cycle->size(); ++j)
        {
            if(fBeta_Time->at(j) > maxTime)
            {
                fBeta_Bad_Cycle.push_back(fBeta_Cycle->at(j));
            }
        }

        for(ULong_t j = 0; j < fTetra_Cycle->size(); ++j)
        {
            if(fTetra_Time->at(j) > maxTime)
            {
                fTetra_Bad_Cycle.push_back(fTetra_Cycle->at(j));
            }
        }
    }

    for(Long64_t i = 0; i < entries; i++)
    {   
        fAlignement.clear();

        data_tree->GetEntry(i);

        for(ULong_t j = 0; j < fBeta_Cycle->size(); ++j)
        {
            if(std::find(fBeta_Bad_Cycle.begin(), fBeta_Bad_Cycle.end(), fBeta_Cycle->at(j)) != fBeta_Bad_Cycle.end())
            {
                fAlignement.push_back(0);
            }

            else
            {
                fAlignement.push_back(1);
            }
        }

        for(ULong_t j = 0; j < fTetra_Cycle->size(); ++j)
        {
            if(std::find(fTetra_Bad_Cycle.begin(), fTetra_Bad_Cycle.end(), fTetra_Cycle->at(j)) != fTetra_Bad_Cycle.end())
            {
                fAlignement.push_back(0);
            }

            else
            {
                fAlignement.push_back(1);
            }
        }

        branch->Fill();
    }

    cout << "Building Aligned time histograms" << endl;

    for(Long64_t i = 0; i < entries; i++)
    {   
        data_tree->GetEntry(i);

        for(ULong_t j = 0; j < fBeta_Time->size(); ++j)
        {
            if(fAlignement.at(j) == 1)
            {
                AlignedBeta_Time_single->Fill(fBeta_Time->at(j));
            }
        }

        for(ULong_t j = 0; j < fTetra_Time->size(); ++j)
        {
            if(fAlignement.at(j) == 1)
            {
                AlignedTetra_Time_single->Fill(fTetra_Time->at(j));
            }
        }
    }

    data_file->Write();
}