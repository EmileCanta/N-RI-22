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

void Aligner(const char* inputFileName, Double_t maxTime, Double_t Offset)
{   
    cout << "Removing missaligned cycles" << endl;

    std::vector<Double_t> *fBeta_Time = 0;
    std::vector<Double_t> *fBeta_Cycle = 0;

    std::vector<Double_t> *fTetra_Time = 0;
    std::vector<Double_t> *fTetra_Cycle = 0;

    std::vector<Double_t> fBeta_Bad_Cycle_single;
    std::vector<Double_t> fTetra_Bad_Cycle_single;
    
    std::vector<UInt_t> fAlignementSingle;

    TFile* data_file = TFile::Open(inputFileName, "UPDATE");

    TTree* data_tree_single = (TTree*)data_file->Get("tsingle");

    TBranch* branch_single = data_tree_single->Branch("AlignementSingle", &fAlignementSingle);

    TH1D* Beta_tSingle = new TH1D("AlignedBeta_tSingle","Beta_tSingle",30000,0,30000);
    TH1D* Tetra_tSingle = new TH1D("AlignedTetra_tSingle","Tetra_tSingle",30000,0,30000);

    Beta_tSingle->Reset();
    Tetra_tSingle->Reset();

    data_tree_single->SetBranchAddress("Beta_tSingle", &fBeta_Time);
    data_tree_single->SetBranchAddress("Beta_Cycle", &fBeta_Cycle);

    data_tree_single->SetBranchAddress("Tetra_tSingle", &fTetra_Time);
    data_tree_single->SetBranchAddress("Tetra_Cycle", &fTetra_Cycle);

    Long64_t entries_single = data_tree_single->GetEntries();

    for(Long64_t i = 0; i < entries_single; i++)
    {   
        data_tree_single->GetEntry(i);

        for(ULong_t j = 0; j < fBeta_Cycle->size(); ++j)
        {
            if(fBeta_Time->at(j) > maxTime)
            {
                fBeta_Bad_Cycle_single.push_back(fBeta_Cycle->at(j));
            }
        }

        for(ULong_t j = 0; j < fTetra_Cycle->size(); ++j)
        {
            if(fTetra_Time->at(j) > maxTime)
            {
                fTetra_Bad_Cycle_single.push_back(fTetra_Cycle->at(j));
            }
        }
    }

    for(Long64_t i = 0; i < entries_single; i++)
    {   
        fAlignementSingle.clear();

        data_tree_single->GetEntry(i);

        for(ULong_t j = 0; j < fBeta_Cycle->size(); ++j)
        {
            if(std::find(fBeta_Bad_Cycle_single.begin(), fBeta_Bad_Cycle_single.end(), fBeta_Cycle->at(j)) != fBeta_Bad_Cycle_single.end())
            {
                fAlignementSingle.push_back(0);
            }

            else
            {
                fAlignementSingle.push_back(1);
            }
        }

        for(ULong_t j = 0; j < fTetra_Cycle->size(); ++j)
        {
            if(std::find(fTetra_Bad_Cycle_single.begin(), fTetra_Bad_Cycle_single.end(), fTetra_Cycle->at(j)) != fTetra_Bad_Cycle_single.end())
            {
                fAlignementSingle.push_back(0);
            }

            else
            {
                fAlignementSingle.push_back(1);
            }
        }

        branch_single->Fill();
    }

    cout << "Building Aligned time histograms" << endl;

    for(Long64_t i = 0; i < entries_single; i++)
    {   
        data_tree_single->GetEntry(i);

        for(ULong_t j = 0; j < fBeta_Time->size(); ++j)
        {
            if(fAlignementSingle.at(j) == 1)
            {
                Beta_tSingle->Fill(fBeta_Time->at(j));
            }

            else
            {
                Beta_tSingle->Fill(fBeta_Time->at(j)-Offset);
            }
        }

        for(ULong_t j = 0; j < fTetra_Time->size(); ++j)
        {
            if(fAlignementSingle.at(j) == 1)
            {
                Tetra_tSingle->Fill(fTetra_Time->at(j));
            }

            else
            {
                Tetra_tSingle->Fill(fTetra_Time->at(j)-Offset);
            }
        }
    }

    data_file->Write();
}