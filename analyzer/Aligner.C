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

    std::vector<Double_t> *fOneNeutron_Time = 0;
    std::vector<Double_t> *fOneNeutron_Cycle = 0;

    std::vector<Double_t> fBeta_Bad_Cycle_single;
    std::vector<Double_t> fTetra_Bad_Cycle_single;

    std::vector<Double_t> fOneNeutron_Bad_Cycle;
    
    std::vector<UInt_t> fAlignementSingle;
    std::vector<UInt_t> fAlignementCoinc;

    TFile* data_file = TFile::Open(inputFileName, "UPDATE");

    TTree* data_tree_single = (TTree*)data_file->Get("tsingle");
    TTree* data_tree_coinc = (TTree*)data_file->Get("tcoinc");

    TBranch* branch_single = data_tree_single->Branch("AlignementSingle", &fAlignementSingle);
    TBranch* branch_coinc = data_tree_coinc->Branch("AlignementCoinc", &fAlignementCoinc);

    TH1D* AlignedBeta_Time_single = new TH1D("AlignedBeta_Time_single","AlignedBeta_Time_single",30000,0,30000);
    TH1D* AlignedTetra_Time_single = new TH1D("AlignedTetra_Time_single","AlignedTetra_Time_single",30000,0,30000);

    TH1D* UnalignedBeta_Time_single = new TH1D("UnalignedBeta_Time_single","UnalignedBeta_Time_single",30000,0,30000);
    TH1D* UnalignedTetra_Time_single = new TH1D("UnalignedTetra_Time_single","UnalignedTetra_Time_single",30000,0,30000);

    TH1D* AlignedOneNeutron_Time = new TH1D("AlignedOneNeutron_Time","AlignedOneNeutron_Time",30000,0,30000);

    AlignedBeta_Time_single->Reset();
    AlignedTetra_Time_single->Reset();
    UnalignedBeta_Time_single->Reset();
    UnalignedTetra_Time_single->Reset();
    AlignedOneNeutron_Time->Reset();

    data_tree_single->SetBranchAddress("Beta_Time_single", &fBeta_Time);
    data_tree_single->SetBranchAddress("Beta_Cycle", &fBeta_Cycle);

    data_tree_single->SetBranchAddress("Tetra_Time_single", &fTetra_Time);
    data_tree_single->SetBranchAddress("Tetra_Cycle", &fTetra_Cycle);

    data_tree_coinc->SetBranchAddress("OneNeutron_Time", &fOneNeutron_Time);
    data_tree_coinc->SetBranchAddress("OneNeutron_Cycle", &fOneNeutron_Cycle);

    Long64_t entries_single = data_tree_single->GetEntries();
    Long64_t entries_coinc = data_tree_coinc->GetEntries();

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

    for(Long64_t i = 0; i < entries_coinc; i++)
    {   
        data_tree_coinc->GetEntry(i);

        for(ULong_t j = 0; j < fOneNeutron_Cycle->size(); ++j)
        {
            if(fOneNeutron_Time->at(j) > maxTime)
            {
                fOneNeutron_Bad_Cycle.push_back(fOneNeutron_Cycle->at(j));
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

    for(Long64_t i = 0; i < entries_coinc; i++)
    {   
        fAlignementCoinc.clear();

        data_tree_coinc->GetEntry(i);

        for(ULong_t j = 0; j < fOneNeutron_Cycle->size(); ++j)
        {
            if(std::find(fOneNeutron_Bad_Cycle.begin(), fOneNeutron_Bad_Cycle.end(), fOneNeutron_Cycle->at(j)) != fOneNeutron_Bad_Cycle.end())
            {
                fAlignementCoinc.push_back(0);
            }

            else
            {
                fAlignementCoinc.push_back(1);
            }
        }

        branch_coinc->Fill();
    }

    cout << "Building Aligned time histograms" << endl;

    for(Long64_t i = 0; i < entries_single; i++)
    {   
        data_tree_single->GetEntry(i);

        for(ULong_t j = 0; j < fBeta_Time->size(); ++j)
        {
            if(fAlignementSingle.at(j) == 1)
            {
                AlignedBeta_Time_single->Fill(fBeta_Time->at(j));
            }

            else
            {
                AlignedBeta_Time_single->Fill(fBeta_Time->at(j)-Offset);
                UnalignedBeta_Time_single->Fill(fBeta_Time->at(j));
            }
        }

        for(ULong_t j = 0; j < fTetra_Time->size(); ++j)
        {
            if(fAlignementSingle.at(j) == 1)
            {
                AlignedTetra_Time_single->Fill(fTetra_Time->at(j));
            }

            else
            {
                AlignedTetra_Time_single->Fill(fTetra_Time->at(j)-Offset);
                UnalignedTetra_Time_single->Fill(fTetra_Time->at(j));
            }
        }
    }

    for(Long64_t i = 0; i < entries_coinc; i++)
    {   
        data_tree_coinc->GetEntry(i);

        for(ULong_t j = 0; j < fOneNeutron_Time->size(); ++j)
        {
            if(fAlignementCoinc.at(j) == 1)
            {
                AlignedOneNeutron_Time->Fill(fOneNeutron_Time->at(j));
            }

            else
            {
                AlignedOneNeutron_Time->Fill(fOneNeutron_Time->at(j)-Offset);
            }
        }
    }

    data_file->Write();
}