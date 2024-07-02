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

void Aligner2n(const char* inputFileName, Double_t maxTime, Double_t Offset)
{   
    cout << "Removing missaligned cycles" << endl;

    std::vector<Double_t> *f2n_tCond = 0;
    std::vector<Double_t> *f2n_Cycle = 0;

    std::vector<Double_t> f2n_BadCycle;
    
    std::vector<UInt_t> fAlignement2n;

    TFile* data_file = TFile::Open(inputFileName, "UPDATE");

    TTree* data_tree_coinc = (TTree*)data_file->Get("tcoinc");

    TBranch* branch_coinc = data_tree_coinc->Branch("Alignement2n", &fAlignement2n);

    TH1D* Aligned2n_tCond = new TH1D("Aligned2n_tCond","Aligned2n_tCond",30000,0,30000);

    Aligned2n_tCond->Reset();

    data_tree_coinc->SetBranchAddress("SecondNeut_tCond", &f2n_tCond);
    data_tree_coinc->SetBranchAddress("SecondNeutronCycle", &f2n_Cycle);

    Long64_t entries_coinc = data_tree_coinc->GetEntries();

    for(Long64_t i = 0; i < entries_coinc; i++)
    {   
        data_tree_coinc->GetEntry(i);

        for(ULong_t j = 0; j < f2n_Cycle->size(); ++j)
        {
            if(f2n_tCond->at(j) > maxTime)
            {
                f2n_BadCycle.push_back(f2n_Cycle->at(j));
            }
        }
    }

    for(Long64_t i = 0; i < entries_coinc; i++)
    {   
        fAlignement2n.clear();

        data_tree_coinc->GetEntry(i);

        for(ULong_t j = 0; j < f2n_Cycle->size(); ++j)
        {
            if(std::find(f2n_BadCycle.begin(), f2n_BadCycle.end(), f2n_Cycle->at(j)) != f2n_BadCycle.end())
            {
                fAlignement2n.push_back(0);
            }

            else
            {
                fAlignement2n.push_back(1);
            }
        }

        branch_coinc->Fill();
    }

    cout << "Building Aligned time histograms" << endl;

    for(Long64_t i = 0; i < entries_coinc; i++)
    {   
        data_tree_coinc->GetEntry(i);

        for(ULong_t j = 0; j < f2n_tCond->size(); ++j)
        {
            if(fAlignement2n.at(j) == 1)
            {
                Aligned2n_tCond->Fill(f2n_tCond->at(j));
            }

            else
            {
                Aligned2n_tCond->Fill(f2n_tCond->at(j)-Offset);
            }
        }
    }

    data_file->Write();
}