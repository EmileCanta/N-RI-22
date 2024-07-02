void hist_substract()
{

    TFile *file;
    TTree *tsingle;

    file = TFile::Open("/Users/cantacuzene/data/n-ri-22/runs/sorted_runs/84Ga/All.root", "READ");

    tsingle = (TTree*)file->Get("tcoinc");

    TH1I* h1 = new TH1I("h1", "h1", 7000,0,7000);
    TH1I* h2 = new TH1I("h2", "h2", 7000,0,7000);

    std::vector<Double_t> *rings = 0;
    std::vector<Double_t> *time = 0;
    std::vector<Double_t> *diff = 0;

    tsingle->SetBranchAddress("TwoNeutronsGamma_ECond",&rings);
    tsingle->SetBranchAddress("TwoNeutronsGamma_tCond",&time);
    tsingle->SetBranchAddress("TwoNeutronsGamma_tDiff",&diff);

    for(Long64_t i = 0; i < tsingle->GetEntries(); ++i)
    {
        tsingle->GetEntry(i);
    
        for(ULong_t j = 0; j < rings->size(); ++j)
        {   
            if(diff->at(j) <= 40)
            {
                if(time->at(j) >= 1000 && time->at(j) <= 1900)
                {
                    h1->Fill(rings->at(j));
                }

                if(time->at(j) >= 2300 && time->at(j) <= 3200)
                {
                    h2->Fill(rings->at(j));
                }
            }
        }
    }

    h1->Add(h2,-1);

    h1->Draw();
}