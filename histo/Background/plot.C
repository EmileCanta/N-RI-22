void plot()
{   
    UInt_t energy;
    UChar_t det;

    TFile *file = TFile::Open("../../raw_runs/TETRA23_RUN4.root");

    TTree *tree = (TTree*)file->Get("Narval_tree");

    TH1I *hist = new TH1I("hist", "hist", 12000, 0, 12000);

    tree->SetBranchAddress("Energy", &energy);
    tree->SetBranchAddress("Det_nbr", &det);

    for(int i = 0; i < tree->GetEntries(); i++)
    {   
        tree->GetEntry(i);

        if(det == 15)
        {
            hist->Fill(0.32819*energy+49.87);
        }
    }

    hist->Draw();
}