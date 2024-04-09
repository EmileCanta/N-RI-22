#include "Sorter.h"
 
Sorter::Sorter(const char* InputFileName, const char* OutputFileName, Double_t collection_time, Double_t acquisition_time, Double_t background_time)                        
{
    cout << "Starting Sorter" << endl;
    
    gROOT->ProcessLine("#include <vector>");

    coll_time = collection_time;
    acqu_time = acquisition_time;
    bgd_time = background_time;

    Load_rawdata(InputFileName);

    SetTreesAndBranches(OutputFileName);

    FillSingleBranches();

    FillNeutronGammaCoincBranches(200e6);

    FillBetaGammaCoincBranches(10e6);

    FillBetaNeutronCoincBranches(200e6);
    FillBetaNeutronBackwardCoincBranches(200e6);

    FillNeutronNeutronCoincBranches(200e6);

    FillOneNeutronBranches(200e6);

    output_file->Close();

    Histogrammer(OutputFileName);
    
    gROOT->ProcessLine(".q");
}

//******************************************************************************
//******************************************************************************

Sorter::~Sorter()
{
    ;
}

//******************************************************************************
//******************************************************************************

Double_t Sorter::Ge_alignement(UInt_t channel)
{
    Double_t channel_aligned = 0;
    Double_t ch = (channel + gRandom->Uniform(1.0) - 0.5);
    //channel_aligned = ch * 0.6264 + 96.469; /*Calibration runs production*/
    channel_aligned = ch * 0.17701043557 + 26.646713014; /*Calibration runs manip*/
    return channel_aligned;
}

//******************************************************************************
//******************************************************************************

UChar_t Sorter::GetStatus(Double_t time)
{
    UChar_t status = -10;
    if(time / 1e6 > bgd_time && time / 1e6 <= bgd_time + coll_time) status = 1; //Collection
    else if(time / 1e6 > acqu_time) status = 1; //Acquisition
    return status;
}

//******************************************************************************
//******************************************************************************

void Sorter::ResetVar()
{
    fTetra_Id = 0;
    fTetra_Time = 0.;

    fBeta_Id = 0;
    fBeta_Time = 0.;

    fGe_E = 0.;
    fGe_Id = 0;
    fGe_Time = 0.;

    fCoding = 0;
    fStatus = -10;
}

//******************************************************************************
//******************************************************************************

void Sorter::ClearVectors()
{
    fTetra_Time_single.clear();
    
    fBeta_Time_single.clear();

    fGe_E_single.clear();
    fGe_Time_single.clear();

    fTetra_Rings.clear();

    fGeBeta_E_coinc.clear();
    fGeTetra_E_coinc.clear();
        
    fGeBeta_Time_coinc.clear();
    fGeTetra_Time_coinc.clear();
    fBeta1n_Time_coinc.clear();
    fBeta2n_Time_coinc.clear();
    fBeta3n_Time_coinc.clear();
    fBeta1nBackward_Time_coinc.clear();
    fBeta2nBackward_Time_coinc.clear();
        
    fGeBeta_TimeDiff.clear();
    fGeTetra_TimeDiff.clear();

    fBeta1n_TimeDiff.clear();
    fBeta2n_TimeDiffFirst.clear();
    fBeta2n_TimeDiffSecond.clear();
    fBeta1nBackward_TimeDiff.clear();
    fBeta2nBackward_TimeDiffFirst.clear();
    fBeta2nBackward_TimeDiffSecond.clear();

    fGeTetra_Index.clear();
    fGeBeta_Index.clear();

    fStoring1n_Time.clear();
    fStoring1n_TimeDiff.clear();
    fStoring2n_Time.clear();
    fStoring2n_TimeDiff.clear();
    fStoring3n_Time.clear();
    fStoring3n_TimeDiff.clear();

    f1n_Time.clear();
    f1n_Cycle.clear();

    f2n_TimeDiff.clear();
    f2n_Time_coinc.clear();

    fStoringFirstNeutronCellGroup.clear();
    fStoringSecondNeutronCellGroup.clear();

    fFirstNeutronCellGroup.clear();
    fSecondNeutronCellGroup.clear();

    fTetra_Cycle.clear();
    fBeta_Cycle.clear();
    fGe_Cycle.clear();
}

//******************************************************************************
//******************************************************************************

void Sorter::SetVar(UChar_t det)
{   
    if(det >= 1 && det <= 12)
    {
        fTetra_Id = det;
        fTetra_Time = raw_time;
        fCoding = raw_coding;
        fStatus = GetStatus(raw_time);
    }

    else if(det == 14)
    { 
        fBeta_Id = det;
        fBeta_Time = raw_time;
        fCoding = raw_coding;
        fStatus = GetStatus(raw_time);
    }
        
    else if(det == 15)
    {   
        fGe_E = raw_energy;
        fGe_Id = det;
        fGe_Time = raw_time;
        fCoding = raw_coding;
        fStatus = GetStatus(raw_time);
    }
}

//******************************************************************************
//******************************************************************************

void Sorter::Load_rawdata (const char* InputFileName)
{
	cout << "Loading Narval Tree" << endl;
	
	data_file = TFile::Open(InputFileName);

    raw_data_tree = (TTree*)data_file->Get("Narval_tree");
	
	raw_data_tree->SetBranchAddress("Energy", &raw_energy);
    raw_data_tree->SetBranchAddress("Time", &raw_time);
    raw_data_tree->SetBranchAddress("Coding", &raw_coding);
    raw_data_tree->SetBranchAddress("Det_nbr", &raw_det_nbr);

    fEntries = raw_data_tree->GetEntries();
    
    std::cout << "Entries to be sorted : " << fEntries << std::endl;
}

//******************************************************************************
//******************************************************************************

void Sorter::SetTreesAndBranches(const char* OutputFileName)
{
    std::cout << "Initialization of output TTree" << std::endl;

    std::cout << "Creating branches" << std::endl;
    
    output_file = new TFile(OutputFileName, "recreate");
    
    output_tree_single = new TTree("tsingle", "Tree with singles");
    output_tree_coinc = new TTree("tcoinc", "Tree with coincidences");
      
    //Single branches
    
    output_tree_single->Branch("Tetra_Time_single", &fTetra_Time_single);
    output_tree_single->Branch("Tetra_Cycle", &fTetra_Cycle);

    output_tree_single->Branch("Ge_E_single", &fGe_E_single);
    output_tree_single->Branch("Ge_Time_single", &fGe_Time_single);
    output_tree_single->Branch("Beta_Cycle", &fBeta_Cycle);

    output_tree_single->Branch("Beta_Time_single", &fBeta_Time_single);
    output_tree_single->Branch("Ge_Cycle", &fGe_Cycle);

    output_tree_single->Branch("Tetra_Rings", &fTetra_Rings);
    
    //Coinc branches

    //Energy
    output_tree_coinc->Branch("GeBeta_E_coinc", &fGeBeta_E_coinc);
    output_tree_coinc->Branch("GeTetra_E_coinc", &fGeTetra_E_coinc);
    
    //Time
    output_tree_coinc->Branch("GeBeta_Time_coinc", &fGeBeta_Time_coinc);
    output_tree_coinc->Branch("GeTetra_Time_coinc", &fGeTetra_Time_coinc);
    output_tree_coinc->Branch("Beta1n_Time_coinc", &fBeta1n_Time_coinc);
    output_tree_coinc->Branch("Beta2n_Time_coinc", &fBeta2n_Time_coinc);
    output_tree_coinc->Branch("Beta3n_Time_coinc", &fBeta3n_Time_coinc);
    output_tree_coinc->Branch("TwoNeutrons_Time_coinc", &f2n_Time_coinc);
    output_tree_coinc->Branch("OneNeutron_Time", &f1n_Time);
    output_tree_coinc->Branch("OneNeutron_Cycle", &f1n_Cycle);
    output_tree_coinc->Branch("Beta1nBackward_Time_coinc", &fBeta1nBackward_Time_coinc);
    output_tree_coinc->Branch("Beta2nBackward_Time_coinc", &fBeta2nBackward_Time_coinc);
    
    //Time Diff
    output_tree_coinc->Branch("GeBeta_TimeDiff", &fGeBeta_TimeDiff);
    output_tree_coinc->Branch("GeTetra_TimeDiff", &fGeTetra_TimeDiff);
    output_tree_coinc->Branch("Beta1n_TimeDiff", &fBeta1n_TimeDiff);
    output_tree_coinc->Branch("Beta2n_TimeDiffFirst", &fBeta2n_TimeDiffFirst);
    output_tree_coinc->Branch("Beta2n_TimeDiffSecond", &fBeta2n_TimeDiffSecond);
    output_tree_coinc->Branch("Beta1nBackward_TimeDiff", &fBeta1nBackward_TimeDiff);
    output_tree_coinc->Branch("Beta2nBackward_TimeDiffFirst", &fBeta2nBackward_TimeDiffFirst);
    output_tree_coinc->Branch("Beta2nBackward_TimeDiffSecond", &fBeta2nBackward_TimeDiffSecond);
    output_tree_coinc->Branch("TwoNeutrons_TimeDiff", &f2n_TimeDiff);

    //Raw index
    output_tree_coinc->Branch("GeTetra_Index", &fGeTetra_Index);
    output_tree_coinc->Branch("GeBeta_Index", &fGeBeta_Index);

    //Cell groups
    output_tree_coinc->Branch("FirstNeutronCellGroup", &fFirstNeutronCellGroup);
    output_tree_coinc->Branch("SecondNeutronCellGroup", &fSecondNeutronCellGroup);

}

//******************************************************************************
//******************************************************************************

void Sorter::FillSingleBranches()
{   
    std::cout << "Sorting single data" << std::endl;

    Double_t last_time = 0.;
    Double_t fCycle = 0.;

    ResetVar();
    
    for(fEntry = 0; fEntry < fEntries; fEntry++)
	{
        raw_data_tree->GetEntry(fEntry);

        if(fEntry == 0)
        {
            last_time = raw_time;
            fCycle = 0;
	    }

        if(raw_time - last_time < 0) fCycle++;
        
        last_time = raw_time;

        ClearVectors();

        if(raw_det_nbr >= 1 && raw_det_nbr <= 12)
        {
            fTetra_Time_single.push_back(raw_time / 1e9); //millisecond
            fTetra_Cycle.push_back(fCycle);
            
            for(int i = 1; i <= 12; i++)
            {
                if(raw_det_nbr == i)
                {
                    fTetra_Rings.push_back(i);
                }
            }
        }

        else if(raw_det_nbr == 14)
        {
            fBeta_Time_single.push_back(raw_time / 1e9); //millisecond
            fBeta_Cycle.push_back(fCycle);
        }

        else if(raw_det_nbr == 15)
        {
            fGe_E_single.push_back(Ge_alignement(raw_energy));
            fGe_Time_single.push_back(raw_time / 1e9); //millisecond
            fGe_Cycle.push_back(fCycle);
        }
        
        output_tree_single->Fill();
        
    }
}

//******************************************************************************
//******************************************************************************

void Sorter::FillNeutronGammaCoincBranches(Double_t coincwindow)
{
    std::cout << "Sorting neutron-gamma coinc data" << std::endl;
    
    lastevent = 0;

    ResetVar();

	for(fEntry = fEntries; fEntry > 0; fEntry--)
	{
        raw_data_tree->GetEntry(fEntry);

        ClearVectors();
        
        SetVar(raw_det_nbr);
        
        eventafter = 0;

        if(fTetra_Id >= 1 && fTetra_Id <= 12)
        {   
            if(fEntry > lastevent) continue;

            while(TMath::Abs(raw_time - fTetra_Time) < coincwindow)
            {
	      		if(raw_det_nbr >= 1 && raw_det_nbr <= 12) break;

                raw_data_tree->GetEntry(fEntry - eventafter++);

                if(fEntry - eventafter < 0) break;
                if(fCoding != raw_coding) break;
                if(fStatus != GetStatus(raw_time)) break;
                
                if(TMath::Abs(raw_time - fTetra_Time) > coincwindow) continue;
                if(TMath::Abs(raw_time - fTetra_Time) < coincwindow)
                {   
                    if(raw_det_nbr == 15)
                    {   
                        fGeTetra_TimeDiff.push_back(TMath::Abs(raw_time - fTetra_Time) / 1e6); //microsecondes
				        fGeTetra_E_coinc.push_back(Ge_alignement(raw_energy));
				        fGeTetra_Time_coinc.push_back(raw_time / 1e9); //milisecondes;
                        fGeTetra_Index.push_back(fEntry - eventafter);
                    }
                 }
            }
	    }
	    
	    lastevent = fEntry - eventafter;
	    
	    output_tree_coinc->Fill();
    }
}

//******************************************************************************
//******************************************************************************

void Sorter::FillBetaGammaCoincBranches(Double_t coincwindow)
{
    std::cout << "Sorting beta-gamma coinc data" << std::endl;
    
    Int_t lastevent = 0;

    ResetVar();

	for(fEntry = 0; fEntry < fEntries; fEntry++)
	{
        raw_data_tree->GetEntry(fEntry);

        ClearVectors();
        
        SetVar(raw_det_nbr);
        
        eventafter = 0;

        if(fBeta_Id == 14)
        {   
            if(fEntry < lastevent) continue;
            
            while(TMath::Abs(raw_time - fBeta_Time) < coincwindow)
            {
	      		if(raw_det_nbr == 14) break;

                raw_data_tree->GetEntry(fEntry + eventafter++);

                if(fEntry + eventafter > fEntries) break;
                if(fCoding != raw_coding) break;
                if(fStatus != GetStatus(raw_time)) break;
                
                if(TMath::Abs(raw_time - fBeta_Time) > coincwindow) continue;
                if(TMath::Abs(raw_time - fBeta_Time) < coincwindow)
                {   
                    if(raw_det_nbr == 15)
                    {
                        fGeBeta_TimeDiff.push_back(TMath::Abs(raw_time - fBeta_Time) / 1e6); //microsecondes
                        fGeBeta_E_coinc.push_back(Ge_alignement(raw_energy));
                        fGeBeta_Time_coinc.push_back(raw_time / 1e9); //milisecondes;
                        fGeBeta_Index.push_back(fEntry + eventafter);
                    }
                 }
            }
	    }
	    
	    lastevent = fEntry + eventafter;
	    
	    output_tree_coinc->Fill();
    }
}

//******************************************************************************
//******************************************************************************

void Sorter::FillBetaNeutronCoincBranches(Double_t coincwindow)
{
    std::cout << "Sorting beta-neutron coinc data" << std::endl;
    
    lastevent = 0;

    ResetVar();

	for(fEntry = 0; fEntry < fEntries; fEntry++)
	{
        raw_data_tree->GetEntry(fEntry);

        ClearVectors();
        
        SetVar(raw_det_nbr);
        
        eventafter = 0;

        int neutcount = 0;

        if(fBeta_Id == 14)
        {   
            if(fEntry < lastevent) continue;
            
            while(TMath::Abs(raw_time - fBeta_Time) < coincwindow)
            {
	      		if(raw_det_nbr == 14) break;

                raw_data_tree->GetEntry(fEntry + eventafter++);

                if(fEntry + eventafter > fEntries) break;
                if(fCoding != raw_coding) break;
                if(fStatus != GetStatus(raw_time)) break;
                
                if(TMath::Abs(raw_time - fBeta_Time) > coincwindow) continue;
                if(TMath::Abs(raw_time - fBeta_Time) < coincwindow)
                {   
                    if(raw_det_nbr >= 1 && raw_det_nbr <= 12)
                    {
                        neutcount++;

                        if(neutcount == 1)
                        {
                            fStoring1n_Time.push_back(fBeta_Time / 1e9);
                            fStoring1n_TimeDiff.push_back(TMath::Abs(raw_time - fBeta_Time) / 1e6);
                        }

                        if(neutcount == 2)
                        {
                            fStoring2n_Time.push_back(fBeta_Time / 1e9);
                            fStoring2n_TimeDiff.push_back(TMath::Abs(raw_time - fBeta_Time) / 1e6);
                        }

                        if(neutcount == 3)
                        {
                            fStoring3n_Time.push_back(fBeta_Time / 1e9);
                        }
                    }
                 }
            }

            if(neutcount == 1)
            {
                fBeta1n_Time_coinc.push_back(fStoring1n_Time.front());
                fBeta1n_TimeDiff.push_back(fStoring1n_TimeDiff.front());
            }

            if(neutcount == 2)
            {
                fBeta2n_Time_coinc.push_back(fStoring2n_Time.front());
                fBeta2n_TimeDiffFirst.push_back(fStoring1n_TimeDiff.front());
                fBeta2n_TimeDiffSecond.push_back(fStoring2n_TimeDiff.front());
            }

            if(neutcount == 3)
            {
                fBeta3n_Time_coinc.push_back(fStoring3n_Time.front());
            }
	    }

	    lastevent = fEntry + eventafter;
	    
	    output_tree_coinc->Fill();
    }
}

//******************************************************************************
//******************************************************************************

void Sorter::FillBetaNeutronBackwardCoincBranches(Double_t coincwindow)
{
    std::cout << "Sorting beta-neutron backward coinc data" << std::endl;
    
    lastevent = 0;

    ResetVar();

	for(fEntry = fEntries; fEntry > 0; fEntry--)
	{
        raw_data_tree->GetEntry(fEntry);

        ClearVectors();
        
        SetVar(raw_det_nbr);
        
        eventafter = 0;

        int neutcount = 0;

        if(fBeta_Id == 14)
        {   
            if(fEntry > lastevent) continue;
            
            while(TMath::Abs(raw_time - fBeta_Time) < coincwindow)
            {
	      		if(raw_det_nbr == 14) break;

                raw_data_tree->GetEntry(fEntry - eventafter++);

                if(fEntry - eventafter < 0) break;
                if(fCoding != raw_coding) break;
                if(fStatus != GetStatus(raw_time)) break;
                
                if(TMath::Abs(raw_time - fBeta_Time) > coincwindow) continue;
                if(TMath::Abs(raw_time - fBeta_Time) < coincwindow)
                {   
                    if(raw_det_nbr >= 1 && raw_det_nbr <= 12)
                    {
                        neutcount++;

                        if(neutcount == 1)
                        {
                            fStoring1n_Time.push_back(fBeta_Time / 1e9);
                            fStoring1n_TimeDiff.push_back((raw_time - fBeta_Time) / 1e6);
                        }

                        if(neutcount == 2)
                        {
                            fStoring2n_Time.push_back(fBeta_Time / 1e9);
                            fStoring2n_TimeDiff.push_back((raw_time - fBeta_Time) / 1e6);
                        }

                        if(neutcount == 3)
                        {
                            fStoring3n_Time.push_back(fBeta_Time / 1e9);
                        }
                    }
                 }
            }

            if(neutcount == 1)
            {
                fBeta1nBackward_Time_coinc.push_back(fStoring1n_Time.front());
                fBeta1nBackward_TimeDiff.push_back(fStoring1n_TimeDiff.front());
            }

            if(neutcount == 2)
            {
                fBeta2nBackward_Time_coinc.push_back(fStoring2n_Time.front());
                fBeta2nBackward_TimeDiffFirst.push_back(fStoring1n_TimeDiff.front());
                fBeta2nBackward_TimeDiffSecond.push_back(fStoring2n_TimeDiff.front());
            }
	    }

	    lastevent = fEntry - eventafter;
	    
	    output_tree_coinc->Fill();
    }
}

//******************************************************************************
//******************************************************************************

void Sorter::FillNeutronNeutronCoincBranches(Double_t coincwindow)
{
    std::cout << "Sorting neutron-neutron coinc data" << std::endl;

    lastevent = 0;

    ResetVar();

	for(fEntry = 0; fEntry < fEntries; fEntry++)
	{
        raw_data_tree->GetEntry(fEntry);

        ClearVectors();
        
        SetVar(raw_det_nbr);
        
        eventafter = 0;

        int neutcount = 0;

        if(fTetra_Id >= 1 && fTetra_Id <= 12)
        {  
            if(fEntry < lastevent) continue;

            fStoringFirstNeutronCellGroup.push_back(fTetra_Id);
            
            while(TMath::Abs(raw_time - fTetra_Time) < coincwindow)
            {
                raw_data_tree->GetEntry(fEntry + eventafter++);

                if(fEntry + eventafter > fEntries) break;
                if(fCoding != raw_coding) break;
                if(fStatus != GetStatus(raw_time)) break;
                
                if(TMath::Abs(raw_time - fTetra_Time) > coincwindow) continue;
                if(TMath::Abs(raw_time - fTetra_Time) < coincwindow)
                {
                    if(raw_det_nbr >= 1 && raw_det_nbr <= 12)
                    {
                        if(TMath::Abs(raw_time - fTetra_Time != 0)) 
                        {
                            neutcount++;

                            if(neutcount == 1)
                            {
                                fStoring2n_Time.push_back(raw_time / 1e9);
                                fStoring2n_TimeDiff.push_back(TMath::Abs(raw_time - fTetra_Time) / 1e6);
                                fStoringSecondNeutronCellGroup.push_back(raw_det_nbr);
                            }
                        }
                    }
                }
	        }

            if(neutcount == 1)
            {
                f2n_TimeDiff.push_back(fStoring2n_TimeDiff.front());
                f2n_Time_coinc.push_back(fStoring2n_Time.front());

                fFirstNeutronCellGroup.push_back(fStoringFirstNeutronCellGroup.front());
                fSecondNeutronCellGroup.push_back(fStoringSecondNeutronCellGroup.front());
            }

	    lastevent = fEntry + eventafter;
	    
	    output_tree_coinc->Fill();
        }
    }
}

//******************************************************************************
//******************************************************************************

void Sorter::FillOneNeutronBranches(Double_t coincwindow)
{
    std::cout << "Sorting one-neutron coinc data" << std::endl;

    Double_t last_time = 0.;
    Double_t fCycle = 0.;
    
    lastevent = 0;

    ResetVar();

	for(fEntry = 0; fEntry < fEntries; fEntry++)
	{
        raw_data_tree->GetEntry(fEntry);

        if(fEntry == 0)
        {
            last_time = raw_time;
            fCycle = 0;
        }

        if(raw_time - last_time < 0) fCycle++;
        
        last_time = raw_time;

        ClearVectors();
        
        SetVar(raw_det_nbr);
        
        eventafter = 0;

        if(fTetra_Id >= 1 && fTetra_Id <= 12)
        {  
            if(fEntry < lastevent) continue;    
            
            while(TMath::Abs(raw_time - fTetra_Time) < coincwindow)
            {
                if(raw_det_nbr >= 1 && raw_det_nbr <= 12 && TMath::Abs(raw_time - fTetra_Time != 0))
                {   
                    break;
                }

                raw_data_tree->GetEntry(fEntry + eventafter++);

                if(fEntry + eventafter > fEntries) break;
                if(fCoding != raw_coding) break;
                if(fStatus != GetStatus(raw_time)) break;
                
                if(TMath::Abs(raw_time - fTetra_Time) > coincwindow) 
                {
                    f1n_Time.push_back(fTetra_Time / 1e9);
                    f1n_Cycle.push_back(fCycle);
                }
	        }

	    lastevent = fEntry + eventafter;
	    
	    output_tree_coinc->Fill();
        }
    }
}

//******************************************************************************
//******************************************************************************

void Sorter::Histogrammer(const char* OutputFileName)
{
    std::cout << "Building histograms" << std::endl;

    TFile *file;
    TTree *tsingle;
    TTree *tcoinc;

    file = TFile::Open(OutputFileName, "UPDATE");

    tsingle = (TTree*)file->Get("tsingle");
    tcoinc = (TTree*)file->Get("tcoinc");

    //Single
    TH1D* Ge_E_single = new TH1D("Ge_E_single", "Ge_E_single", 7000,0,7000);

    TH1D* Ge_Time_single = new TH1D("Ge_Time_single", "Ge_Time_single", 30000,0,30000);
    TH1D* Beta_Time_single = new TH1D("Beta_Time_single", "Beta_Time_single", 30000,0,30000);
    TH1D* Tetra_Time_single = new TH1D("Tetra_Time_single", "Tetra_Time_single", 30000,0,30000);

    TH1D* Tetra_Cycle = new TH1D("Tetra_Cycle", "Tetra_Cycle", 1000,0,1000);
    TH1D* Beta_Cycle = new TH1D("Beta_Cycle", "Beta_Cycle", 1000,0,1000);
    TH1D* Ge_Cycle = new TH1D("Ge_Cycle", "Ge_Cycle", 1000,0,1000);

    //Coinc
    TH1D* GeBeta_TimeDiff = new TH1D("GeBeta_TimeDiff","GeBeta_TimeDiff", 5000, 0, 10);
    TH1D* GeTetra_TimeDiff = new TH1D("GeTetra_TimeDiff", "GeTetra_TimeDiff", 1000, 0, 200);
    TH1D* Beta1n_TimeDiff = new TH1D("Beta1n_TimeDiff", "Beta1n_TimeDiff", 20000, -2000, 2000);
    TH1D* Beta2n_TimeDiffFirst = new TH1D("Beta2n_TimeDiffFirst", "Beta2n_TimeDiffFirst", 20000, -2000, 2000);
    TH1D* Beta2n_TimeDiffSecond = new TH1D("Beta2n_TimeDiffSecond", "Beta2n_TimeDiffSecond", 20000, -2000, 2000);
    TH1D* Beta1nBackward_TimeDiff = new TH1D("Beta1nBackward_TimeDiff", "Beta1nBackward_TimeDiff", 20000, -2000, 2000);
    TH1D* Beta2nBackward_TimeDiffFirst = new TH1D("Beta2nBackward_TimeDiffFirst", "Beta2nBackward_TimeDiffFirst", 20000, -2000, 2000);
    TH1D* Beta2nBackward_TimeDiffSecond = new TH1D("Beta2nBackward_TimeDiffSecond", "Beta2nBackward_TimeDiffSecond", 20000, -2000, 2000);
    TH1D* TwoNeutrons_TimeDiff = new TH1D("TwoNeutrons_TimeDiff", "TwoNeutrons_TimeDiff", 10000, 0, 2000);

    TH1D* GeBeta_Time_coinc = new TH1D("GeBeta_Time_coinc","GeBeta_Time_coinc", 30000, 0, 30000);
    TH1D* GeTetra_Time_coinc = new TH1D("GeTetra_Time_coinc","GeTetra_Time_coinc", 30000, 0, 30000);
    TH1D* Beta1n_Time_coinc = new TH1D("Beta1n_Time_coinc","Beta1n_Time_coinc", 30000, 0, 30000);
    TH1D* Beta2n_Time_coinc = new TH1D("Beta2n_Time_coinc","Beta2n_Time_coinc", 30000, 0, 30000);
    TH1D* Beta3n_Time_coinc = new TH1D("Beta3n_Time_coinc","Beta3n_Time_coinc", 30000, 0, 30000);
    TH1D* TwoNeutrons_Time_coinc = new TH1D("TwoNeutrons_Time_coinc", "TwoNeutrons_Time_coinc", 30000, 0, 30000);
    TH1D* OneNeutron_Time = new TH1D("OneNeutron_Time", "OneNeutron_Time", 30000, 0, 30000);
    TH1D* OneNeutron_Cycle = new TH1D("OneNeutron_Cycle", "OneNeutron_Cycle", 1000, 0, 1000);
    TH1D* Beta1nBackward_Time_coinc = new TH1D("Beta1nBackward_Time_coinc","Beta1nBackward_Time_coinc", 30000, 0, 30000);
    TH1D* Beta2nBackward_Time_coinc = new TH1D("Beta2nBackward_Time_coinc","Beta2nBackward_Time_coinc", 30000, 0, 30000);

    TH1D* GeBeta_E_coinc = new TH1D("GeBeta_E_coinc", "GeBeta_E_coinc", 7000, 0, 7000);
    TH1D* GeTetra_E_coinc = new TH1D("GeTetra_E_coinc", "GeTetra_E_coinc", 7000, 0, 7000);

    TH1D* FirstNeutronCellGroup = new TH1D("FirstNeutronCellGroup", "FirstNeutronCellGroup", 28, 0, 14);
    TH1D* SecondNeutronCellGroup = new TH1D("SecondNeutronCellGroup", "SecondNeutronCellGroup", 28, 0, 14);

    TH1D* Tetra_Rings = new TH1D("Tetra_Rings", "Tetra_Rings", 28, 0, 14);

    //Bidim
    TH2D* ESvsBTD = new TH2D("ESvsBTD", "ESvsBTD", 5000, 0, 10, 7000, 0, 7000);
    TH2D* ESvsTTD = new TH2D("ESvsTTD", "ESvsTTD", 1000, 0, 200, 7000, 0, 7000);

    TH2D* CellGroups = new TH2D("CellGroups", "CellGroups", 28, 0, 14, 28, 0, 14);

    std::vector<Double_t> *ffGe_E_single = 0;

    std::vector<Double_t> *ffGe_Time_single = 0;
    std::vector<Double_t> *ffBeta_Time_single = 0;
    std::vector<Double_t> *ffTetra_Time_single = 0;

    std::vector<Double_t> *ffTetra_Cycle = 0;
    std::vector<Double_t> *ffBeta_Cycle = 0;
    std::vector<Double_t> *ffGe_Cycle = 0;

    std::vector<Double_t> *ffGeBeta_TimeDiff = 0;
    std::vector<Double_t> *ffGeTetra_TimeDiff = 0;
    std::vector<Double_t> *ffBeta1n_TimeDiff = 0;
    std::vector<Double_t> *ffBeta2n_TimeDiffFirst = 0;
    std::vector<Double_t> *ffBeta2n_TimeDiffSecond = 0;
    std::vector<Double_t> *ffBeta1nBackward_TimeDiff = 0;
    std::vector<Double_t> *ffBeta2nBackward_TimeDiffFirst = 0;
    std::vector<Double_t> *ffBeta2nBackward_TimeDiffSecond = 0;
    std::vector<Double_t> *ff2n_TimeDiff = 0;

    std::vector<Double_t> *ffGeBeta_Time_coinc = 0;
    std::vector<Double_t> *ffGeTetra_Time_coinc = 0;
    std::vector<Double_t> *ffBeta1n_Time_coinc = 0;
    std::vector<Double_t> *ffBeta2n_Time_coinc = 0;
    std::vector<Double_t> *ffBeta3n_Time_coinc = 0;
    std::vector<Double_t> *ff2n_Time_coinc = 0;
    std::vector<Double_t> *ff1n_Time = 0;
    std::vector<Double_t> *ff1n_Cycle = 0;
    std::vector<Double_t> *ffBeta1nBackward_Time_coinc = 0;
    std::vector<Double_t> *ffBeta2nBackward_Time_coinc = 0;

    std::vector<Double_t> *ffGeBeta_E_coinc = 0;
    std::vector<Double_t> *ffGeTetra_E_coinc = 0;

    std::vector<Double_t> *ffFirstNeutronCellGroup = 0;
    std::vector<Double_t> *ffSecondNeutronCellGroup = 0;

    std::vector<Double_t> *ffTetra_Rings = 0;

    tsingle->SetBranchAddress("Ge_E_single",&ffGe_E_single);

    tsingle->SetBranchAddress("Ge_Time_single",&ffGe_Time_single);
    tsingle->SetBranchAddress("Beta_Time_single",&ffBeta_Time_single);
    tsingle->SetBranchAddress("Tetra_Time_single",&ffTetra_Time_single);

    tsingle->SetBranchAddress("Tetra_Cycle",&ffTetra_Cycle);
    tsingle->SetBranchAddress("Beta_Cycle",&ffBeta_Cycle);
    tsingle->SetBranchAddress("Ge_Cycle",&ffGe_Cycle);

    tsingle->SetBranchAddress("Tetra_Rings",&ffTetra_Rings);

    tcoinc->SetBranchAddress("GeBeta_TimeDiff",&ffGeBeta_TimeDiff);
    tcoinc->SetBranchAddress("GeTetra_TimeDiff",&ffGeTetra_TimeDiff);
    tcoinc->SetBranchAddress("Beta1n_TimeDiff",&ffBeta1n_TimeDiff);
    tcoinc->SetBranchAddress("Beta2n_TimeDiffFirst",&ffBeta2n_TimeDiffFirst);
    tcoinc->SetBranchAddress("Beta2n_TimeDiffSecond",&ffBeta2n_TimeDiffSecond);
    tcoinc->SetBranchAddress("Beta1nBackward_TimeDiff",&ffBeta1nBackward_TimeDiff);
    tcoinc->SetBranchAddress("Beta2nBackward_TimeDiffFirst",&ffBeta2nBackward_TimeDiffFirst);
    tcoinc->SetBranchAddress("Beta2nBackward_TimeDiffSecond",&ffBeta2nBackward_TimeDiffSecond);
    tcoinc->SetBranchAddress("TwoNeutrons_TimeDiff",&ff2n_TimeDiff);

    tcoinc->SetBranchAddress("GeBeta_Time_coinc",&ffGeBeta_Time_coinc);
    tcoinc->SetBranchAddress("GeTetra_Time_coinc",&ffGeTetra_Time_coinc);
    tcoinc->SetBranchAddress("Beta1n_Time_coinc",&ffBeta1n_Time_coinc);
    tcoinc->SetBranchAddress("Beta2n_Time_coinc",&ffBeta2n_Time_coinc);
    tcoinc->SetBranchAddress("Beta3n_Time_coinc",&ffBeta3n_Time_coinc);
    tcoinc->SetBranchAddress("TwoNeutrons_Time_coinc",&ff2n_Time_coinc);
    tcoinc->SetBranchAddress("OneNeutron_Time",&ff1n_Time);
    tcoinc->SetBranchAddress("OneNeutron_Cycle",&ff1n_Cycle);
    tcoinc->SetBranchAddress("Beta1nBackward_Time_coinc",&ffBeta1nBackward_Time_coinc);
    tcoinc->SetBranchAddress("Beta2nBackward_Time_coinc",&ffBeta2nBackward_Time_coinc);

    tcoinc->SetBranchAddress("GeBeta_E_coinc",&ffGeBeta_E_coinc);
    tcoinc->SetBranchAddress("GeTetra_E_coinc",&ffGeTetra_E_coinc);

    tcoinc->SetBranchAddress("FirstNeutronCellGroup",&ffFirstNeutronCellGroup);
    tcoinc->SetBranchAddress("SecondNeutronCellGroup",&ffSecondNeutronCellGroup);

    for(Long64_t i = 0; i < tsingle->GetEntries(); ++i)
    {
        tsingle->GetEntry(i);
    
        for(ULong_t j = 0; j < ffGe_E_single->size(); ++j)
        {
            Ge_E_single->Fill(ffGe_E_single->at(j));
            Ge_Time_single->Fill(ffGe_Time_single->at(j));
        }
        
        for(ULong_t j = 0; j < ffGe_Cycle->size(); ++j)
        {
            Ge_Cycle->Fill(ffGe_Cycle->at(j));
        }

        for(ULong_t j = 0; j < ffBeta_Time_single->size(); ++j)
        {
            Beta_Time_single->Fill(ffBeta_Time_single->at(j));
        }

        for(ULong_t j = 0; j < ffTetra_Cycle->size(); ++j)
        {
            Tetra_Cycle->Fill(ffTetra_Cycle->at(j));
        }

        for(ULong_t j = 0; j < ffTetra_Time_single->size(); ++j)
        {
            Tetra_Time_single->Fill(ffTetra_Time_single->at(j));
            Tetra_Rings->Fill(ffTetra_Rings->at(j));
        }

        for(ULong_t j = 0; j < ffBeta_Cycle->size(); ++j)
        {
            Beta_Cycle->Fill(ffBeta_Cycle->at(j));
        }
    }

    for(Long64_t i = 0; i < tcoinc->GetEntries(); ++i)
    {
        tcoinc->GetEntry(i);
    
        for(ULong_t j = 0; j < ffGeBeta_TimeDiff->size(); ++j)
        {
            GeBeta_TimeDiff->Fill(ffGeBeta_TimeDiff->at(j));
            GeBeta_Time_coinc->Fill(ffGeBeta_Time_coinc->at(j));
            GeBeta_E_coinc->Fill(ffGeBeta_E_coinc->at(j));
            ESvsBTD->Fill(ffGeBeta_TimeDiff->at(j), ffGeBeta_E_coinc->at(j));
        }

        for(ULong_t j = 0; j < ffGeTetra_TimeDiff->size(); ++j)
        {
            GeTetra_TimeDiff->Fill(ffGeTetra_TimeDiff->at(j));
            GeTetra_Time_coinc->Fill(ffGeTetra_Time_coinc->at(j));
            GeTetra_E_coinc->Fill(ffGeTetra_E_coinc->at(j));
            ESvsTTD->Fill(ffGeTetra_TimeDiff->at(j), ffGeTetra_E_coinc->at(j));
        }

        for(ULong_t j = 0; j < ffBeta1n_TimeDiff->size(); ++j)
        {
            Beta1n_TimeDiff->Fill(ffBeta1n_TimeDiff->at(j));
            Beta1n_Time_coinc->Fill(ffBeta1n_Time_coinc->at(j));
        }

        for(ULong_t j = 0; j < ffBeta2n_Time_coinc->size(); ++j)
        {
            Beta2n_Time_coinc->Fill(ffBeta2n_Time_coinc->at(j));
            Beta2n_TimeDiffFirst->Fill(ffBeta2n_TimeDiffFirst->at(j));
            Beta2n_TimeDiffSecond->Fill(ffBeta2n_TimeDiffSecond->at(j));
        }

        for (ULong_t j = 0; j < ffBeta1nBackward_TimeDiff->size(); ++j)
        {
            Beta1nBackward_TimeDiff->Fill(ffBeta1nBackward_TimeDiff->at(j));
            Beta1nBackward_Time_coinc->Fill(ffBeta1nBackward_Time_coinc->at(j));
        }

        for(ULong_t j = 0; j < ffBeta2nBackward_Time_coinc->size(); ++j)
        {
            Beta2nBackward_Time_coinc->Fill(ffBeta2nBackward_Time_coinc->at(j));
            Beta2nBackward_TimeDiffFirst->Fill(ffBeta2nBackward_TimeDiffFirst->at(j));
            Beta2nBackward_TimeDiffSecond->Fill(ffBeta2nBackward_TimeDiffSecond->at(j));
        }

        for(ULong_t j = 0; j < ffBeta3n_Time_coinc->size(); ++j)
        {
            Beta3n_Time_coinc->Fill(ffBeta3n_Time_coinc->at(j));
        }

        for(ULong_t j = 0; j < ff2n_TimeDiff->size(); ++j)
        {
            TwoNeutrons_TimeDiff->Fill(ff2n_TimeDiff->at(j));
            TwoNeutrons_Time_coinc->Fill(ff2n_Time_coinc->at(j));

            FirstNeutronCellGroup->Fill(ffFirstNeutronCellGroup->at(j));
            SecondNeutronCellGroup->Fill(ffSecondNeutronCellGroup->at(j));

            CellGroups->Fill(ffFirstNeutronCellGroup->at(j), ffSecondNeutronCellGroup->at(j));
        }

        for(ULong_t j = 0; j < ff1n_Time->size(); ++j)
        {
            OneNeutron_Time->Fill(ff1n_Time->at(j));
            OneNeutron_Cycle->Fill(ff1n_Cycle->at(j));
        }
    }

    //Normalization of cell groups

    Tetra_Rings->SetBinContent(3,Tetra_Rings->GetBinContent(3)/(7.*28.6));
    Tetra_Rings->SetBinContent(5,Tetra_Rings->GetBinContent(5)/(4.*29.7));
    Tetra_Rings->SetBinContent(7,Tetra_Rings->GetBinContent(7)/(7.*34.3));
    Tetra_Rings->SetBinContent(9,Tetra_Rings->GetBinContent(9)/(7.3*34.1));
    Tetra_Rings->SetBinContent(11,Tetra_Rings->GetBinContent(11)/(3.*34.7));
    Tetra_Rings->SetBinContent(13,Tetra_Rings->GetBinContent(13)/(7.*35.6));
    Tetra_Rings->SetBinContent(15,Tetra_Rings->GetBinContent(15)/(7.*36.1));
    Tetra_Rings->SetBinContent(17,Tetra_Rings->GetBinContent(17)/(8.*34.75));
    Tetra_Rings->SetBinContent(19,Tetra_Rings->GetBinContent(19)/(8.*37.9));
    Tetra_Rings->SetBinContent(21,Tetra_Rings->GetBinContent(21)/(7.*35.1));
    Tetra_Rings->SetBinContent(23,Tetra_Rings->GetBinContent(23)/(7.*35));
    Tetra_Rings->SetBinContent(25,Tetra_Rings->GetBinContent(25)/(7.*34.3));

    for(int i = 2; i <= 13; i++)
    {
        CellGroups->SetBinContent(2*i-1,3,CellGroups->GetBinContent(2*i-1,3)/(7.*28.6));
        CellGroups->SetBinContent(2*i-1,5,CellGroups->GetBinContent(2*i-1,5)/(4.*29.7));
        CellGroups->SetBinContent(2*i-1,7,CellGroups->GetBinContent(2*i-1,7)/(7.*34.3));
        CellGroups->SetBinContent(2*i-1,9,CellGroups->GetBinContent(2*i-1,9)/(7.3*34.1));
        CellGroups->SetBinContent(2*i-1,11,CellGroups->GetBinContent(2*i-1,11)/(3.*34.7));
        CellGroups->SetBinContent(2*i-1,13,CellGroups->GetBinContent(2*i-1,13)/(7.*35.6));
        CellGroups->SetBinContent(2*i-1,15,CellGroups->GetBinContent(2*i-1,15)/(7.*36.1));
        CellGroups->SetBinContent(2*i-1,17,CellGroups->GetBinContent(2*i-1,17)/(8.*34.75));
        CellGroups->SetBinContent(2*i-1,19,CellGroups->GetBinContent(2*i-1,19)/(8.*37.9));
        CellGroups->SetBinContent(2*i-1,21,CellGroups->GetBinContent(2*i-1,21)/(7.*35.1));
        CellGroups->SetBinContent(2*i-1,23,CellGroups->GetBinContent(2*i-1,23)/(7.*35));
        CellGroups->SetBinContent(2*i-1,25,CellGroups->GetBinContent(2*i-1,25)/(7.*34.3));
    }

    file->Write();
}