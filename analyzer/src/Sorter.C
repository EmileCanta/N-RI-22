#include "../include/Sorter.h"
 
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
    else if(time / 1e6 > bgd_time + coll_time) status = 1; //Acquisition
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
        //Single

    //Time
    fTetra_tSingle.clear();
    fBeta_tSingle.clear();
    fGe_tSingle.clear();

    //Energy
    fGe_ESingle.clear();

    //Cycles
    fTetra_Cycle.clear();
    fBeta_Cycle.clear();
    fGe_Cycle.clear();
    
    //Tetra cell groups and rings
    fTetra_CellGroups.clear();
    fTetra_Rings.clear();

        //Coincidences

    //Ge-Beta
    fGeBeta_ECond.clear();
    fGeBeta_tCond.clear();
    fGeBeta_tDiff.clear();

    //Ge-Tetra
    fGeTetra_ECond.clear();
    fGeTetra_tCond.clear();
    fGeTetra_tDiff.clear();
        
    //Beta-Tetra
    fBeta1n_tCond.clear();
    fBeta2n_tCond.clear();
    fBeta1nBackward_tCond.clear();
    fBeta2nBackward_tCond.clear();
    fBeta1n_tDiff.clear();
    fBeta2n_tDiffFirst.clear();
    fBeta2n_tDiffSecond.clear();
    fBeta1nBackward_tDiff.clear();
    fBeta2nBackward_tDiffFirst.clear();
    fBeta2nBackward_tDiffSecond.clear();

    //Tetra-Tetra
    fFirstNeutronCellGroup.clear();
    fSecondNeutronCellGroup.clear();
    fSecondNeut_tDiff.clear();
    fSecondNeut_tCond.clear();

    //Storing vectors
    fStoring1n_Time.clear();
    fStoring1n_TimeDiff.clear();
    fStoring2n_Time.clear();
    fStoring2n_TimeDiff.clear();
    fStoringFirstNeutronCellGroup.clear();
    fStoringSecondNeutronCellGroup.clear();
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

    if(det == 14)
    { 
        fBeta_Id = det;
        fBeta_Time = raw_time;
        fCoding = raw_coding;
        fStatus = GetStatus(raw_time);
    }
        
    if(det == 15)
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
      
        //Single
    
    //Time
    output_tree_single->Branch("Tetra_tSingle", &fTetra_tSingle);
    output_tree_single->Branch("Ge_tSingle", &fGe_tSingle);
    output_tree_single->Branch("Beta_tSingle", &fBeta_tSingle);

    //Energy
    output_tree_single->Branch("Ge_ESingle", &fGe_ESingle);

    //Cycles
    output_tree_single->Branch("Tetra_Cycle", &fTetra_Cycle);    
    output_tree_single->Branch("Ge_Cycle", &fGe_Cycle);
    output_tree_single->Branch("Beta_Cycle", &fBeta_Cycle);

    //Tetra cells and rings
    output_tree_single->Branch("Tetra_Rings", &fTetra_Rings);
    output_tree_single->Branch("Tetra_CellGroups", &fTetra_CellGroups);
    
        //Coincidences

    //Ge-Beta
    output_tree_coinc->Branch("GeBeta_ECond", &fGeBeta_ECond);
    output_tree_coinc->Branch("GeBeta_tCond", &fGeBeta_tCond);
    output_tree_coinc->Branch("GeBeta_tDiff", &fGeBeta_tDiff);
    
    //Ge-Tetra
    output_tree_coinc->Branch("GeTetra_ECond", &fGeTetra_ECond);
    output_tree_coinc->Branch("GeTetra_tCond", &fGeTetra_tCond);
    output_tree_coinc->Branch("GeTetra_tDiff", &fGeTetra_tDiff);

    //Beta-Tetra
    output_tree_coinc->Branch("Beta1n_tCond", &fBeta1n_tCond);
    output_tree_coinc->Branch("Beta2n_tCond", &fBeta2n_tCond);
    output_tree_coinc->Branch("Beta1nBackward_tCond", &fBeta1nBackward_tCond);
    output_tree_coinc->Branch("Beta2nBackward_tCond", &fBeta2nBackward_tCond);
    output_tree_coinc->Branch("Beta1n_tDiff", &fBeta1n_tDiff);
    output_tree_coinc->Branch("Beta2n_tDiffFirst", &fBeta2n_tDiffFirst);
    output_tree_coinc->Branch("Beta2n_tDiffSecond", &fBeta2n_tDiffSecond);
    output_tree_coinc->Branch("Beta1nBackward_tDiff", &fBeta1nBackward_tDiff);
    output_tree_coinc->Branch("Beta2nBackward_tDiffFirst", &fBeta2nBackward_tDiffFirst);
    output_tree_coinc->Branch("Beta2nBackward_tDiffSecond", &fBeta2nBackward_tDiffSecond);

    //Tetra-Tetra
    output_tree_coinc->Branch("SecondNeut_tCond", &fSecondNeut_tCond);
    output_tree_coinc->Branch("SecondNeut_tDiff", &fSecondNeut_tDiff);
    output_tree_coinc->Branch("FirstNeutronCellGroup", &fFirstNeutronCellGroup);
    output_tree_coinc->Branch("SecondNeutronCellGroup", &fSecondNeutronCellGroup);
}

//******************************************************************************
//******************************************************************************

void Sorter::FillSingleBranches()
{   
    std::cout << "Sorting single data" << std::endl;

    Double_t last_time = 0.;
    UInt_t fCycle = 1;

    ResetVar();
    
    for(fEntry = 0; fEntry < fEntries; fEntry++)
	{
        raw_data_tree->GetEntry(fEntry);

        if(fEntry == 0)
        {
            last_time = raw_time;
            fCycle = 1;
	    }

        ClearVectors();

        if(raw_det_nbr >= 1 && raw_det_nbr <= 12)
        {
            fTetra_tSingle.push_back(raw_time / 1.e9); //millisecond
            fTetra_Cycle.push_back(fCycle);

            Tetra_tSingle->Fill(raw_time / 1.e9);
            Tetra_Cycle->Fill(fCycle);
            
            for(UInt_t i = 1; i <= 13; i++)
            {
                if(raw_det_nbr == i)
                {
                    fTetra_CellGroups.push_back(i);
                    Tetra_CellGroups->Fill(i);
                }
            }

            if(raw_det_nbr >= 1 && raw_det_nbr <= 2)
            {
                fTetra_Rings.push_back(1);
                Tetra_Rings->Fill(1);
            }

            if(raw_det_nbr >= 3 && raw_det_nbr <= 5)
            {
                fTetra_Rings.push_back(2);
                Tetra_Rings->Fill(2);
            }

            if(raw_det_nbr >= 6 && raw_det_nbr <= 8)
            {
                fTetra_Rings.push_back(3);
                Tetra_Rings->Fill(3);
            }

            if(raw_det_nbr >= 9 && raw_det_nbr <= 12)
            {
                fTetra_Rings.push_back(4);
                Tetra_Rings->Fill(4);
            }

            if(raw_det_nbr == 13)
            {
                fTetra_Rings.push_back(0);
                Tetra_Rings->Fill(0);
            }
        }

        if(raw_det_nbr == 14)
        {
            fBeta_tSingle.push_back(raw_time / 1.e9); //millisecond
            fBeta_Cycle.push_back(fCycle);

            Beta_tSingle->Fill(raw_time / 1.e9);
            Beta_Cycle->Fill(fCycle);
        }

        if(raw_det_nbr == 15)
        {
            fGe_ESingle.push_back(Ge_alignement(raw_energy));
            fGe_tSingle.push_back(raw_time / 1.e9); //millisecond
            fGe_Cycle.push_back(fCycle);

            Ge_ESingle->Fill(Ge_alignement(raw_energy));
            Ge_tSingle->Fill(raw_time / 1.e9);
            Ge_Cycle->Fill(fCycle);
        }

        if(raw_time - last_time < 0) fCycle++;
        
        last_time = raw_time;
        
        output_tree_single->Fill();

        std::cout << std::setprecision(3) << std::setw(5) << (100.*fEntry/fEntries) << " %\r";
    }

    Ge_ESingle->Write();
    Ge_tSingle->Write();
    Ge_Cycle->Write();
    Beta_tSingle->Write();
    Beta_Cycle->Write();
    Tetra_tSingle->Write();
    Tetra_Cycle->Write();
    Tetra_Rings->Write();
    Tetra_CellGroups->Write();
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
                        fGeTetra_tDiff.push_back(TMath::Abs(raw_time - fTetra_Time) / 1e6); //microsecondes
				        fGeTetra_ECond.push_back(Ge_alignement(raw_energy));
				        fGeTetra_tCond.push_back(raw_time / 1.e9); //milisecondes;

                        GeTetra_tDiff->Fill(TMath::Abs(raw_time - fTetra_Time) / 1.e6);
                        GeTetra_ECond->Fill(Ge_alignement(raw_energy));
                        GeTetra_tCond->Fill(raw_time / 1.e9);
                        ESvsTTD->Fill(TMath::Abs(raw_time - fTetra_Time) / 1e6, Ge_alignement(raw_energy));
                    }
                 }
            }
	    }
	    
	    lastevent = fEntry - eventafter;
	    
	    output_tree_coinc->Fill();

        std::cout << std::setprecision(3) << std::setw(5) << (100.*fEntry/fEntries) << " %\r";
    }

    GeTetra_tDiff->Write();
    GeTetra_ECond->Write();
    GeTetra_tCond->Write();
    ESvsTTD->Write();
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
                        fGeBeta_tDiff.push_back(TMath::Abs(raw_time - fBeta_Time) / 1.e6); //microsecondes
                        fGeBeta_ECond.push_back(Ge_alignement(raw_energy));
                        fGeBeta_tCond.push_back(raw_time / 1.e9); //milisecondes;

                        GeBeta_tDiff->Fill(TMath::Abs(raw_time - fBeta_Time) / 1.e6); //microsecondes
                        GeBeta_ECond->Fill(Ge_alignement(raw_energy));
                        GeBeta_tCond->Fill(raw_time / 1.e9); //milisecondes;

                        ESvsBTD->Fill(TMath::Abs(raw_time - fBeta_Time) / 1.e6, Ge_alignement(raw_energy));
                    }
                 }
            }
	    }
	    
	    lastevent = fEntry + eventafter;
	    
	    output_tree_coinc->Fill();

        std::cout << std::setprecision(3) << std::setw(5) << (100.*fEntry/fEntries) << " %\r";
    }

    GeBeta_tDiff->Write();
    GeBeta_ECond->Write();
    GeBeta_tCond->Write();
    ESvsBTD->Write();
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
                    }
                 }
            }

            if(neutcount == 1)
            {
                fBeta1n_tCond.push_back(fStoring1n_Time.front());
                fBeta1n_tDiff.push_back(fStoring1n_TimeDiff.front());

                Beta1n_tCond->Fill(fStoring1n_Time.front());
                Beta1n_tDiff->Fill(fStoring1n_TimeDiff.front());
            }

            if(neutcount == 2)
            {
                fBeta2n_tCond.push_back(fStoring2n_Time.front());
                fBeta2n_tDiffFirst.push_back(fStoring1n_TimeDiff.front());
                fBeta2n_tDiffSecond.push_back(fStoring2n_TimeDiff.front());

                Beta2n_tCond->Fill(fStoring2n_Time.front());
                Beta2n_tDiffFirst->Fill(fStoring1n_TimeDiff.front());
                Beta2n_tDiffSecond->Fill(fStoring2n_TimeDiff.front());
            }
	    }

	    lastevent = fEntry + eventafter;
	    
	    output_tree_coinc->Fill();

        std::cout << std::setprecision(3) << std::setw(5) << (100.*fEntry/fEntries) << " %\r";
    }

    Beta1n_tCond->Write();
    Beta1n_tDiff->Write();
    Beta2n_tCond->Write();
    Beta2n_tDiffFirst->Write();
    Beta2n_tDiffSecond->Write();
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

        neutcount = 0;

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
                    }
                 }
            }

            if(neutcount == 1)
            {
                fBeta1nBackward_tCond.push_back(fStoring1n_Time.front());
                fBeta1nBackward_tDiff.push_back(fStoring1n_TimeDiff.front());

                Beta1nBackward_tCond->Fill(fStoring1n_Time.front());
                Beta1nBackward_tDiff->Fill(fStoring1n_TimeDiff.front());
            }

            if(neutcount == 2)
            {
                fBeta2nBackward_tCond.push_back(fStoring2n_Time.front());
                fBeta2nBackward_tDiffFirst.push_back(fStoring1n_TimeDiff.front());
                fBeta2nBackward_tDiffSecond.push_back(fStoring2n_TimeDiff.front());

                Beta2nBackward_tCond->Fill(fStoring2n_Time.front());
                Beta2nBackward_tDiffFirst->Fill(fStoring1n_TimeDiff.front());
                Beta2nBackward_tDiffSecond->Fill(fStoring2n_TimeDiff.front());
            }
	    }

	    lastevent = fEntry - eventafter;
	    
	    output_tree_coinc->Fill();

        std::cout << std::setprecision(3) << std::setw(5) << (100.*fEntry/fEntries) << " %\r";
    }

    Beta1nBackward_tCond->Write();
    Beta1nBackward_tDiff->Write();
    Beta2nBackward_tCond->Write();
    Beta2nBackward_tDiffFirst->Write();
    Beta2nBackward_tDiffSecond->Write();
}

//******************************************************************************
//******************************************************************************

void Sorter::FillNeutronNeutronCoincBranches(Double_t coincwindow)
{
    std::cout << "Sorting neutron-neutron coinc data" << std::endl;

    lastevent = 0;

    lastneutron = 0;

    ResetVar();

	for(fEntry = 0; fEntry < fEntries; fEntry++)
	{
        raw_data_tree->GetEntry(fEntry);

        ClearVectors();
        
        SetVar(raw_det_nbr);
        
        eventafter = 0;

        neutcount = 0;

        neut_id = 0;

        if(fTetra_Id >= 1 && fTetra_Id <= 12)
        {  
            if(fEntry < lastneutron - 1) continue;

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

                            neut_id = eventafter;

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
                fSecondNeut_tDiff.push_back(fStoring2n_TimeDiff.front());
                fSecondNeut_tCond.push_back(fStoring2n_Time.front());
                fFirstNeutronCellGroup.push_back(fStoringFirstNeutronCellGroup.front());
                fSecondNeutronCellGroup.push_back(fStoringSecondNeutronCellGroup.front());

                SecondNeut_tDiff->Fill(fStoring2n_TimeDiff.front());
                SecondNeut_tCond->Fill(fStoring2n_Time.front());
                FirstNeutronCellGroup->Fill(fStoringFirstNeutronCellGroup.front());
                SecondNeutronCellGroup->Fill(fStoringSecondNeutronCellGroup.front());

                CellGroups->Fill(fStoringFirstNeutronCellGroup.front(), fStoringSecondNeutronCellGroup.front());
            }
        }

        lastevent = fEntry + eventafter;

        lastneutron = fEntry + neut_id;

        output_tree_coinc->Fill();

        std::cout << std::setprecision(3) << std::setw(5) << (100.*fEntry/fEntries) << " %\r";
    }

    SecondNeut_tDiff->Write();
    SecondNeut_tCond->Write();
    FirstNeutronCellGroup->Write();
    SecondNeutronCellGroup->Write();
    CellGroups->Write();
}

//******************************************************************************
//******************************************************************************

/*void Sorter::Histogrammer(const char* OutputFileName)
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

    TH1I* Tetra_Cycle = new TH1I("Tetra_Cycle", "Tetra_Cycle", 1000,0,1000);
    TH1I* Beta_Cycle = new TH1I("Beta_Cycle", "Beta_Cycle", 1000,0,1000);
    TH1I* Ge_Cycle = new TH1I("Ge_Cycle", "Ge_Cycle", 1000,0,1000);

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

    TH1D* Beta1nBackward_Time_coinc = new TH1D("Beta1nBackward_Time_coinc","Beta1nBackward_Time_coinc", 30000, 0, 30000);
    TH1D* Beta2nBackward_Time_coinc = new TH1D("Beta2nBackward_Time_coinc","Beta2nBackward_Time_coinc", 30000, 0, 30000);

    TH1D* GeBeta_E_coinc = new TH1D("GeBeta_E_coinc", "GeBeta_E_coinc", 7000, 0, 7000);
    TH1D* GeTetra_E_coinc = new TH1D("GeTetra_E_coinc", "GeTetra_E_coinc", 7000, 0, 7000);

    TH1D* FirstNeutronCellGroup = new TH1D("FirstNeutronCellGroup", "FirstNeutronCellGroup", 28, 0, 14);
    TH1D* SecondNeutronCellGroup = new TH1D("SecondNeutronCellGroup", "SecondNeutronCellGroup", 28, 0, 14);

    TH1D* Tetra_Rings = new TH1D("Tetra_Rings", "Tetra_Rings", 28, 0, 14);
    TH1D* Tetra_CellGroups = new TH1D("Tetra_CellGroups", "Tetra_CellGroups", 28, 0, 14);

    //Bidim
    TH2D* ESvsBTD = new TH2D("ESvsBTD", "ESvsBTD", 5000, 0, 10, 7000, 0, 7000);
    TH2D* ESvsTTD = new TH2D("ESvsTTD", "ESvsTTD", 1000, 0, 200, 7000, 0, 7000);

    TH2D* CellGroups = new TH2D("CellGroups", "CellGroups", 28, 0, 14, 28, 0, 14);

    std::vector<Double_t> *ffGe_E_single = 0;

    std::vector<Double_t> *ffGe_Time_single = 0;
    std::vector<Double_t> *ffBeta_Time_single = 0;
    std::vector<Double_t> *ffTetra_Time_single = 0;

    std::vector<UInt_t> *ffTetra_Cycle = 0;
    std::vector<UInt_t> *ffBeta_Cycle = 0;
    std::vector<UInt_t> *ffGe_Cycle = 0;

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

    std::vector<Double_t> *ffBeta1nBackward_Time_coinc = 0;
    std::vector<Double_t> *ffBeta2nBackward_Time_coinc = 0;

    std::vector<Double_t> *ffGeBeta_E_coinc = 0;
    std::vector<Double_t> *ffGeTetra_E_coinc = 0;

    std::vector<Double_t> *ffFirstNeutronCellGroup = 0;
    std::vector<Double_t> *ffSecondNeutronCellGroup = 0;

    std::vector<Double_t> *ffTetra_Rings = 0;
     std::vector<Double_t> *ffTetra_CellGroups = 0;

    tsingle->SetBranchAddress("Ge_E_single",&ffGe_E_single);

    tsingle->SetBranchAddress("Ge_Time_single",&ffGe_Time_single);
    tsingle->SetBranchAddress("Beta_Time_single",&ffBeta_Time_single);
    tsingle->SetBranchAddress("Tetra_Time_single",&ffTetra_Time_single);

    tsingle->SetBranchAddress("Tetra_Cycle",&ffTetra_Cycle);
    tsingle->SetBranchAddress("Beta_Cycle",&ffBeta_Cycle);
    tsingle->SetBranchAddress("Ge_Cycle",&ffGe_Cycle);

    tsingle->SetBranchAddress("Tetra_Rings",&ffTetra_Rings);
    tsingle->SetBranchAddress("Tetra_CellGroups",&ffTetra_CellGroups);

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
            Tetra_CellGroups->Fill(ffTetra_CellGroups->at(j));
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
    }

    //Normalization of cell groups

    Tetra_CellGroups->SetBinContent(3,Tetra_CellGroups->GetBinContent(3)/(7.*28.6));
    Tetra_CellGroups->SetBinContent(5,Tetra_CellGroups->GetBinContent(5)/(4.*29.7));
    Tetra_CellGroups->SetBinContent(7,Tetra_CellGroups->GetBinContent(7)/(7.*34.3));
    Tetra_CellGroups->SetBinContent(9,Tetra_CellGroups->GetBinContent(9)/(7.3*34.1));
    Tetra_CellGroups->SetBinContent(11,Tetra_CellGroups->GetBinContent(11)/(3.*34.7));
    Tetra_CellGroups->SetBinContent(13,Tetra_CellGroups->GetBinContent(13)/(7.*35.6));
    Tetra_CellGroups->SetBinContent(15,Tetra_CellGroups->GetBinContent(15)/(7.*36.1));
    Tetra_CellGroups->SetBinContent(17,Tetra_CellGroups->GetBinContent(17)/(8.*34.75));
    Tetra_CellGroups->SetBinContent(19,Tetra_CellGroups->GetBinContent(19)/(8.*37.9));
    Tetra_CellGroups->SetBinContent(21,Tetra_CellGroups->GetBinContent(21)/(7.*35.1));
    Tetra_CellGroups->SetBinContent(23,Tetra_CellGroups->GetBinContent(23)/(7.*35));
    Tetra_CellGroups->SetBinContent(25,Tetra_CellGroups->GetBinContent(25)/(7.*34.3));

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
}*/