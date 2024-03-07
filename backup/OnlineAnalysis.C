/*
 *  OnlineAnalysis.C
 *
 */
#include "OnlineAnalysis.h"




OnlineAnalysis::OnlineAnalysis(const char* InputFileName, const char* OutputFileName,
                               const char* CalibFileGe,
                               UChar_t collection_time, UChar_t acquisition_time)
{
    cout<< "hello im here" << endl;
    gROOT->ProcessLine("#include <vector>");
    //	std::cout << "Analysis code" << std::endl;
    timer.Reset();
    timer.Start();

    coll_time = collection_time;
    acq_time = acquisition_time;

    //Set_Nb_Ge_det(2);                    // Set Number of Germanium detector
    //Set_Nb_Veto_det(1);                  // Set Number of Particle detector
    //Ini_align_param_ge(CalibFileGe);     // Reading of calibration files

    Load_rawdata(InputFileName);         // Loading Narval Tree
    SetBranches(OutputFileName);   // Initialization of output Tree
    FillBranches();          // Fill output Tree
    output_tree->Print();
    timer.Stop();
    cout <<endl<<"\t\t";
    timer.Print();
    gROOT->ProcessLine(".q");
}

//******************************************************************************
//******************************************************************************

OnlineAnalysis::~OnlineAnalysis()
{
    ;
}


//******************************************************************************
//******************************************************************************



//******************************************************************************
//******************************************************************************




void OnlineAnalysis::Load_rawdata (const char* InputFileName){

    cout <<'\t' << "Loading Narval Tree" << endl;

    data_file   = TFile::Open(InputFileName);
    while (!data_file) {
        std::cerr << "Input file doesn't exist, please check the path or press q to exit!" << std::endl;
        string input;
        std::cin >> input;
        if (strcmp(input.c_str(),"q")==0) {
            gROOT->ProcessLine(".q");
        }
        data_file   = TFile::Open(input.c_str());
    }

    data_file->Print();

    raw_data_tree = (TTree*)data_file->Get("Narval_tree");
    if (!raw_data_tree) {
        std::cerr << "Input Narval Tree doesn't exist, please check the input file!" << std::endl;
        gROOT->ProcessLine(".q");
    }

    raw_data_tree->SetBranchAddress("Energy",&raw_energy);
    raw_data_tree->SetBranchAddress("Time",  &raw_time);
    raw_data_tree->SetBranchAddress("Marker",&raw_marker);
    raw_data_tree->SetBranchAddress("Coding",&raw_coding);
    raw_data_tree->SetBranchAddress("Det_nbr",&raw_det_nbr);



    fEntries = raw_data_tree -> GetEntries();

    std::cout << "Entries to be sorted : " << fEntries << std::endl;
}

//******************************************************************************
//******************************************************************************


void OnlineAnalysis::Ini_align_param_ge(const char* CalibFileGe){

    // Parameters for the alignement of ge detector in energy
    ifstream align_ge_param;
    align_ge_param.open(CalibFileGe);
    if (align_ge_param.good()) {
        std::cout << "Calibration parameters :" << std::endl;
        Int_t det = 0, cal_order = 0;
        for (Int_t j=0; j<Get_Nb_Ge_det(); j++) {
            align_ge_param >> det >> cal_order;
            std::cout << det << "\t" << cal_order << "\t";
            for (Int_t k=0; k<cal_order; k++) {
                Double_t tmp;
                align_ge_param >> tmp;
                std::cout << tmp << "\t" ;
                tmp_vector.push_back(tmp);
            }
            std::cout << std::endl;
            ge_align_param_e.push_back(tmp_vector);
            tmp_vector.clear();
        }
    }
}


//******************************************************************************
//******************************************************************************


Float_t OnlineAnalysis::Ge_alignement(Int_t det_id, UInt_t channel){

    Float_t channel_aligned=0;
    Float_t ch = ((float)channel+gRandom->Uniform(1.0)-.5);
    channel_aligned= -8.267+0.19731438802456427*ch;

    return channel_aligned;
}


//******************************************************************************
//******************************************************************************


Double_t OnlineAnalysis::DopplerCorrection(Double_t angle, Float_t part_velocity, Float_t gammaE){

    return gammaE*1./sqrt(1-part_velocity*part_velocity)*(1.-part_velocity*TMath::Cos(angle));
}

//******************************************************************************
//******************************************************************************

void OnlineAnalysis::SetBranches(const char* out_name)
{
    cout << "Initialization of output TTree :" << endl;
    cout << '\t' << "Creating branches ..."; ;

    output_file = new TFile(out_name,"recreate");
    output_tree = new TTree("tvec","Tree with vectors");
    if (1==1){
    //Single event branches
    output_tree->Branch("Ge_E_single"   , &fGe_E_single    , "Ge_E_single/F"     );
    output_tree->Branch("Ge_Id_single"  , &fGe_Id_single   , "Ge_Id_single/b"    );
    output_tree->Branch("Ge_Time_single", &fGe_Time_single , "Ge_Time_single/D"  );
    output_tree->Branch("TETRA_E_single"   , &fTETRA_E_single    , "TETRA_E_single/F"     );
    output_tree->Branch("TETRA_Time_single", &fTETRA_Time_single , "TETRA_Time_single/D"  );
    output_tree->Branch("Veto_single"   , &fVeto_single    , "Veto_E_single/i"   );
    output_tree->Branch("Veto_Id"       , &fVeto_Id        , "Veto_Id_single/b"  );
    output_tree->Branch("Veto_Time"     , &fVeto_Time      , "Veto_Time_single/D");
    output_tree->Branch("Marker"        , &fMarker         , "Marker/b"          );
    output_tree->Branch("Coding"        , &fCoding         , "Coding/i"	         );
    output_tree->Branch("Status"        , &fStatus         , "Status/i"          );
    output_tree->Branch("Cycle"         , &fCycle           , "Cycle/I"          );

    //Coinc branches
    output_tree->Branch("Ge_E_coinc"    , &fGe_E_coinc			        );//plastic gamma coincidence
    output_tree->Branch("Ge_Id_coinc"   , &fGe_Id_coinc			        );
    output_tree->Branch("Ge_Time_coinc" , &fGe_Time_coinc		        );
    output_tree->Branch("TETRA_E_coinc"    , &fTETRA_E_coinc			        );//plastic gamma coincidence
    output_tree->Branch("TETRA_Time_coinc" , &fTETRA_Time_coinc		        );
    output_tree->Branch("Veto_E_coinc"  , &fVeto_E_coinc                        );//plastic energy coincidence with gamma
    output_tree->Branch("Veto_Id_coinc" , &fVeto_Id_coinc                       );
    output_tree->Branch("Veto_Time_coinc", &fVeto_Time_coinc                    );
    output_tree->Branch("Nb_gamma_coinc", &fnb_gamma_coinc, "Nb_gamma_coinc/I"  );
    output_tree->Branch("Time_diff"     , &fTimeDiff                            );
    output_tree->Branch("TETRATime_diff"     , &fTETRATimeDiff                            );
    output_tree->Branch("Marker_coinc"  , &fMarker_coinc                        );
    output_tree->Branch("Coding_coinc"  , &fCoding_coinc                        );
    output_tree->Branch("Status_coinc"  , &fStatus_coinc                        );
    output_tree->Branch("Cycle_coinc"   , &fCycle_coinc                         );
    output_tree->Branch("TETRAhit"   , &TETRAhit                         );
    output_tree->Branch("GeTETRA_E_coinc"   , &fGeTETRA_E_coinc                       );
    output_tree->Branch("GeTETRA_Id_coinc"  , &fGeTETRA_Id_coinc                      );
    output_tree->Branch("GeTETRA_Time_coinc", &fGeTETRA_Time_coinc                    );
    }

}


//******************************************************************************
//******************************************************************************

UChar_t OnlineAnalysis::GetStatus(Double_t time)
{
    UChar_t status = -10;
    if(time > 0 && time<=coll_time*1.E9) status = 1; //Collection
    else if (time>coll_time*1.E9) status = 1 ; //Acquisition
    return status;

}

//******************************************************************************
//******************************************************************************
void OnlineAnalysis::FillBranches()
{
    Int_t counter=0;
    Int_t lastevent=0;
    int BETADET = 3;

    // bool first = true;
    std::cout << "Sorting data: " << std::endl;
    UInt_t previous_coding = 0 ;
    Double_t last_time     = 0.;


    if (1==1){
    for(fEntry = 0 ; fEntry < fEntries ; fEntry++)	{
        raw_data_tree->GetEntry(fEntry);
	  if (fEntry==0) {
            last_time = raw_time ;
            fCycle=0;
	    }

        // clearing of all std::vector
        fGe_E_coinc.clear()       ;
        fGe_Id_coinc.clear()      ;
        fGe_Time_coinc.clear()    ;
        fVeto_E_coinc.clear()     ;
        fVeto_Id_coinc.clear()    ;
        fVeto_Time_coinc.clear()  ;
        fTimeDiff.clear()         ;
        fTETRATimeDiff.clear()         ;
        fMarker_coinc.clear()     ;
        fCoding_coinc.clear()     ;
        fStatus_coinc.clear()     ;
        fCycle_coinc.clear()      ;
        fTETRA_E_coinc.clear()       ;
        fTETRA_Time_coinc.clear()    ;
        fGeTETRA_E_coinc.clear()     ;
        fGeTETRA_Id_coinc.clear()    ;
        fGeTETRA_Time_coinc.clear()  ;
        //TETRAhit =false;



        fnb_gamma_coinc = 0       ;


        // Single Event
        if (raw_det_nbr==1) { // Ge Event
            fGe_E_single    =  Ge_alignement(raw_det_nbr,raw_energy) ;
            fGe_Id_single   =  raw_det_nbr        ;
            fGe_Time_single =  raw_time/1000000000.           ;
            fMarker         =  raw_marker         ;
            fCoding         =  raw_coding         ;
            fStatus         =  GetStatus(raw_time);
        }
        /*if (raw_det_nbr==0) { // LaBr event
            fTETRA_E_single    = raw_energy ;//LaBr_alignement(raw_energy) ;
            fTETRA_Time_single =  raw_time/1000000000.                  ;
            fStatus         =  GetStatus(raw_time)       ;
        }*/

        if (raw_det_nbr==BETADET) { // Veto Event
            fVeto_single    =  raw_energy          ;
            fVeto_Id        =  raw_det_nbr         ;
            fVeto_Time      =  raw_time/1000000000.            ;
            fMarker         =  raw_marker          ;
            fCoding         =  raw_coding          ;
            fStatus         =  GetStatus(raw_time) ;
        }



        if (raw_time - last_time < 0) fCycle++;
        last_time = raw_time ;
        // Coinc Event
        int eventafter  =0;
        int eventbefore =0;



        if (fVeto_Id==BETADET ) { // BETA
            if(fEntry<lastevent)  continue;
            //std::cout << "qsfghj " << raw_det_nbr << std::endl;
            while(TMath::Abs(raw_time-fVeto_Time*1000000000.) < 90000000){



	      if (raw_det_nbr == BETADET) break;
                // Search backwards
               // std::cout << "Event before " << eventbefore << std::endl;
	      raw_data_tree->GetEntry(fEntry-eventbefore++);

	      if(fEntry-eventbefore==0 || fEntry-eventbefore == lastevent) break;
	      if(fCoding != raw_coding) break;
	      if(fStatus != GetStatus(raw_time)) break;
	      if(TMath::Abs(raw_time-fVeto_Time*1000000000.) > 90000000) continue;
	      if(TMath::Abs(raw_time-fVeto_Time*1000000000.) < 90000000){
		if(TMath::Abs(raw_time-fVeto_Time*1000000000.)<90000000)
		  {
		    if (raw_det_nbr==1) { // Ge in coinc
		      //  cout <<"evviva"<< endl;

		      fGe_E_coinc.push_back(Ge_alignement(raw_det_nbr,raw_energy))  ;
		      fGe_Id_coinc.push_back(raw_det_nbr);
		      fGe_Time_coinc.push_back(raw_time/1000000000.);
		      fTimeDiff.push_back(TMath::Abs(raw_time-fVeto_Time*1000000000.));

		      //fGeTETRA_E_coinc.push_back(Ge_alignement(raw_det_nbr,raw_energy))  ;
		      // fGeTETRA_Id_coinc.push_back(raw_det_nbr);
		      //fGeTETRA_Time_coinc.push_back(raw_time);

		    }

        /*else if (raw_det_nbr==0) { // LaBr in coinc
		      fTETRA_E_coinc.push_back(raw_energy)  ;
		      fTETRA_Time_coinc.push_back(raw_time/1000000000.);
          fTETRA_Id_coinc.push_back(raw_det_nbr);
          fTETRATimeDiff.push_back(TMath::Abs(raw_time-fVeto_Time*1000000000.));
          TETRAhit = true;
		      //TETRA_coinc_cal->Fill(TETRA_alignement(raw_energy));
		      //TETRA_coinc->Fill(raw_energy);
		      //fGeTETRA_E_coinc.push_back(TETRA_alignement(raw_energy))  ;
		      //fGeTETRA_Id_coinc.push_back(raw_det_nbr);
		      // fGeTETRA_Time_coinc.push_back(raw_time);
        }*/


                    fMarker_coinc.push_back(raw_marker);
                    fCoding_coinc.push_back(raw_coding);
                    fStatus_coinc.push_back(GetStatus(raw_time));
                    fnb_gamma_coinc++;
		  }
	      }
	    }
            while(TMath::Abs(raw_time-fVeto_Time*1000000000.)<90000000){
                // Search forwards
	      	if (raw_det_nbr == BETADET) break;
                raw_data_tree->GetEntry(fEntry+eventafter++);

                if(fEntry+eventafter>fEntries) break;
                if(fCoding != raw_coding) break;
                if(fStatus != GetStatus(raw_time)) break;
                if(TMath::Abs(raw_time-fVeto_Time*1000000000.) > 30000000) continue;
                if(TMath::Abs(raw_time-fVeto_Time*1000000000.) < 30000000){
		  if(TMath::Abs(raw_time-fVeto_Time*1000000000.)<3000000)
      {
        if (raw_det_nbr==1) { // Ge in coinc
          //  cout <<"evviva"<< endl;

          fGe_E_coinc.push_back(Ge_alignement(raw_det_nbr,raw_energy))  ;
          fGe_Id_coinc.push_back(raw_det_nbr);
          fGe_Time_coinc.push_back(raw_time/1000000000.);
          fTimeDiff.push_back(TMath::Abs(raw_time-fVeto_Time*1000000000.));

          //fGeTETRA_E_coinc.push_back(Ge_alignement(raw_det_nbr,raw_energy))  ;
          // fGeTETRA_Id_coinc.push_back(raw_det_nbr);
          //fGeTETRA_Time_coinc.push_back(raw_time);
        }

        /*else if (raw_det_nbr==0) { // LaBr in coinc
          fTETRA_E_coinc.push_back(raw_energy)  ;
          fTETRA_Time_coinc.push_back(raw_time/1000000000.);
          fTETRA_Id_coinc.push_back(raw_det_nbr);
          fTETRATimeDiff.push_back(TMath::Abs(raw_time-fVeto_Time*1000000000.));
          TETRAhit = true;
          //TETRA_coinc_cal->Fill(TETRA_alignement(raw_energy));
          //TETRA_coinc->Fill(raw_energy);
          //fGeTETRA_E_coinc.push_back(TETRA_alignement(raw_energy))  ;
          //fGeTETRA_Id_coinc.push_back(raw_det_nbr);
          // fGeTETRA_Time_coinc.push_back(raw_time);
        }*/

                    fMarker_coinc.push_back(raw_marker);
                    fCoding_coinc.push_back(raw_coding);
                    fStatus_coinc.push_back(GetStatus(raw_time));
                    fnb_gamma_coinc++;
      }
		}

            }
}


        lastevent=fEntry+eventafter;

	/*
	if (CL1_Id.size()>1 && CL1_4==0) {
          //std::cout << "Clover 2 hitted " << cluster2hit.size() << " times" << std::endl;
          Float_t tmp_energy=0;
	  UChar_t tmp_Marker=0;
          for (UInt_t k=0; k<CL1_Id.size(); k++)
           {
      tmp_energy+=fGe_E_coinc[CL1_Id[k]];
      fGe_E_coinc[CL1_Id[k]]=-10.;
      if(tmp_Marker==1 || fMarker_coinc[CL1_Id[k]]==1){tmp_Marker=1;}
      fMarker_coinc[CL1_Id[k]]=1;
          }
      fGe_E_coinc[CL1_Id[0]]=tmp_energy;
      fMarker_coinc[CL1_Id[0]]=tmp_Marker;
        }
	if(CL1_4==1 && CL1_Id.size()>1){
	  for (UInt_t k=0; k<CL1_Id.size(); k++)
	    {fGe_E_coinc[CL1_Id[k]]=-10; }
	}

        //Clover 2
        if (CL2_Id.size()>1) {
          //std::cout << "Clover 2 hitted " << cluster2hit.size() << " times" << std::endl;
          Float_t tmp_energy=0;
	  UChar_t tmp_Marker=0;
          for (UInt_t k=0; k<CL2_Id.size(); k++)
           {
      tmp_energy+=fGe_E_coinc[CL2_Id[k]];
      fGe_E_coinc[CL2_Id[k]]=-10.;
      if(tmp_Marker==1 || fMarker_coinc[CL2_Id[k]]==1){tmp_Marker=1;}
      fMarker_coinc[CL2_Id[k]]=1;
          }
      fGe_E_coinc[CL2_Id[0]]=tmp_energy;
      fMarker_coinc[CL2_Id[0]]=tmp_Marker;
      }*/
	/*

        if(fVeto_Id_coinc.size()>1)
        {
          for(int i=0;i<fVeto_Id_coinc.size();i++)
            {
              if(fVeto_Id_coinc[i]==11)
                {
                  for(int k; fGe_E_coinc.size();k++)
                    {
                      if(fGe_Id_coinc[k]==10)
                        {
                          fGe_E_coinc[k]=-10;
                        }
                    }
                }
              if(fVeto_Id_coinc[i]==13)
                {
                  for(int k;fGe_E_coinc.size();k++)
                    {
                      if(fGe_Id_coinc[k]==1 || fGe_Id_coinc[k]==2 || fGe_Id_coinc[k]==3 || fGe_Id_coinc[k]==4)
                        {
                          fGe_E_coinc[k]=-10;
                        }
                    }
                }
              if(fVeto_Id_coinc[i]==14)
                {
                  for(int k; fGe_E_coinc.size();k++)
                    {
                      if(fGe_Id_coinc[k]==11)
                        {
                          fGe_E_coinc[k]=-10;
                        }
                    }
                }
              if(fVeto_Id_coinc[i]==15)
                {
                  for(int k; fGe_E_coinc.size();k++)
                    {
                      if(fGe_Id_coinc[k]==5 || fGe_Id_coinc[k]==6 || fGe_Id_coinc[k]==7 || fGe_Id_coinc[k]==8)
                        {
                          fGe_E_coinc[k]=-10;
                        }
                    }
                }
              if(fVeto_Id_coinc[i]==16){continue;}
              if(fVeto_Id_coinc[i]==17){continue;}
            }
        }
*/

        output_tree->Fill();
        if (1==1) {
            std::cout << std::setprecision(3) << std::setw(5)
            << (100.*fEntry/fEntries) << " %\r";
            if (counter == 5E4)
                counter=0;
        }
        counter++;
	// previous_coding = raw_coding;
    }
    }

    /*  for (int h = 0; h < 2; h++)
       {
         delete[] Coinc_Windows[h];
       }
    delete[] Coinc_Windows;*/
  }









//******************************************************************************
//******************************************************************************