/*
 *  OnlineAnalysis.C
 *
 */

#include <vector>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <Riostream.h>

#include "TFile.h"
#include "TTree.h"
#include "TStopwatch.h"
#include "TVector3.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TROOT.h"
#include "TDirectory.h"
#include "TCutG.h"
#include "TRandom.h"


using namespace std;

class OnlineAnalysis{

 private:
  TStopwatch timer;

 private:

  TFile *data_file        ; // Input File
  TTree *raw_data_tree    ; // Input TTree


  Long64_t fEntries       ; // Total number of entries
  Long64_t fEntry         ; // Entry index

  //TFile *calib_file	  ; // Calibration file

  TFile *output_file      ; // Output File
  TTree *output_tree      ; // Output TTree




 private: // single events

  // Ge events
  Float_t  fGe_E_single    ; // Energy
  UChar_t  fGe_Id_single   ; // #Detector firing
  Double_t fGe_Time_single ; // Time

  Float_t  fParne_E_single    ; // Energy
  UChar_t  fParne_Id_single   ; // #Detector firing
  Double_t fParne_Time_single ; // Time

   Float_t  fParne_ESi_single    ; // Energy
  UChar_t  fParne_IdSi_single   ; // #Detector firing
  Double_t fParne_TimeSi_single ; // Time

  // TETRA
  Float_t  fTETRA_E_single    ;  // Energy
  Double_t fTETRA_Time_single ;  //Time

  // Veto events
  UInt_t   fVeto_single   ; // Energy
  UChar_t  fVeto_Id       ; // #Detector firing
  Double_t fVeto_Time     ; // Time
  Int_t fCycle            ;


 protected:


  Int_t N_Ge_det          ; // Number of Ge detectors
  Int_t N_Veto_det        ; // Number of Particle detectors
  Double_t start_cycle    ;
  UChar_t  coll_time      ;
  UChar_t  acq_time       ;

  // Marker of BGO
  UChar_t     fMarker     ; //
  // Coding
  UInt_t      fCoding     ; //
  // Status (Collection or Acquisition)
  UChar_t     fStatus     ;

 private: // coincidence events
  // Ge Events
  std::vector<Double_t> fGe_E_coinc      ; // Energy
  std::vector<UInt_t>   fGe_Id_coinc     ; // # Detector firing
  std::vector<Double_t> fTimeDiff        ; // Index diff of the events in coincidence in the TTree
  std::vector<Double_t> fTETRATimeDiff        ; // Index diff of the events in coincidence in the TTree
  std::vector<Double_t> fGe_Time_coinc   ;
  Int_t                 fnb_gamma_coinc  ; // Number of gamma event
  std::vector<Double_t> CL1_E ;
  Int_t nb_1;
  Int_t nb_2;
  bool TETRAhit;

  std::vector<Double_t> CL1_Time ;
  std::vector<UInt_t>   CL1_Id ;
  std::vector<Double_t> CL2_E ;
  std::vector<Double_t> CL2_Time ;
  std::vector<UInt_t>   CL2_Id ;

  // Veto events
  std::vector<UInt_t>    fVeto_E_coinc    ; //
  std::vector<UInt_t>    fVeto_Id_coinc   ; // # Detector firing
  std::vector<Double_t>  fVeto_Time_coinc ; //

  // TETRA Events
  std::vector<Double_t> fTETRA_E_coinc      ; // Energy
  std::vector<Double_t> fTETRA_Time_coinc   ; // Index diff of the events in coincidence in the TTree
  std::vector<Double_t> fTETRA_Id_coinc   ; // Index diff of the events in coincidence in the TTree

// TETRA Events + Ge Events in coinc
  std::vector<Double_t> fGeTETRA_E_coinc      ; // Energy
  std::vector<UInt_t>   fGeTETRA_Id_coinc     ; // # Detector firing
  std::vector<Double_t> fGeTETRA_Time_coinc   ;

  std::vector<UChar_t>   fMarker_coinc   ;
  std::vector<UInt_t>    fCoding_coinc   ;
  std::vector<UChar_t>   fStatus_coinc   ;
  std::vector<Int_t>     fCycle_coinc    ;


 private: //raw data tree
  UInt_t   raw_energy  ;
  Double_t raw_time    ;
  UChar_t  raw_marker  ;
  UInt_t   raw_coding  ;
  UChar_t  raw_det_nbr ;

 private:
  std::vector<std::vector<Double_t> > ge_align_param_e  ;
  std::vector<Double_t> parne_align_param_e ;
  std::vector<Double_t> tmp_vector;


 public:

  OnlineAnalysis(const char*, const char*, const char*, UChar_t, UChar_t);
  //OnlineAnalysis(const char* InputFileName, const char* OutputFileName, const char* CalibFile, bool UseParticleGate==true);
  virtual ~OnlineAnalysis();

  void Load_rawdata (const char* InputFileName);
  //virtual void SetCoinc_window(Double_t time_diff){Coinc_Windows=time_diff;};
  //double GetCoinc_window(){return Coinc_Windows;};

  //double Get_Min(Int_t det_nbr){return Coinc_Windows[0][det_nbr-1];}
  //double Get_Max(Int_t det_nbr){return Coinc_Windows[1][det_nbr-1];}
  //double Get_Mingg(Int_t det_nbr1,Int_t det_nbr2){return AddbackWindows[0][det_nbr1-1][det_nbr2-1];}
  //double Get_Maxgg(Int_t det_nbr1,Int_t det_nbr2){return AddbackWindows[1][det_nbr1-1][det_nbr2-1];}

  virtual void Set_Nb_Ge_det(Int_t Ge_det){N_Ge_det=Ge_det;};
  Int_t Get_Nb_Ge_det(){return N_Ge_det;};

  virtual void Set_Nb_Veto_det(Int_t Veto_det){N_Veto_det=Veto_det;};
  Int_t Get_Nb_Veto_det(){return N_Veto_det;};

 protected:
  TH1I **Ge_coinc         ;
  TH1F **Ge_coinc_cal     ;
  TH1F *TETRA_coinc_cal      ;
  TH1I *TETRA_coinc          ;
  TH1F ***ge_spec         ;
  TH1F **ge_time_spec     ;
  TH2F *Mat_gg            ;

 protected:
  virtual void SetBranches(const char* out_name)  ; // Initialize TTree
  virtual void FillBranches() ; // Filling the TTree
  virtual Float_t Ge_alignement(Int_t , UInt_t); // Alignement of Ge detector
  virtual void Ini_align_param_ge(const char*) ;
  virtual Double_t DopplerCorrection(Double_t angle, Float_t part_velocity, Float_t gammaE); //
  virtual UChar_t GetStatus(Double_t time);
  virtual void SetStartCycle(Double_t start) {start_cycle = start;};
  virtual void SetCurrentCycle(UInt_t cycle) {fCycle=cycle;};


  ClassDef(OnlineAnalysis,0);

};