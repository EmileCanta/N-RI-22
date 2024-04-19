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

using namespace std;

class Sorter{
  
	private:

	TFile *data_file;
	TFile *output_file;

	TTree *raw_data_tree;
	TTree *output_tree_single;
	TTree *output_tree_coinc;

	UInt_t raw_energy;
	Double_t raw_time;
	UChar_t raw_det_nbr;
	UInt_t raw_coding;

	Long64_t fEntries;
	Long64_t fEntry;

	UInt_t lastevent;
	UInt_t neutcount;
	UInt_t eventafter;

	int neut_id;
	int lastneutron;

	Double_t fGe_E;
	UChar_t fGe_Id;
	Double_t fGe_Time;

	UChar_t fBeta_Id;
	Double_t fBeta_Time;

	UChar_t fTetra_Id;
	Double_t fTetra_Time;

		//Single

	//Time
	std::vector<Double_t> fGe_tSingle;
	std::vector<Double_t> fBeta_tSingle;
	std::vector<Double_t> fTetra_tSingle;

	//Energy
	std::vector<Double_t> fGe_ESingle;

	//Cycles
	std::vector<UInt_t> fTetra_Cycle;
	std::vector<UInt_t> fBeta_Cycle;
	std::vector<UInt_t> fGe_Cycle;

	//Tetra cells and rings
	std::vector<UInt_t> fTetra_Rings;
	std::vector<UInt_t> fTetra_CellGroups;

		//Coincidences

	//Ge-Beta
	std::vector<Double_t> fGeBeta_ECond;
	std::vector<Double_t> fGeBeta_tCond;
	std::vector<Double_t> fGeBeta_tDiff;

	//Ge-Tetra
	std::vector<Double_t> fGeTetra_ECond;
	std::vector<Double_t> fGeTetra_tCond;
	std::vector<Double_t> fGeTetra_tDiff;

	//Beta-Tetra
	std::vector<Double_t> fBeta1n_tCond;
	std::vector<Double_t> fBeta2n_tCond;
	std::vector<Double_t> fBeta1nBackward_tCond;
	std::vector<Double_t> fBeta2nBackward_tCond;
	std::vector<Double_t> fBeta1n_tDiff;
	std::vector<Double_t> fBeta2n_tDiffFirst;
	std::vector<Double_t> fBeta2n_tDiffSecond;
	std::vector<Double_t> fBeta1nBackward_tDiff;
	std::vector<Double_t> fBeta2nBackward_tDiffFirst;
	std::vector<Double_t> fBeta2nBackward_tDiffSecond;

	//Tetra-Tetra
	std::vector<Double_t> fFirstNeutronCellGroup;
	std::vector<Double_t> fSecondNeutronCellGroup;
	std::vector<Double_t> fStoringFirstNeutronCellGroup;
	std::vector<Double_t> fStoringSecondNeutronCellGroup;

	//Storing vectors
	std::vector<Double_t> fStoring1n_Time;
	std::vector<Double_t> fStoring1n_TimeDiff;
	std::vector<Double_t> fStoring2n_Time;
	std::vector<Double_t> fStoring2n_TimeDiff;
	std::vector<Double_t> fSecondNeut_tDiff;
	std::vector<Double_t> fSecondNeut_tCond;

    TH1D* Ge_tSingle = new TH1D("Ge_tSingle", "Ge_tSingle", 30000,0,30000);
    TH1D* Beta_tSingle = new TH1D("Beta_tSingle", "Beta_tSingle", 30000,0,30000);
    TH1D* Tetra_tSingle = new TH1D("Tetra_tSingle", "Tetra_tSingle", 30000,0,30000);

	TH1D* Ge_ESingle = new TH1D("Ge_ESingle", "Ge_ESingle", 7000,0,7000);

    TH1I* Tetra_Cycle = new TH1I("Tetra_Cycle", "Tetra_Cycle", 1000,0,1000);
    TH1I* Beta_Cycle = new TH1I("Beta_Cycle", "Beta_Cycle", 1000,0,1000);
    TH1I* Ge_Cycle = new TH1I("Ge_Cycle", "Ge_Cycle", 1000,0,1000);

	TH1I* Tetra_Rings = new TH1I("Tetra_Rings", "Tetra_Rings", 28, 0, 14);
    TH1I* Tetra_CellGroups = new TH1I("Tetra_CellGroups", "Tetra_CellGroups", 28, 0, 14);

	TH1D* GeBeta_ECond = new TH1D("GeBeta_ECond", "GeBeta_ECond", 7000, 0, 7000);
	TH1D* GeBeta_tCond = new TH1D("GeBeta_tCond","GeBeta_tCond", 30000, 0, 30000);
    TH1D* GeBeta_tDiff = new TH1D("GeBeta_tDiff","GeBeta_tDiff", 5000, 0, 10);

	TH1D* GeTetra_ECond = new TH1D("GeTetra_ECond", "GeTetra_ECond", 7000, 0, 7000);
	TH1D* GeTetra_tCond = new TH1D("GeTetra_tCond","GeTetra_tCond", 30000, 0, 30000);
    TH1D* GeTetra_tDiff = new TH1D("GeTetra_tDiff", "GeTetra_tDiff", 1000, 0, 200);

	TH1D* Beta1n_tCond = new TH1D("Beta1n_tCond","Beta1n_tCond", 30000, 0, 30000);
    TH1D* Beta2n_tCond = new TH1D("Beta2n_tCond","Beta2n_tCond", 30000, 0, 30000);
	TH1D* Beta1nBackward_tCond = new TH1D("Beta1nBackward_tCond","Beta1nBackward_tCond", 30000, 0, 30000);
    TH1D* Beta2nBackward_tCond = new TH1D("Beta2nBackward_tCond","Beta2nBackward_tCond", 30000, 0, 30000);

    TH1D* Beta1n_tDiff = new TH1D("Beta1n_tDiff", "Beta1n_tDiff", 20000, -2000, 2000);
    TH1D* Beta2n_tDiffFirst = new TH1D("Beta2n_tDiffFirst", "Beta2n_tDiffFirst", 20000, -2000, 2000);
    TH1D* Beta2n_tDiffSecond = new TH1D("Beta2n_tDiffSecond", "Beta2n_tDiffSecond", 20000, -2000, 2000);
    TH1D* Beta1nBackward_tDiff = new TH1D("Beta1nBackward_tDiff", "Beta1nBackward_tDiff", 20000, -2000, 2000);
    TH1D* Beta2nBackward_tDiffFirst = new TH1D("Beta2nBackward_tDiffFirst", "Beta2nBackward_tDiffFirst", 20000, -2000, 2000);
    TH1D* Beta2nBackward_tDiffSecond = new TH1D("Beta2nBackward_tDiffSecond", "Beta2nBackward_tDiffSecond", 20000, -2000, 2000);

    TH1D* SecondNeut_tDiff = new TH1D("SecondNeut_tDiff", "SecondNeut_tDiff", 10000, 0, 2000);
    TH1D* SecondNeut_tCond = new TH1D("SecondNeut_tCond", "SecondNeut_tCond", 30000, 0, 30000);
    TH1D* FirstNeutronCellGroup = new TH1D("FirstNeutronCellGroup", "FirstNeutronCellGroup", 28, 0, 14);
    TH1D* SecondNeutronCellGroup = new TH1D("SecondNeutronCellGroup", "SecondNeutronCellGroup", 28, 0, 14);

    //Bidim
    TH2D* ESvsBTD = new TH2D("ESvsBTD", "ESvsBTD", 5000, 0, 10, 7000, 0, 7000);
    TH2D* ESvsTTD = new TH2D("ESvsTTD", "ESvsTTD", 1000, 0, 200, 7000, 0, 7000);
    TH2D* CellGroups = new TH2D("CellGroups", "CellGroups", 28, 0, 14, 28, 0, 14);

	protected:

	UInt_t fCoding;
	UChar_t fStatus;
	Double_t coll_time;
	Double_t acqu_time;
	Double_t bgd_time;

	public:

	Sorter(const char*, const char*, Double_t, Double_t, Double_t);
	virtual ~Sorter();

	void Load_rawdata(const char*);

	protected:

	virtual void SetTreesAndBranches(const char*);
	virtual void FillSingleBranches();
	virtual void FillNeutronGammaCoincBranches(Double_t);
	virtual void FillBetaGammaCoincBranches(Double_t);
	virtual void FillBetaNeutronCoincBranches(Double_t);
	virtual void FillBetaNeutronBackwardCoincBranches(Double_t);
	virtual void FillNeutronNeutronCoincBranches(Double_t);
	virtual Double_t Ge_alignement(UInt_t);
	virtual void ResetVar();
	virtual void ClearVectors();
	virtual void SetVar(UChar_t);
	//virtual void Histogrammer(const char*);
	virtual UChar_t GetStatus(Double_t time);
};
