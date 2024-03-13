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
	TTree *raw_data_tree;

	Long64_t fEntries;
	Long64_t fEntry;

	UInt_t lastevent;

	UInt_t eventbefore;
	UInt_t eventafter;

	TFile *output_file;
	TTree *output_tree_single;
	TTree *output_tree_coinc;

	Double_t fGe_E;
	UChar_t fGe_Id;
	Double_t fGe_Time;

	UChar_t fBeta_Id;
	Double_t fBeta_Time;

	UChar_t fTetra_Id;
	Double_t fTetra_Time;

	std::vector<Double_t> fGe_E_single;
	std::vector<Double_t> fGe_Time_single;

	std::vector<Double_t> fBeta_Time_single;

	std::vector<Double_t> fTetra_Time_single;

	std::vector<Double_t> fTetra_Rings;

	std::vector<Double_t> fGeBeta_E_coinc;
	std::vector<Double_t> fGeTetra_E_coinc;

	std::vector<Double_t> fGeBeta_Time_coinc;
	std::vector<Double_t> fGeTetra_Time_coinc;

	std::vector<Double_t> fGeBeta_TimeDiff;
	std::vector<Double_t> fGeTetra_TimeDiff;

	std::vector<UInt_t> fGeBeta_Index;
	std::vector<UInt_t> fGeTetra_Index;

	std::vector<Double_t> fBeta1n_Time_coinc;
	std::vector<Double_t> fBeta1n_TimeDiff;
	std::vector<Double_t> fBeta1nBackward_Time_coinc;
	std::vector<Double_t> fBeta1nBackward_TimeDiff;

	std::vector<Double_t> fBeta2n_Time_coinc;
	std::vector<Double_t> fBeta2n_TimeDiffFirst;
	std::vector<Double_t> fBeta2n_TimeDiffSecond;
	std::vector<Double_t> fBeta2nBackward_Time_coinc;
	std::vector<Double_t> fBeta2nBackward_TimeDiffFirst;
	std::vector<Double_t> fBeta2nBackward_TimeDiffSecond;

	std::vector<Double_t> fBeta3n_Time_coinc;

	std::vector<Double_t> fStoring1n_Time;
	std::vector<Double_t> fStoring1n_TimeDiff;

	std::vector<Double_t> fStoring2n_Time;
	std::vector<Double_t> fStoring2n_TimeDiff;

	std::vector<Double_t> fStoring3n_Time;
	std::vector<Double_t> fStoring3n_TimeDiff;

	std::vector<Double_t> fnn_TimeDiff;
	std::vector<Double_t> fnn_Time_coinc;

	std::vector<Double_t> fStoringFirstNeutronCellGroup;
	std::vector<Double_t> fStoringSecondNeutronCellGroup;

	std::vector<Double_t> fFirstNeutronCellGroup;
	std::vector<Double_t> fSecondNeutronCellGroup;

	std::vector<Double_t> fTetra_Cycle;
	std::vector<Double_t> fBeta_Cycle;
	std::vector<Double_t> fGe_Cycle;

	UInt_t raw_energy;
	Double_t raw_time;
	UChar_t raw_det_nbr;
	UInt_t raw_coding;

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
	virtual void FillTetraCoincBranches(Double_t);
	virtual void FillBetaCoincBranches(Double_t);
	virtual void FillBetaXnCoincBranches(Double_t);
	virtual void FillBetaXnBackwardCoincBranches(Double_t);
	virtual void FillnnCoincBranches(Double_t);
	virtual Double_t Ge_alignement(UInt_t);
	virtual void ResetVar();
	virtual void ClearVectors();
	virtual void SetVar(UChar_t);
	virtual void Histogrammer(const char*);
	virtual UChar_t GetStatus(Double_t time);
};
