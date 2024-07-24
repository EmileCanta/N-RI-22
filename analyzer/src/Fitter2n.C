#include "../include/Fitter84Ga2n.h"

void Fitter2n(const char* inputFile, Double_t min, Double_t max)
{
	TCanvas* c1 = new TCanvas();
	
	TFile* rootFile = TFile::Open(inputFile, "READ");
	
	TH1D *hist_tetra = (TH1D*)rootFile->Get("Aligned2n_tCond");

	TF1 *FitBat_tetra = new TF1("Bat_tetra", Bat_tetra, min, max, 2);

	c1->cd();
	hist_tetra->Fit("Bat_tetra", "R");
}