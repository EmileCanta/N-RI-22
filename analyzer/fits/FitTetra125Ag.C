#include "../include/Fitter125Ag.h"

Double_t Bateman_tot(Double_t *x, Double_t *par)
{
	
	if(x[0] < t0)
    {
      	return (41.8472);
    }
    
  	if(x[0] >= t0 && x[0] < tc+t0)
    {	
      	return (41.8472)
		+(-exp(-x[0]*l1)*(exp(t0*l1)-exp(x[0]*l1))*8.57509);
    }
    
  	if(x[0] >= tc+t0 && x[0] <= ta)
    {
     	return (41.8472)
		+(-exp(-x[0]*l1)*(exp(t0*l1)-exp((t0+tc)*l1))*8.57509);
    }
    
  	if(x[0] > ta)
    {
    	return 0;
    }
    
    return 0;
}

Double_t bgd(Double_t *x, Double_t *par)
{
	
	if(0 < x[0] && x[0] < ta)
    {
      	return (41.8472);
    }
    
    return 0;

}


Double_t Bateman_A1(Double_t *x, Double_t *par)
{  
	if(x[0] < t0)
    {
      	return 0;
    }
    
  	if(x[0] >= t0 && x[0] < tc+t0)
    {	
      	return (-exp(-x[0]*l1)*(exp(t0*l1)-exp(x[0]*l1))*8.57509);
    }
    
  	if(x[0] >= tc+t0 && x[0] <= ta)
    {
     	return (-exp(-x[0]*l1)*(exp(t0*l1)-exp((t0+tc)*l1))*8.57509);
    }
    
  	if(x[0] > ta)
    {
    	return 0;
    }
    
    return 0;
}

void FitTetra125Ag()
{
	TCanvas *c1 = new TCanvas();
	
	TFile *input = new TFile("/Users/cantacuzene/data/n-ri-22/runs/sorted_runs/125Ag/All.root");
	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0001);
	
	TH1D *hist_neutron = (TH1D*)input->Get("AlignedTetra_tSingle");
  	
  	TF1 *FitBatemanTot = new TF1("Bateman_tot", Bateman_tot, 0.0e3, 4.5e3);
  	TF1 *FitBgd = new TF1("bgd", bgd,  0.0e3, 4.5e3);
  	TF1 *FitA1 = new TF1("Bateman_A1", Bateman_A1,  0.0e3, 4.5e3);
  	
  	hist_neutron->Draw();

	FitBatemanTot->SetLineColor(kOrange);
  	FitBatemanTot->Draw("SAME");
  	
  	FitBgd->SetLineColor(kGreen);
  	FitBgd->Draw("SAME");

  	FitA1->SetLineColor(kCyan);
  	FitA1->Draw("SAME");
  	
  	TLegend *legend=new TLegend(0.65,0.65,0.80,0.85);
	legend->SetTextFont(72);
    legend->SetTextSize(0.02);
    legend->AddEntry(hist_neutron,"Data","lpe");
    legend->AddEntry(FitBatemanTot,"Bateman fit","l");
    legend->AddEntry(FitBgd,"Background","l");
    legend->AddEntry(FitA1,"Argent 125","l");
    legend->Draw();
  	
  	Double_t IntA1 = FitA1->Integral(0.0e3, 4.5e3);
  	Double_t IntBgd = FitBgd->Integral(0.0e3, 4.5e3);
  	
  	cout << "IntegralBgd:" << IntBgd << endl;
	
  	cout << "IntegralA1:" << IntA1 << endl;
}
