#include "Fitter84Ga.h"

Double_t Bateman_tot(Double_t *x, Double_t *par)
{
	
	if(x[0] < t0)
    {
      	return (157.81);
    }
    
  	if(x[0] >= t0 && x[0] < tc+t0)
    {	
      	return (157.81)
		+(-exp(-x[0]*l1)*(exp(t0*l1)-exp(x[0]*l1))*105.052)
		+((-exp(-x[0]*l1-x[0]*l2)*(-1+0.742436)*(exp(x[0]*l1+t0*l2)*l1-exp(x[0]*l1+x[0]*l2)*l1-exp(t0*l1+x[0]*l2)*l2+exp(x[0]*l1+x[0]*l2)*l2)*105.052)/(-l1+l2));
    }
    
  	if(x[0] >= tc+t0 && x[0] <= ta)
    {
     	return (157.81)
		+(-exp(-x[0]*l1)*(exp(t0*l1)-exp((t0+tc)*l1))*105.052)
		+((-exp(-x[0]*l1-x[0]*l2)*(-1+0.742436)*(-exp(x[0]*l1+t0*l2)*l1+exp(x[0]*l1+(t0+tc)*l2)*l1+exp(t0*l1+x[0]*l2)*l2-exp((t0+tc)*l1+x[0]*l2)*l2)*105.052)/(l1-l2));
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
      	return (157.81);
    }
    
    return 0;

}


Double_t Bateman_A1(Double_t *x, Double_t *par)
{  
	if(x[0] < t0)
    {
      	return (157.81);
    }
    
  	if(x[0] >= t0 && x[0] < tc+t0)
    {	
      	return (157.81)+(-exp(-x[0]*l1)*(exp(t0*l1)-exp(x[0]*l1))*105.052);
    }
    
  	if(x[0] >= tc+t0 && x[0] <= ta)
    {
     	return (157.81)+(-exp(-x[0]*l1)*(exp(t0*l1)-exp((t0+tc)*l1))*105.052);
    }
    
  	if(x[0] > ta)
    {
    	return 0;
    }
    
    return 0;
}

Double_t Bateman_A2(Double_t *x, Double_t *par)
{  
	if(x[0] < t0)
    {
      	return (157.81);
    }
    
  	if(x[0] >= t0 && x[0] < tc+t0)
    {	
      	return (157.81)+((-exp(-x[0]*l1-x[0]*l2)*(-1+0.742436)*(exp(x[0]*l1+t0*l2)*l1-exp(x[0]*l1+x[0]*l2)*l1-exp(t0*l1+x[0]*l2)*l2+exp(x[0]*l1+x[0]*l2)*l2)*105.052)/(-l1+l2));
    }
    
  	if(x[0] >= tc+t0 && x[0] <= ta)
    {
     	return (157.81)+((-exp(-x[0]*l1-x[0]*l2)*(-1+0.742436)*(-exp(x[0]*l1+t0*l2)*l1+exp(x[0]*l1+(t0+tc)*l2)*l1+exp(t0*l1+x[0]*l2)*l2-exp((t0+tc)*l1+x[0]*l2)*l2)*105.052)/(l1-l2));
    }
    
  	if(x[0] > ta)
    {
    	return 0;
    }
    
    return 0;
}

void FitTetra84Ga()
{
	TCanvas *c1 = new TCanvas();
	
	TFile *input = new TFile("/Users/cantacuzene/n-ri-22/runs/sorted_runs/84Ga/AllRuns.root");
	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0001);
	
	TH1D *hist_neutron = (TH1D*)input->Get("AlignedTetra_Time_single");
  	
  	TF1 *FitBatemanTot = new TF1("Bateman_tot", Bateman_tot, 0.0e3, 3.3e3);
  	TF1 *FitBgd = new TF1("bgd", bgd,  0.0e3, 3.3e3);
  	TF1 *FitA1 = new TF1("Bateman_A1", Bateman_A1,  0.0e3, 3.3e3);
	TF1 *FitA2 = new TF1("Bateman_A2", Bateman_A2,  0.0e3, 3.3e3);
  	
  	hist_neutron->Draw();

	FitBatemanTot->SetLineColor(kOrange);
  	FitBatemanTot->Draw("SAME");
  	
  	FitBgd->SetLineColor(kGreen);
  	FitBgd->Draw("SAME");

  	FitA1->SetLineColor(kCyan);
  	FitA1->Draw("SAME");

	FitA2->SetLineColor(kBlack);
  	FitA2->Draw("SAME");
  	
  	TLegend *legend=new TLegend(0.65,0.65,0.80,0.85);
	legend->SetTextFont(72);
    legend->SetTextSize(0.02);
    legend->AddEntry(hist_neutron,"Data","lpe");
    legend->AddEntry(FitBatemanTot,"Bateman fit","l");
    legend->AddEntry(FitBgd,"Background","l");
    legend->AddEntry(FitA1,"Galium 84","l");
	legend->AddEntry(FitA2,"Germanium 84","l");
    legend->Draw();
  	
  	Double_t IntA1 = FitA1->Integral(0.0e3, 3.3e3);
  	Double_t IntBgd = FitBgd->Integral(0.0e3, 3.3e3);
	Double_t IntA2 = FitA2->Integral(0.0e3, 3.3e3);
  	
  	cout << "IntegralBgd:" << IntBgd << endl;
	
  	cout << "IntegralA1-Bgd:" << (IntA1 - IntBgd) << endl;
	cout << "IntegralA2-Bgd:" << (IntA2 - IntBgd) << endl;
	cout << "Total integral with background:" << ((IntA1 - IntBgd) + (IntA2 - IntBgd) + IntBgd) << endl;
	cout << "Total integral without background:" << ((IntA1 - IntBgd) + (IntA2 - IntBgd)) << endl;
  	
}
