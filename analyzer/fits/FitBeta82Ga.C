#include "Fitter82Ga.h"

Double_t Bateman_tot(Double_t *x, Double_t *par)
{
	
	if(x[0] < t0)
    {
      	return (129.287);
    }
    
  	if(x[0] >= t0 && x[0] < tc+t0)
    {	
      	return (129.287)
		+(-exp(-x[0]*l1)*(exp(t0*l1)-exp(x[0]*l1))*1223.22)
		+((-exp(-x[0]*l1-x[0]*l2)*(-1+0.717522)*(exp(x[0]*l1+t0*l2)*l1-exp(x[0]*l1+x[0]*l2)*l1-exp(t0*l1+x[0]*l2)*l2+exp(x[0]*l1+x[0]*l2)*l2)*1223.22)/(-l1+l2))
		+((exp(-x[0]*l1-x[0]*l3)*(0.717522)*(exp(x[0]*l1+t0*l3)*l1-exp(x[0]*l1+x[0]*l3)*l1-exp(t0*l1+x[0]*l3)*l3+exp(x[0]*l1+x[0]*l3)*l3)*1223.22)/(-l1+l3));
    }
    
  	if(x[0] >= tc+t0 && x[0] <= ta)
    {
     	return (129.287)
		+(-exp(-x[0]*l1)*(exp(t0*l1)-exp((t0+tc)*l1))*1223.22)
		+((-exp(-x[0]*l1-x[0]*l2)*(-1+0.717522)*(-exp(x[0]*l1+t0*l2)*l1+exp(x[0]*l1+(t0+tc)*l2)*l1+exp(t0*l1+x[0]*l2)*l2-exp((t0+tc)*l1+x[0]*l2)*l2)*1223.22)/(l1-l2))
		+((exp(-x[0]*l1-x[0]*l3)*(0.717522)*(-exp(x[0]*l1+t0*l3)*l1+exp(x[0]*l1+(t0+tc)*l3)*l1+exp(t0*l1+x[0]*l3)*l3-exp((t0+tc)*l1+x[0]*l3)*l3)*1223.22)/(l1-l3));
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
      	return (129.287);
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
      	return (-exp(-x[0]*l1)*(exp(t0*l1)-exp(x[0]*l1))*1223.22);
    }
    
  	if(x[0] >= tc+t0 && x[0] <= ta)
    {
     	return (-exp(-x[0]*l1)*(exp(t0*l1)-exp((t0+tc)*l1))*1223.22);
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
      	return 0;
    }
    
  	if(x[0] >= t0 && x[0] < tc+t0)
    {	
      	return ((-exp(-x[0]*l1-x[0]*l2)*(-1+0.717522)*(exp(x[0]*l1+t0*l2)*l1-exp(x[0]*l1+x[0]*l2)*l1-exp(t0*l1+x[0]*l2)*l2+exp(x[0]*l1+x[0]*l2)*l2)*1223.22)/(-l1+l2));
    }
    
  	if(x[0] >= tc+t0 && x[0] <= ta)
    {	
     	return ((-exp(-x[0]*l1-x[0]*l2)*(-1+0.717522)*(-exp(x[0]*l1+t0*l2)*l1+exp(x[0]*l1+(t0+tc)*l2)*l1+exp(t0*l1+x[0]*l2)*l2-exp((t0+tc)*l1+x[0]*l2)*l2)*1223.22)/(l1-l2));
    }
    
  	if(x[0] > ta)
    {
    	return 0;
    }
    
    return 0;
}

Double_t Bateman_A3(Double_t *x, Double_t *par)
{		
	if(x[0] < t0)
    {
      	return 0;
    }
    
  	if(x[0] >= t0 && x[0] < tc+t0)
    {	
      	return ((exp(-x[0]*l1-x[0]*l3)*(0.717522)*(exp(x[0]*l1+t0*l3)*l1-exp(x[0]*l1+x[0]*l3)*l1-exp(t0*l1+x[0]*l3)*l3+exp(x[0]*l1+x[0]*l3)*l3)*1223.22)/(-l1+l3));
    }
    
  	if(x[0] >= tc+t0 && x[0] <= ta)
    {	
     	return ((exp(-x[0]*l1-x[0]*l3)*(0.717522)*(-exp(x[0]*l1+t0*l3)*l1+exp(x[0]*l1+(t0+tc)*l3)*l1+exp(t0*l1+x[0]*l3)*l3-exp((t0+tc)*l1+x[0]*l3)*l3)*1223.22)/(l1-l3));
    }
    
  	if(x[0] > ta)
    {
    	return 0;
    }
    
    return 0;
}

void FitBeta82Ga()
{
	TCanvas *c1 = new TCanvas();
	
	TFile *input = new TFile("/Users/cantacuzene/data/n-ri-22/runs/sorted_runs/RUN121.root");
	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0001);
	
	TH1D *hist_beta = (TH1D*)input->Get("AlignedBeta_Time_single");
  	
  	TF1 *FitBatemanTot = new TF1("Bateman_tot", Bateman_tot, 0.0e3, 6.5e3);
  	TF1 *FitBgd = new TF1("bgd", bgd, 0.0e3, 6.5e3);
  	TF1 *FitA1 = new TF1("Bateman_A1", Bateman_A1, 0.0e3, 6.5e3);
  	TF1 *FitA2 = new TF1("Bateman_A2", Bateman_A2, 0.0e3, 6.5e3);
	TF1 *FitA3 = new TF1("Bateman_A3", Bateman_A3, 0.0e3, 6.5e3);

  	hist_beta->Draw();

	FitBatemanTot->SetLineColor(kOrange);
  	FitBatemanTot->Draw("SAME");
  	
  	FitBgd->SetLineColor(kGreen);
  	FitBgd->Draw("SAME");

  	FitA1->SetLineColor(kCyan);
  	FitA1->Draw("SAME");
  	
  	FitA2->SetLineColor(kBlack);
  	FitA2->Draw("SAME");

	FitA3->SetLineColor(kRed);
  	FitA3->Draw("SAME");
  	
  	TLegend *legend = new TLegend(0.65,0.65,0.80,0.85);
	legend->SetTextFont(72);
    legend->SetTextSize(0.02);
    legend->AddEntry(hist_beta,"Data","lpe");
    legend->AddEntry(FitBatemanTot,"Bateman fit","l");
    legend->AddEntry(FitBgd,"Background","l");
    legend->AddEntry(FitA1,"Galium 84","l");
    legend->AddEntry(FitA2,"Germanium 84","l");
	legend->AddEntry(FitA3,"Germanium 83","l");
    legend->Draw();
  	
  	Double_t IntA1 = FitA1->Integral(0.0e3, 6.5e3);
  	Double_t IntBgd = FitBgd->Integral(0.0e3, 6.5e3);
	Double_t IntA2 = FitA2->Integral(0.0e3, 6.5e3);
	Double_t IntA3 = FitA3->Integral(0.0e3, 6.5e3);
  	
	cout << "IntegralBgd:" << IntBgd << endl;
	
  	cout << "IntegralA1:" << IntA1 << endl;
	cout << "IntegralA2:" << IntA2 << endl;
	cout << "IntegralA3:" << IntA3 << endl;
  	
}
