#include "Ga84RandomRate.h"

Double_t the_function(Double_t *x, Double_t *par)
{
	
	if(x[0] < t0)
    {
      	return 20*((Bgd_rate* N_cycles)+(pow((197.721/N_cycles), 2)*theta*exp(-theta*(197.721/N_cycles)))*N_cycles);
		
    }
    
  	if(x[0] >= t0 && x[0] < tc+t0)
    {	
      	return 
		20*((Bgd_rate* N_cycles)
		+(pow(((
			(197.721)
			+(-exp(-x[0]*l1)*(exp(t0*l1)-exp(x[0]*l1))*129.058)
			+((-exp(-x[0]*l1-x[0]*l2)*(-1+0.752108)*(exp(x[0]*l1+t0*l2)*l1-exp(x[0]*l1+x[0]*l2)*l1-exp(t0*l1+x[0]*l2)*l2+exp(x[0]*l1+x[0]*l2)*l2)*129.058)/(-l1+l2)))
			/N_cycles),2)*theta*exp(-theta*((197.721)
											+(-exp(-x[0]*l1)*(exp(t0*l1)-exp(x[0]*l1))*129.058)
											+((-exp(-x[0]*l1-x[0]*l2)*(-1+0.752108)*(exp(x[0]*l1+t0*l2)*l1-exp(x[0]*l1+x[0]*l2)*l1-exp(t0*l1+x[0]*l2)*l2+exp(x[0]*l1+x[0]*l2)*l2)*129.058)/(-l1+l2)))
											/N_cycles))*N_cycles);
    }
    
  	if(x[0] >= tc+t0 && x[0] <= ta)
    {
     	return
		20*((Bgd_rate* N_cycles)
		+(pow(((
			(197.721)
			+(-exp(-x[0]*l1)*(exp(t0*l1)-exp((t0+tc)*l1))*129.058)
			+((-exp(-x[0]*l1-x[0]*l2)*(-1+0.752108)*(-exp(x[0]*l1+t0*l2)*l1+exp(x[0]*l1+(t0+tc)*l2)*l1+exp(t0*l1+x[0]*l2)*l2-exp((t0+tc)*l1+x[0]*l2)*l2)*129.058)/(l1-l2)))
			/N_cycles),2)*theta*exp(-theta*((197.721)
											+(-exp(-x[0]*l1)*(exp(t0*l1)-exp((t0+tc)*l1))*129.058)
											+((-exp(-x[0]*l1-x[0]*l2)*(-1+0.752108)*(-exp(x[0]*l1+t0*l2)*l1+exp(x[0]*l1+(t0+tc)*l2)*l1+exp(t0*l1+x[0]*l2)*l2-exp((t0+tc)*l1+x[0]*l2)*l2)*129.058)/(l1-l2)))
											/N_cycles))*N_cycles);
    }
    
  	if(x[0] > ta)
    {
    	return 0;
    }
    
    return 0;
}

void Ga84RandomRate()
{
	TCanvas *c1 = new TCanvas();
	
	TFile *input = new TFile("/Users/cantacuzene/data/n-ri-22/runs/sorted_runs/84Ga/AllBut97.root");
	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0001);
	
	TH1D *hist_neutron = (TH1D*)input->Get("Aligned2n_tCond");
  	
  	TF1 *f1 = new TF1("the_function", the_function, 0.0e3, 3.3e3);
  	
  	hist_neutron->Draw();

	f1->SetLineColor(kOrange);
  	f1->Draw("SAME");
  	
  	TLegend *legend = new TLegend(0.65,0.65,0.80,0.85);
	legend->SetTextFont(72);
    legend->SetTextSize(0.02);
    legend->AddEntry(hist_neutron,"2n data","lpe");
    legend->AddEntry(f1,"Calculated Random rate","l");
    legend->Draw();  	
}