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
#include <root/TF1.h>
#include <root/TCanvas.h>

using namespace std;

//125Ag

Double_t t0 = 0.5e3;
Double_t tc = 2.0e3;
Double_t td = 2.0e3;
Double_t ta = tc+td+t0;

Double_t l1 = 0.003938; //Production Argent non iso
Double_t l2 = 2.937e-4; //Production Indium non iso


//par0 is Bgd
//par1 is Phi1
//par2 is Phi2

Double_t Bat_tetra(Double_t *x_2, Double_t *par_2)
{	
    
	if(x_2[0] <= t0)
    {
      	return (par_2[0]);
    }
    
  	if(x_2[0] >= t0 && x_2[0] <= tc+t0)
    {  	
      	return (par_2[0])
		+(-exp(-x_2[0]*l1)*(exp(t0*l1)-exp(x_2[0]*l1))*par_2[1]);
    }
    
  	if(x_2[0] >= tc+t0 && x_2[0] <= ta)
    {
     	return (par_2[0])
		+(-exp(-x_2[0]*l1)*(exp(t0*l1)-exp((t0+tc)*l1))*par_2[1]);
    }
    
  	if(x_2[0] > ta)
    {
    	return 0;
    }
    
    return 0;
}

Double_t Bat_beta(Double_t *x_1, Double_t *par_1)
{
	
	if(x_1[0] <= t0)
    {
      	return (par_1[0]);
    }
    
  	if(x_1[0] >= t0 && x_1[0] <= tc+t0)
    {	
      	return (par_1[0])
		+(-exp(-x_1[0]*l1)*(exp(t0*l1)-exp(x_1[0]*l1))*par_1[1])
		+(-exp(-x_1[0]*l2)*(exp(t0*l2)-exp(x_1[0]*l2))*par_1[2]); 
    }
    
  	if(x_1[0] >= tc+t0 && x_1[0] <= ta)
    {
     	return (par_1[0])
		+(-exp(-x_1[0]*l1)*(exp(t0*l1)-exp((t0+tc)*l1))*par_1[1])
		+(-exp(-x_1[0]*l2)*(exp(t0*l2)-exp((t0+tc)*l2))*par_1[2]);
    }
    
  	if(x_1[0] > ta)
    {
    	return 0;
    }
    
    return 0;
}