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

//82Ga

Double_t t0 = 0.5e3;
Double_t tc = 3.0e3;
Double_t td = 3.0e3;
Double_t ta = tc+td+t0;

Double_t l1 = 0.0011552453;
Double_t l2 = 1.732867951e-4;
Double_t l3 = 9.120357639e-5;

Double_t theta = 0.2;
Double_t N_cycles = 500;
Double_t Bgd_rate = 3.2e-4;