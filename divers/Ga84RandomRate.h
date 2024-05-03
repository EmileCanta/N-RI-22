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

//84Ga

Double_t t0 = 0.3e3;
Double_t tc = 2.0e3;
Double_t td = 1.0e3;
Double_t ta = tc+td+t0;

Double_t l1 = 0.0071458472; //84Ga
Double_t l2 = 7.358250324e-4; //84Ge
Double_t l3 = 3.746741517e-4; //83Ge
Double_t l4 = 1.724246718e-4; //84As

Double_t pn2 = 0.102;

Double_t theta = 0.2;
Double_t N_cycles = 720*29;
Double_t Bgd_rate = 3.2e-4;