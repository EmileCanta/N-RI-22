#include "TFile.h"
#include "TH1.h"
#include "fstream"

void RootToAscii(const char* inputfile, const char* outputfile)
{

TFile *data_file = TFile::Open(inputfile);

TH1D *hist = (TH1D*)data_file->Get("h1");

fstream output;

output.open(outputfile, ios::out);

for(int i = 0; i < 7000; i++)
{
output << i << " " << hist->GetBinContent(i) << endl;
}

output.close();

}