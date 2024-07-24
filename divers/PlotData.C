#include "TFile.h"
#include "TGraph.h"
#include "fstream"

void PlotData()
{   
    TGraph* graph = new TGraph();

    fstream data;

    data.open("data.dat", ios::in);

    double energy;
    double weight;

    while(1)
    {
        data >> energy >> weight;
        graph->AddPoint(energy,weight);
        if(data.eof()) break;
    }

    graph->Draw();
}