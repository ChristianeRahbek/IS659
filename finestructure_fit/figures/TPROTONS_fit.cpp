//
// Created by Christiane Rahbek on 31-03-2023.
//

#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TStyle.h"

using namespace std;

void TPROTONS_fit() {

    string filePath = "/mnt/d/IS659/finestructure_fit/analysis1/output/doubleDSSD/mergedRuns.root";
    string treeName = "a";
    string branchName = "TPROTONS";
    string selectionCrit = "hitAng>130 && ((id0==0 && id1==2) ||(id0==2 && id1==0) || (id0==1 && id1 ==3) || (id0==3 && id1==1)) && abs(FT0-FT1)<1500 && abs(Edep0-Edep1)<500";
    string title = "TPROTONS fit";
    string xLabel = "T [ms]";
    int noOfBins = 1000;
    int xMin = 0;
    int xMax = 3500;

    string saveFileName = "TPROTONS_fit.pdf";

    /* ********************** Making fit **************************** */

    auto canv = new TCanvas("", "", 1000, 800);
    //canv->SetLogy();

    /* Figure */
    auto *f= new TFile(filePath.c_str());
    auto *tr=(TTree*)f->Get(treeName.c_str());
    auto hist = new TH1F("hist", title.c_str(), noOfBins, xMin, xMax);
    tr->Draw((branchName +">> hist").c_str(), selectionCrit.c_str());

    double lambda8He = 1/119.1; // 1/halflife of 8He
    double lambda8Li = 1/839.9; // 1/halflife of 8Li

    //string B = to_string(lambda8He) + "/("+to_string(lambda8Li)+"-"+ to_string(lambda8He)+")";
    //string formula = "[0]*exp(-x*" + to_string(lambda8He) + ") + [1]*exp(-x*" + to_string(lambda8Li) + ")";
    string formula = "[0]*exp(-x*[1]) + [2]*exp(-x*[3])";


    auto fitFunc = new TF1("fitFunc", formula.c_str());

    hist->Fit("fitFunc", "", "", xMin+200, xMax-1200);
}
