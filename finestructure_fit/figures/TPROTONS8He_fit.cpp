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

void TPROTONS8He_fit() {

    string filePath = "/mnt/d/IS659/TPROTONS_id8_980.root";
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
    auto *hist = (TH1F*)f->Get("h");
    /*
    auto *tr=(TTree*)f->Get(treeName.c_str());
    auto hist = new TH1F("hist", title.c_str(), noOfBins, xMin, xMax);
    tr->Draw((branchName +">> hist").c_str(), selectionCrit.c_str());
     */

    hist->GetXaxis()->SetRangeUser(xMin,xMax);

    double lambda8He = 1/119.1; // 1/halflife of 8He
    double lambda8Li = 1/839.9; // 1/halflife of 8Li

    //string formula = "[0]*exp(-x*" + to_string(lambda8He) + ") + [1]*exp(-x*" + to_string(lambda8Li) + ")";
    //string formula = "[0]*exp(-log(2)*x/119.1)*(1-exp(-log(2)*x/(0.010*pow(10,-3))))*(0.38*exp(-log(2)*x/(0.020*pow(10,-3)))+(1-0.38)*exp(-log(2)*x/(0.080*pow(10,-3))))";
    string formula = "[0]*exp(-log(2)*x/119.1)*(1-exp(-log(2)*x/(0.010*pow(10,3))))*(0.38*exp(-log(2)*x/(0.020*pow(10,3)))+(1-0.38)*exp(-log(2)*x/(0.080*pow(10,3))))";


    auto fitFunc = new TF1("fitFunc", formula.c_str());

    hist->Fit("fitFunc", "", "", xMin+20, 1000);
}
