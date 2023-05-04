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

    string filePath = "/mnt/d/IS659/finestructure_fit/analysis1/output/doubleDSSD/Run121mliolio.root";
    string treeName = "a";
    string branchName = "TPROTONS";
    string selectionCrit = "hitAng>130 && ((id0==0 && id1==2) ||(id0==2 && id1==0) || (id0==1 && id1 ==3) || (id0==3 && id1==1)) && abs(FT0-FT1)<1500 && abs(Edep0-Edep1)<500";
    string title = "TPROTONS fit";
    string xLabel = "T [ms]";
    int noOfBins = 1000;
    int xMin = 0;
    int xMax = 3500;

    double BG_del = 10;
    double BG_open = 300;

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

    string formula = "[0]*(exp(-log(2)/119.1*(x+[3]))*(-1/(log(2)/119.1+log(2)/11.17)*(exp(-(log(2)/119.1+log(2)/11.17)*(x+[3]))-1)+1/(log(2)/119.1+log(2)/11.17+log(2)/55.12)*(exp(-(log(2)/119.1+log(2)/11.17+log(2)/55.12)*(x+[3]))-1)) + exp(-log(2)/839.9*(x+[3]))*(-1/(2*log(2)/119.1+log(2)/11.17+log(2)/55.12-log(2)/839.9)*(exp(-(2*log(2)/119.1+log(2)/11.17+log(2)/55.12-log(2)/839.9)*(x+[3]))-1)+1/(2*log(2)/119.1+log(2)/11.17-log(2)/839.9)*(exp(-(2*log(2)/119.1+log(2)/11.17-log(2)/839.9)*(x+[3]))-1)))+[1]*(x+[3])+[2]";
    //string formula = "[0]*(exp(-log(2)/119.1*(x+[3]))*(-1/(log(2)/119.1+log(2)/[4])*(exp(-(log(2)/119.1+log(2)/[4])*(x+[3]))-1)+1/(log(2)/119.1+log(2)/[4]+log(2)/[5])*(exp(-(log(2)/119.1+log(2)/[4]+log(2)/[5])*(x+[3]))-1)) + exp(-log(2)/839.9*(x+[3]))*(-1/(2*log(2)/119.1+log(2)/[4]+log(2)/[5]-log(2)/839.9)*(exp(-(2*log(2)/119.1+log(2)/[4]+log(2)/[5]-log(2)/839.9)*(x+[3]))-1)+1/(2*log(2)/119.1+log(2)/[4]-log(2)/839.9)*(exp(-(2*log(2)/119.1+log(2)/[4]-log(2)/839.9)*(x+[3]))-1)))+[1]*(x+[3])+[2]";


    auto fitFunc = new TF1("fitFunc", formula.c_str());
    fitFunc->FixParameter(3, -1.36e+01);
    fitFunc->SetParameters(100, 2.7e-01, -4.79,-1.36e+01,11.17, 55.12);
    fitFunc->SetParNames("N","a", "b", "t'", "t_f", "t_r");


    hist->Fit("fitFunc", "", "", BG_del, BG_del+BG_open);

    cout << "number of bins fitted = " << hist->GetXaxis()->FindBin(BG_open+BG_del) - hist->GetXaxis()->FindBin(xMin+BG_del) << endl;
}
