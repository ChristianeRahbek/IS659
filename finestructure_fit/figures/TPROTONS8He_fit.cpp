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
    int run = 121;

    string filePath, histName;
    double BG_del, BG_open;
    if (run==121) {
        filePath = "/mnt/d/IS659/tprotons_hists.root";
        histName = "h121";
        BG_del = 10;
        BG_open = 300;
    }
    else if(run == 163) {
        filePath = "/mnt/d/IS659/tprotons_hists.root";
        histName = "h163";
        BG_del = 15;
        BG_open = 300;
    }
    else {
        filePath = "/mnt/d/IS659/TPROTONS_id8_980.root";
        histName = "h";
        BG_del = 15;
        BG_open = 300;
    }

    string title = "TPROTONS fit";
    string xLabel = "T [ms]";
    int xMin = 0;
    int xMax = 3500;

    string saveFileName = "TPROTONS_fit.pdf";

    /* ********************** Making fit **************************** */

    auto canv = new TCanvas("", "", 1000, 800);
    //canv->SetLogy();

    /* Figure */
    auto *f= new TFile(filePath.c_str());
    auto *hist = (TH1F*)f->Get(histName.c_str());
    /*
    auto *tr=(TTree*)f->Get(treeName.c_str());
    auto hist = new TH1F("hist", title.c_str(), noOfBins, xMin, xMax);
    tr->Draw((branchName +">> hist").c_str(), selectionCrit.c_str());
     */

    hist->GetXaxis()->SetRangeUser(xMin,xMax);

    double lambda8He = 1/119.1; // 1/halflife of 8He
    double lambda8Li = 1/839.9; // 1/halflife of 8Li

    string formula;
    formula = "[0]*exp(-log(2)*(x+[3])/119.1)*(1/[1]*(1-exp(-[1]*(x+[3]))) + 1/([2]+[1])*(exp(-([2]+[1])*(x+[3]))-1)) + [4]*x + [5]";


    auto fitFunc = new TF1("fitFunc", formula.c_str(), xMin, xMax);
    fitFunc->SetParameters(100, log(2)/10, log(2)/20);
    fitFunc->SetParNames("N", "lambda_r", "lambda_f", "t'", "a", "b");
    fitFunc->SetParLimits(1,pow(10,-3),pow(10,3));
    fitFunc->SetParLimits(2,pow(10,-3),pow(10,3));

    hist->Fit("fitFunc", "V", "", xMin+BG_del, BG_open+BG_del);


    cout << "number of bins fitted = " << hist->GetXaxis()->FindBin(BG_open+BG_del) - hist->GetXaxis()->FindBin(xMin+BG_del) << endl;


    /* DRAWING THE FITTED FUNCTION WITHOUT THE BACKGROUND */
    auto fitfuncCanv = new TCanvas("", "", 1000, 800);

    fitfuncCanv->cd();

    string releaseform;
    releaseform = "[0]*exp(-log(2)*x/119.1)*(1/[1]*(1-exp(-[1]*x)) + 1/([2]+[1])*(exp(-([2]+[1])*x)-1))";


    auto releasefunc = new TF1("releaseFunc", releaseform.c_str(), xMin, xMax); //not really release func
    releasefunc->SetParameters(fitFunc->GetParameter(0), fitFunc->GetParameter(1), fitFunc->GetParameter(2));

    releasefunc->Draw();

    /* DRAWING THE FITTED RELEASE FUNCTION WITHOUT THE BACKGROUND */
    auto releasefuncCanv = new TCanvas("", "", 1000, 800);

    releasefuncCanv->cd();

    string relfuncform = "[0]*exp(-log(2)/119.1*x)*(1-exp(-log(1)/[1]*x))*exp(log(2)/[2]*t)";
    auto relfunc = new TF1("relFunc", relfuncform.c_str(), 0,1000);

    relfunc->SetParameters(1, log(2)/fitFunc->GetParameter(1), log(2)/fitFunc->GetParameter(2));

    relfunc->Draw();

    relfunc->SetParameters(1,10,20);

    relfunc->Draw("same");

}
