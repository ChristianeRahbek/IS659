//
// Created by Christiane Rahbek on 02-06-2023.
//

#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"

using namespace std;

void TPROTONS8Li() {
    string filePath = "/mnt/d/IS659/finestructure_fit/analysis1/output/doubleDSSD/mergedRuns.root";
    string treeName = "a";
    int BG_del = 10;
    int BG_open = 300;

    string title = "";

    string selectionCrit0 = "((id0==0 && id1==2) || (id0==1 && id1==3) || (id0==2 && id1==0) || (id0==3 && id1==1)) && abs(FT0-FT1)<1500 && hitAng > 130 && abs(Edep0-Edep1)<500";
    int xMin = 0;
    int xMax = 3500;
    int noOfBins = 500;

    string saveFileName = "TPROTONS8Li.pdf";
    string xTitle0 = "T_{p} [ms]";

    /* ************************************************** */

    auto canv = new TCanvas("", "", 1000, 800);

    /* Figure */
    auto *f= new TFile(filePath.c_str());
    auto *tr0=(TTree*)f->Get(treeName.c_str());

    /* ***************************************** */
    auto hist0 = new TH1F("hist0", "T_{p} diagonal back-to-back data", noOfBins, xMin, xMax);
    tr0->Draw("TPROTONS >> hist0", selectionCrit0.c_str());

    hist0->GetXaxis()->SetRangeUser(xMin,xMax);

    int fitMin = 1000;
    int fitMax = 2100;

    cout << "Fitting TPROTONS... " << endl;
    hist0->Fit("expo", "V", "", fitMin, fitMax);

    hist0->SetXTitle(xTitle0.c_str());
    hist0->SetYTitle("Entries/bin [ms^{-1}]");
    hist0->GetXaxis()->SetMaxDigits(3);
    hist0->GetYaxis()->SetMaxDigits(4);
    hist0->GetYaxis()->SetTitleOffset(1);
    hist0->GetXaxis()->SetTitleOffset(0.95);
    hist0->GetXaxis()->SetTitleSize(0.045);
    hist0->GetXaxis()->SetLabelSize(0.045);
    hist0->GetYaxis()->SetTitleSize(0.045);
    hist0->GetYaxis()->SetLabelSize(0.045);

    cout << "DOF = " << hist0->GetXaxis()->FindBin(fitMax) - hist0->GetXaxis()->FindBin(fitMin) - 2 << endl;

    gPad->SetTickx();
    gPad->SetTicky();

    hist0->SetStats(kFALSE);

    canv->Update();
    canv->Draw();
    canv->SaveAs(saveFileName.c_str());
}