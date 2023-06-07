//
// Created by Christiane Rahbek on 16-05-2023.
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

void TPROTONS90deg() {
    string filePath = "/mnt/d/IS659/finestructure_fit/analysis1/output/doubleDSSD/mergedRuns.root";
    string treeName = "a";
    int BG_del = 10;
    int BG_open = 300;

    string title = "";

    string selectionCrit0 = "((id0 == 0 && (id1 == 3 || id1 == 1)) || (id0 == 1 && (id1 == 0 || id1 == 2)) || (id0 == 2 && (id1 == 1 || id1 == 3)) || (id0 == 3 && (id1 == 2 || id1 == 0))) && abs(Edep0-Edep1)>500 && abs(FT0-FT1)<2500";
    string selectionCrit1 = "((id0 == 0 && (id1 == 3 || id1 == 1)) || (id0 == 1 && (id1 == 0 || id1 == 2)) || (id0 == 2 && (id1 == 1 || id1 == 3)) || (id0 == 3 && (id1 == 2 || id1 == 0))) && abs(Edep0-Edep1)>500";

    int xMin = 0;
    int xMax = 3500;
    int noOfBins = 500;

    string saveFileName = "times90deg.pdf";
    string xTitle0 = "T_{p} [ms]";
    string xTitle1 = "#DeltaT [ns]";

    /* ************************************************** */

    auto canv = new TCanvas("", "", 1000, 400);
    canv->Divide(2,1);
    //canv->SetLogy();

    /* Figure */
    auto *f= new TFile(filePath.c_str());
    auto *tr0=(TTree*)f->Get(treeName.c_str());
    auto *tr1=(TTree*)f->Get(treeName.c_str());

    /* ***************************************** */
    canv->cd(1);
    auto hist0 = new TH1F("hist0", "T_{p} for reduced 90#circ coincidence data", noOfBins, xMin, xMax);
    tr0->Draw("TPROTONS >> hist0", selectionCrit0.c_str());

    hist0->GetXaxis()->SetRangeUser(xMin,xMax);

    int fitMin = 500;
    int fitMax = 3000;

    cout << "Fitting TPROTONS... " << endl;
    hist0->Fit("expo", "V", "", fitMin, fitMax);

    hist0->SetXTitle(xTitle0.c_str());
    hist0->SetYTitle("Entries/bin [ms^{-1}]");
    hist0->GetXaxis()->SetMaxDigits(3);
    hist0->GetYaxis()->SetMaxDigits(3);
    hist0->GetYaxis()->SetTitleOffset(1);
    hist0->GetXaxis()->SetTitleOffset(0.9);

    cout << "DOF = " << hist0->GetXaxis()->FindBin(fitMax) - hist0->GetXaxis()->FindBin(fitMin) - 2 << endl;

    gPad->SetTickx();
    gPad->SetTicky();

    hist0->SetStats(kFALSE);

    /* ***************************************** */
    canv->cd(2);
    auto hist1 = new TH1F("hist1", "Absolute time difference between hits", 500, 0, 50000);
    tr1->Draw("abs(FT0-FT1) >> hist1", selectionCrit1.c_str());

    hist1->SetXTitle(xTitle1.c_str());
    hist1->SetYTitle("Entries/bin [ns^{-1}]");
    hist1->GetXaxis()->SetMaxDigits(3);
    hist1->GetYaxis()->SetMaxDigits(3);
    hist1->GetYaxis()->SetTitleOffset(0.9);
    hist1->GetXaxis()->SetTitleOffset(0.9);

    gPad->SetTickx();
    gPad->SetTicky();

    hist1->SetStats(kFALSE);

    /* ***************************************** */
    hist0->SetTitleSize(0.06);
    hist1->SetTitleSize(0.06);

    hist0->GetYaxis()->SetTitleSize(0.05);
    hist0->GetXaxis()->SetTitleSize(0.05);
    hist0->GetYaxis()->SetLabelSize(0.05);
    hist0->GetXaxis()->SetLabelSize(0.05);
    hist1->GetYaxis()->SetTitleSize(0.05);
    hist1->GetXaxis()->SetTitleSize(0.05);
    hist1->GetYaxis()->SetLabelSize(0.05);
    hist1->GetXaxis()->SetLabelSize(0.05);



    canv->Update();
    canv->Draw();
    canv->SaveAs(saveFileName.c_str());
}