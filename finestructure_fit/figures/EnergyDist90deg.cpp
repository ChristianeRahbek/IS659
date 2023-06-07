//
// Created by Christiane Rahbek on 16-05-2023.
//

#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"

using namespace std;

void EnergyDist90deg() {
    string filePath = "/mnt/d/IS659/finestructure_fit/analysis1/output/doubleDSSD/mergedRuns.root";
    string treeName = "a";

    bool removeBetas = true;
    string saveFilename, selectionCrit;
    if(removeBetas) {
        saveFilename = "energyDist90degNoBetas.pdf";
        selectionCrit = "((id0 == 0 && (id1 == 3 || id1 == 1)) || (id0 == 1 && (id1 == 0 || id1 == 2)) || (id0 == 2 && (id1 == 1 || id1 == 3)) || (id0 == 3 && (id1 == 2 || id1 == 0))) && abs(FT0-FT1)<2500 && isBetas == 0 && isParticleBeta == 0";

    } else {
        saveFilename = "energyDist90deg.pdf";
        selectionCrit = "((id0 == 0 && (id1 == 3 || id1 == 1)) || (id0 == 1 && (id1 == 0 || id1 == 2)) || (id0 == 2 && (id1 == 1 || id1 == 3)) || (id0 == 3 && (id1 == 2 || id1 == 0))) && abs(FT0-FT1)<2500";

    }
    int noOfBins = 1000;
    int minBin = 0;
    int maxBin = 8;

    auto canv = new TCanvas();
    auto *f= new TFile(filePath.c_str());
    auto *tr=(TTree*)f->Get(treeName.c_str());
    auto hist = new TH2D("hist", "hist", noOfBins, minBin, maxBin, noOfBins, minBin, maxBin);
    tr->Draw("Edep0/1000:Edep1/1000 >> hist", selectionCrit.c_str(), "colz");

    hist->SetXTitle("E1 [MeV]");
    hist->SetYTitle("E2 [MeV]");
    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetYaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetTitleSize(0.05);
    hist->SetTitle("");
    hist->SetStats(kFALSE);
    hist->GetXaxis()->SetTitleOffset(0.8);
    hist->GetYaxis()->SetTitleOffset(0.8);
    //canv->SetLogz();

    gPad->SetTickx();
    gPad->SetTicky();

    canv->Update();
    canv->Draw();
    canv->SaveAs(saveFilename.c_str());
}