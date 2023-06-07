//
// Created by Christiane Rahbek on 07-06-2023.
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

void FrontPageFig() {
    string filePath = "/mnt/d/IS659/finestructure_fit/analysis1/output/doubleDSSD/mergedRuns.root";
    string treeName = "a";

    string saveFilename, selectionCrit;
    saveFilename = "FrontPageFig.png";
    selectionCrit = "";

    int noOfBins = 1000;
    int minBin = 0;
    int maxBin = 8;

    auto canv = new TCanvas("","",1000,800);
    auto *f= new TFile(filePath.c_str());
    auto *tr=(TTree*)f->Get(treeName.c_str());
    auto hist = new TH2D("hist", "hist", noOfBins, minBin, maxBin, noOfBins, minBin, maxBin);
    tr->Draw("Edep0/1000:Edep1/1000 >> hist", selectionCrit.c_str(), "colz");


    gPad->SetTickx();
    gPad->SetTicky();

    hist->GetXaxis()->SetLabelOffset(999);
    hist->GetXaxis()->SetTickLength(0);
    hist->GetXaxis()->SetAxisColor(0);
    hist->GetYaxis()->SetLabelOffset(999);
    hist->GetYaxis()->SetTickLength(0);
    hist->GetYaxis()->SetAxisColor(0);
    hist->GetZaxis()->SetLabelOffset(999);
    hist->GetZaxis()->SetTickLength(0);
    hist->GetZaxis()->SetAxisColor(0);

    hist->SetTitle("");
    hist->SetStats(kFALSE);



    canv->SetLogz();




    canv->Update();
    canv->Draw();
    canv->SaveAs(saveFilename.c_str());
}