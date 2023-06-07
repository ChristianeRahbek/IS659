//
// Created by Christiane Rahbek on 02-06-2023.
//

#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPad.h"
#include "TBox.h"
#include "TLatex.h"

using namespace std;

void Na20Edep() {
    string filePath = "/mnt/d/IS659/finestructure_fit/calibrationAnalysis/output/Run208mlio.root";
    string treeName = "a";
    string branchName = "FE/1000";
    string selectionCrit = "";
    string title = "E_{dep} for the decay of ^{20}Na";
    string xLabel = "E_{dep} [MeV]";
    int noOfBins = 7000;
    double xMin = 0;
    double xMax = 7;

    auto canv = new TCanvas("", "", 1000, 800);

    auto *f= new TFile(filePath.c_str());
    auto *tr=(TTree*)f->Get(treeName.c_str());
    auto hist = new TH1F("hist", title.c_str(), noOfBins, xMin, xMax);
    tr->Draw((branchName +">> hist").c_str(), selectionCrit.c_str());
    gPad->SetTickx();
    gPad->SetTicky();

    hist ->SetXTitle(xLabel.c_str());
    hist->SetYTitle("Events/bin [MeV^{-1}]");
    gStyle ->SetOptStat(kFALSE);

    hist->GetXaxis()->SetLabelSize(0.04);
    hist->GetXaxis()->SetTitleSize(0.04);
    hist->GetYaxis()->SetLabelSize(0.04);
    hist->GetYaxis()->SetTitleSize(0.04);
    hist->GetYaxis()->SetMaxDigits(3);
    hist->GetYaxis()->SetTitleOffset(1.3);
    hist->GetXaxis()->SetTitleOffset(0.95);

    canv->Update();
    canv->SaveAs("Na20Edep.pdf");
}
