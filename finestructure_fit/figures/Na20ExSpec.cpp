//
// Created by Christiane Rahbek on 03-04-2023.
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

void Na20ExSpec(){
    string filePath = "/mnt/d/IS659/finestructure_fit/calibrationAnalysis/output/Run208mlio.root";
    string treeName = "a";
    string branchName = "Ex/1000";
    string selectionCrit = "";
    string title = "Excitation Energies of ^{20}Na";
    string xLabel = "Energy [MeV]";
    int noOfBins = 7000;
    double xMin = 4.500;
    double xMax = 11.500;

    auto canv = new TCanvas("", "", 1000, 800);

    auto *f= new TFile(filePath.c_str());
    auto *tr=(TTree*)f->Get(treeName.c_str());
    auto hist = new TH1F("hist", title.c_str(), noOfBins, xMin, xMax);
    tr->Draw((branchName +">> hist").c_str(), selectionCrit.c_str());

    hist ->SetXTitle(xLabel.c_str());
    hist->SetYTitle("Events");
    gStyle ->SetOptStat(kFALSE);

    TPad *p = new TPad("p", "p", .50, .40, 0.95, 0.92); //Where the overlaying canvas (pad) is placed (numbers between 0 and 1)
    p->Draw();
    p->cd();
    p->DrawFrame(6.80,0,8.0,500);

    auto hist01 = new TH1F("hist01", "hist01", noOfBins, xMin, xMax); //Making this so fit only shows up in the corner figure.
    tr->Draw((branchName +">> hist01").c_str(), selectionCrit.c_str(), "Same");
    hist01->Fit("gaus", "", "", 7.37, 7.48); //fitting
    hist01->GetYaxis()->SetMaxDigits(3);

    canv->Update();
    canv->SaveAs("Na20ExSpec.pdf");
}
