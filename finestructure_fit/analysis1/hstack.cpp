//
//
//

#include "THStack.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"

void hstack() {
    auto *f1= new TFile("/mnt/d/IS659/finestructure_fit/analysis1/output/Run167mlio.root");
    auto *tr1=(TTree*)f1->Get("a");

    auto h1st = new TH1F("h1st", "h1st", 1000, 0 , 8000);
    tr1 -> Draw("Edssd >> h1st");
    h1st ->SetLineColor(kRed);

    auto *f2=new TFile("/mnt/d/IS659/finestructure_fit/analysis1/output/bat/Run167mlio.root");
    auto *tr2=(TTree*)f2->Get("a");

    auto h2st = new TH1F("h2st", "h2st", 1000, 0 , 8000);
    tr2 -> Draw("Edssd >> h2st");
    h2st ->SetLineColor(kBlue);

    auto hs = new THStack("hs","Stacked 1D histograms");
    hs->SetHistogram(h1st);
    hs->Add(h1st, "hist");
    hs->Add(h2st, "hist");

    hs->Draw();

    hs->SaveAs("/mnt/d/IS659/finestructure_fit/analysis1/output/bat/hists.png");
}