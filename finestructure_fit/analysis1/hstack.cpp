//
//
//

#include "THStack.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"

void hstack() {
    TFile *f=new TFile("/mnt/d/IS659/finestructure_fit/analysis1/output/bat/Run167mlio.root");
    TTree *tr=(TTree*)f->Get("a");

    auto hs = new THStack("hs","Stacked 1D histograms");
    //create three 1-d histograms
    //auto h1st = new TH1F("h1st","test hstack",100,-4,4);

    auto h1st = new TH1F("h1st", "h1st", 1000, 0 , 8000);
    tr -> Draw("Edssd >> h1st");
    h1st ->SetLineColor(kRed);
    h1st->SetMarkerStyle(21);
    h1st->SetMarkerColor(kRed);
    hs->Add(h1st);

    auto h2st = new TH1F("h2st", "h2st", 1000, 0 , 8000);
    tr -> Draw("Edep0 >> h2st");
    h2st ->SetLineColor(kBlue);
    h2st->SetMarkerStyle(21);
    h2st->SetMarkerColor(kBlue);
    hs->Add(h2st);

    auto cst = new TCanvas("cst","stacked hists",10,10,700,700);

    hs->Draw();

    hs->SaveAs("/mnt/d/IS659/finestructure_fit/analysis1/output/bat/hists.png");
}