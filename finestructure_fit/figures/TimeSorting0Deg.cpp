//
// Created by Christiane Rahbek on 10-05-2023.
//

#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPad.h"
#include "TArrow.h"
#include <TROOT.h>
#include "TColor.h"
#include "TBox.h"
#include "TLatex.h"

using namespace std;

void TimeSorting0Deg() {
/* INSERT INFORMATION HERE*/
    string filePath = "/mnt/d/IS659/finestructure_fit/analysis1/output/doubleDSSD/mergedRuns.root";
    string treeName = "a";
    string branchName = "abs(FT0-FT1)";
    string selectionCrit = "id0==id1";
    string title = "";
    string xLabel = "#DeltaT [ns]";
    int noOfBins = 500;
    int xMin = 0;
    int xMax = 50000;

    string saveFileName = "TimeSorting0Deg.pdf";

/* ********************** MAKING FIGURE ************************************************** */

    auto canv = new TCanvas("", "", 1500, 800);
    canv->SetLogy();

/* First figure */
    auto *f= new TFile(filePath.c_str());
    auto *tr=(TTree*)f->Get(treeName.c_str());
    auto hist = new TH1F("hist", title.c_str(), noOfBins, xMin, xMax);
    tr->Draw((branchName +">> hist").c_str(), selectionCrit.c_str());

    hist->SetXTitle(xLabel.c_str());
    hist->SetYTitle("Events/bin [ns^{-1}]");
    canv->SetLeftMargin(0.1);
    hist->GetYaxis()->SetTitleOffset(1);
    hist->GetYaxis()->SetMaxDigits(3); //Setting numbers on axis to x*10^y
    hist->GetXaxis()->SetMaxDigits(3); //Setting numbers on axis to x*10^y
    gStyle ->SetOptStat(kFALSE);
    gStyle ->SetPalette(55);

    gPad->SetTickx();
    gPad->SetTicky();

    hist->GetYaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetTitleOffset(0.95);
    hist->GetYaxis()->SetTitleOffset(0.8);

/*
//Making zoom on the graph in the corner
    TPad *p = new TPad("p", "p", .47, .5, 0.97, 0.95); //Where the overlaying canvas (pad) is placed (numbers between 0 and 1)
    p->SetLogy();
    p->Draw();
    p->cd();
    //p->SetLogy();

    auto hist01 = new TH1F("hist01", "hist01", noOfBins, xMin, xMax); //Making this so fit only shows up in the corner figure.
    tr->Draw((branchName +">> hist01").c_str(), selectionCrit.c_str());


    hist01->GetXaxis()->SetRangeUser(0,5000);
    hist01->GetXaxis()->SetMaxDigits(3);
    hist01->GetYaxis()->SetMaxDigits(3);
    hist01->SetTitle("");
    hist01->Draw();

    gPad->SetTickx();
    gPad->SetTicky();


    canv->cd();
*/

    //auto BGpad = new TPad("BGpad", "BGpad", 0,0,1,1);
    //BGpad->cd();

    auto *BGarrow = new TArrow(33000,1700,32000,400, 0.01,"|>");
    BGarrow->Draw();
    auto BGtext = new TLatex(33027.5,1710.34,"Background");
    BGtext->SetTextSize(0.035);
    BGtext->Draw();

    auto *CoinArrow = new TArrow(4250,15942,1100,13000, 0.01,"|>");
    CoinArrow->Draw();
    auto *CoinArrow1 = new TArrow(4250,15942,2500,2000, 0.01,"|>");
    CoinArrow1->Draw();
    auto *CoinArrow2 = new TArrow(4250,15942,13000,3000, 0.01,"|>");
    CoinArrow2->Draw();
    auto CoinText = new TLatex(4306,15942,"Valid coincidences");
    CoinText->SetTextSize(0.035);
    CoinText->Draw();


    canv->Modified();
    canv->Update();
    canv->Draw();

    canv->SaveAs(saveFileName.c_str());
}