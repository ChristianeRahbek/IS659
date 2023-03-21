//
// Created by Christiane Rahbek on 15-03-2023.
// This script makes a figure that describes how i do my time-coincidence analysis
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

using namespace std;

void TimeSorting() {
/* INSERT INFORMATION HERE*/
    string filePath = "/mnt/d/IS659/finestructure_fit/analysis1/output/2coincidence/mergedRuns.root";
    string treeName = "a";
    string branchName = "abs(FT0-FT1)";
    string selectionCrit = "((id0==0 && id1==2) || (id0==1 && id1==3) || (id0==2 && id1==0) || (id0==3 && id1==1))";
    string title = "Process of time sorting";
    string xLabel = "\u0394T";
    int noOfBins = 500;
    int xMin = 0;
    int xMax = 50000;

    string saveFileName = "TimeSorting.pdf";

/* ********************** MAKING FIGURE ************************************************** */

    auto canv = new TCanvas("", "", 1000, 800);
    canv->SetLogy();

/* First figure */
    auto *f= new TFile(filePath.c_str());
    auto *tr=(TTree*)f->Get(treeName.c_str());
    auto hist = new TH1F("hist", title.c_str(), noOfBins, xMin, xMax);
    tr->Draw((branchName +">> hist").c_str(), selectionCrit.c_str());

    hist->SetXTitle(xLabel.c_str());
    hist->SetYTitle("Events");
    canv->SetLeftMargin(0.1);
    hist->GetYaxis()->SetTitleOffset(1);
    hist->GetYaxis()->SetMaxDigits(3); //Setting numbers on axis to x*10^y
    hist->GetXaxis()->SetMaxDigits(3); //Setting numbers on axis to x*10^y
    gStyle ->SetOptStat(kFALSE);
    gStyle ->SetPalette(55);


//Making zoom on the graph in the corner
    TPad *p = new TPad("p", "p", .47, .5, 0.97, 0.95); //Where the overlaying canvas (pad) is placed (numbers between 0 and 1)
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

    canv->cd();
    canv->Update();
    canv->Draw();

    canv->SaveAs(saveFileName.c_str());
}