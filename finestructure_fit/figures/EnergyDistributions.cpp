//
// Created by Christiane Rahbek on 28-03-2023.
//

#include "TH2D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TPad.h"

using namespace std;

void EnergyDistributions() {
    /* INSERT INFORMATION HERE*/
    bool withColz = true;
    string filePath = "/mnt/d/IS659/finestructure_fit/analysis1/output/doubleDSSD/mergedRuns.root";
    string treeName = "a";
    string branchName = "Edep0/1000:Edep1/1000";
    string selectionCrit0 =
            "id0==id1 && abs(FT0-FT1)<1500 && isParticleBeta==0 && isBetas==0";
    string selectionCrit90 =
            "((id0 == 0 && (id1 == 3 || id1 == 1)) || (id0 == 1 && (id1 == 0 || id1 == 2)) || (id0 == 2 && (id1 == 1 || id1 == 3)) || (id0 == 3 && (id1 == 2 || id1 == 0))) && abs(FT0-FT1)<3000 && isParticleBeta==0 && isBetas==0 && abs(Edep0-Edep1)>500";
    string selectionCrit180 =
            "((id0==0 && id1==2) || (id0==1 && id1==3) || (id0==2 && id1==0) || (id0==3 && id1==1)) && abs(FT0-FT1)<1500 && isParticleBeta==0 && isBetas==0 && hitAng < 130";
    string title0 = "0#circ Coincidences";
    string title90 = "90#circ Coincidences";
    string title180 = "180#circ Coincidences";
    string xLabel = "E1 [MeV]";
    string yLabel = "E2 [MeV]";
    string option;

    if (withColz) option = "colz";
    else option = "";

    int noOfBins = 1000;
    int xMin = 0;
    int xMax = 8;

    string saveFileDir = "EnergyDistributions/";

    auto canv0 = new TCanvas("", "", 1000, 800);
    auto canv90 = new TCanvas("", "", 1000, 800);
    auto canv180 = new TCanvas("", "", 1000, 800);

/* ********************** MAKING FIGURE FOR 0 DEG ************************************************** */

    canv0->cd();

    auto *f0= new TFile(filePath.c_str());
    auto *tr0=(TTree*)f0->Get(treeName.c_str());
    auto hist0 = new TH2D("hist", title0.c_str(), noOfBins, xMin, xMax, noOfBins, xMin, xMax);
    tr0->Draw((branchName +">> hist").c_str(), selectionCrit0.c_str(), option.c_str());

    hist0->SetXTitle(xLabel.c_str());
    hist0->SetYTitle(yLabel.c_str());

    hist0->SetStats(kFALSE);

    canv0->Draw();
    if(withColz) {
        canv0->SaveAs((saveFileDir + "0Deg_colz.pdf").c_str());
    } else{
        canv0->SaveAs((saveFileDir + "0Deg.pdf").c_str());
    }

    /* ********************** MAKING FIGURE FOR 90 DEG ************************************************** */

    canv90->cd();

    auto *f90= new TFile(filePath.c_str());
    auto *tr90=(TTree*)f90->Get(treeName.c_str());
    auto hist90 = new TH2D("hist", title90.c_str(), noOfBins, xMin, xMax, noOfBins, xMin, xMax);
    tr90->Draw((branchName +">> hist").c_str(), selectionCrit90.c_str(), option.c_str());

    hist90->SetXTitle(xLabel.c_str());
    hist90->SetYTitle(yLabel.c_str());

    hist90->SetStats(kFALSE);

    canv90->Draw();
    if(withColz) {
        canv90->SaveAs((saveFileDir + "90Deg_colz.pdf").c_str());
    } else{
        canv90->SaveAs((saveFileDir + "90Deg.pdf").c_str());
    }

    /* ********************** MAKING FIGURE FOR 180 DEG ************************************************** */

    canv180->cd();

    auto *f180= new TFile(filePath.c_str());
    auto *tr180=(TTree*)f180->Get(treeName.c_str());
    auto hist180 = new TH2D("hist", title180.c_str(), noOfBins, xMin, xMax, noOfBins, xMin, xMax);
    tr180->Draw((branchName +">> hist").c_str(), selectionCrit180.c_str(), option.c_str());

    hist180->SetXTitle(xLabel.c_str());
    hist180->SetYTitle(yLabel.c_str());

    hist180->SetStats(kFALSE);

    canv180->Draw();
    if(withColz) {
        canv180->SaveAs((saveFileDir + "180Deg_colz.pdf").c_str());
    } else{
        canv180->SaveAs((saveFileDir + "180Deg.pdf").c_str());
    }
}