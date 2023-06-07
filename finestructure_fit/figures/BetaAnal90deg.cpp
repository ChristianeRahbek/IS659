//
// Created by Christiane Rahbek on 02-06-2023.
//

#include "TH2D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TPad.h"
#include "Math/GoFTest.h"
#include <iostream>

using namespace std;

void BetaAnal90deg() {
    /* INSERT INFORMATION HERE*/

    string filePath = "/mnt/d/IS659/finestructure_fit/analysis1/output/doubleDSSD/mergedRuns.root";
    string treeName = "a";
    string branchName = "Edep0/1000:Edep1/1000";
    string selectionCritNoBeta =
            "((id0 == 0 && (id1 == 3 || id1 == 1)) || (id0 == 1 && (id1 == 0 || id1 == 2)) || (id0 == 2 && (id1 == 1 || id1 == 3)) || (id0 == 3 && (id1 == 2 || id1 == 0)))  && abs(FT0-FT1)<2500 && isParticleBeta==0 && isBetas==0";
    string selectionCrit1Beta =
            "((id0 == 0 && (id1 == 3 || id1 == 1)) || (id0 == 1 && (id1 == 0 || id1 == 2)) || (id0 == 2 && (id1 == 1 || id1 == 3)) || (id0 == 3 && (id1 == 2 || id1 == 0)))  && abs(FT0-FT1)<2500 && (isParticleBeta==1 || isBetas==1)";
    string titleNoBeta = "No electron coincidences";
    string title1Beta = "Possible electron coincidence(s)";
    string title2Betas = "Maximum 2 electron coincidences";
    string xLabel = "E1 [MeV]";
    string yLabel = "E2 [MeV]";
    string option = "colz";

    int noOfBins = 1000;
    int xMin = 0;
    int xMax = 8;

    string saveFileDir = "EnergyDistributions/";

    auto canv = new TCanvas("", "", 1200, 550);
    canv->Divide(2,1);

    canv->cd(1);

    auto *fNoBeta= new TFile(filePath.c_str());
    auto *trNoBeta=(TTree*)fNoBeta->Get(treeName.c_str());
    auto histNoBeta = new TH2D("noBeta", titleNoBeta.c_str(), noOfBins, xMin, xMax, noOfBins, xMin, xMax);
    trNoBeta->Draw((branchName +">> noBeta").c_str(), selectionCritNoBeta.c_str(), option.c_str());

    histNoBeta->SetXTitle(xLabel.c_str());
    histNoBeta->SetYTitle(yLabel.c_str());
    histNoBeta->GetXaxis()->SetTitleSize(0.05);
    histNoBeta->GetXaxis()->SetLabelSize(0.05);
    histNoBeta->GetYaxis()->SetTitleSize(0.05);
    histNoBeta->GetYaxis()->SetLabelSize(0.05);

    histNoBeta->SetStats(kFALSE);

    gPad->SetTickx();
    gPad->SetTicky();


/* ********************** MAKING FIGURE FOR 1 BETA ************************************************** */

    canv->cd(2);

    auto *f1Beta= new TFile(filePath.c_str());
    auto *tr1Beta=(TTree*)f1Beta->Get(treeName.c_str());
    auto hist1Beta = new TH2D("1Beta", title1Beta.c_str(), noOfBins, xMin, xMax, noOfBins, xMin, xMax);
    tr1Beta->Draw((branchName +">> 1Beta").c_str(), selectionCrit1Beta.c_str(), option.c_str());

    hist1Beta->SetXTitle(xLabel.c_str());
    hist1Beta->SetYTitle(yLabel.c_str());

    hist1Beta->GetXaxis()->SetTitleSize(0.05);
    hist1Beta->GetXaxis()->SetLabelSize(0.05);
    hist1Beta->GetYaxis()->SetTitleSize(0.05);
    hist1Beta->GetYaxis()->SetLabelSize(0.05);

    hist1Beta->SetStats(kFALSE);

    gPad->SetTickx();
    gPad->SetTicky();

/* ********************** SAVING FIGURE ************************************************** */
    canv->Update();
    canv->SaveAs((saveFileDir + "90degBetaAnalysis.pdf").c_str());
}