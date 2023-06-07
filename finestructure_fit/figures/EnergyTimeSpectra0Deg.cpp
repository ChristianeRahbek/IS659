//
// Created by Christiane Rahbek on 12-05-2023.
//

#include "TH2D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TPad.h"
#include <iostream>

using namespace std;

void EnergyTimeSpectra0Deg(){
    /* INSERT INFORMATION HERE*/
    bool withColz = true;
    string filePath = "/mnt/d/IS659/finestructure_fit/analysis1/output/doubleDSSD/mergedRuns.root";
    string treeName = "a";
    string branchNameEnergy = "Edep0/1000:Edep1/1000";
    string branchNameTime = "abs(FT1-FT0)";
    string selectionCritEnergyBasic =
            "id0==id1 && abs(FT0-FT1)<1500";
    string selectionCritTimeBasic =
            "id0==id1";
    string selectionCritEnergyAdvanced =
            "id0==id1 && abs(FT0-FT1)<1500 && (Edep0<1000 || Edep1<1000)";
    string selectionCritTimeAdvanced =
            "id0==id1 && (Edep0<1000 || Edep1<1000)";
    string titleEnergyBasic = "Energy distribution, high background";
    string titleEnergyAdvanced = "Energy distribution, reduced background";
    string titleTimeBasic = "Time between coincidences, high background";
    string titleTimeAdvanced = "Time between coincidences, reduced background";
    string option = "colz";

    string saveFileName = "EnergyTimeSpectra0Deg.pdf";

    string xLabelEnergy = "E1 [MeV]";
    string yLabelEnergy = "E2 [MeV]";
    string xLabelTime = "#DeltaT [ns]";
    string yLabelTime = "Entries/bin [ns^{-1}]";

    int noOfBinsEnergy = 1000;
    int xMinEnergy = 0;
    int xMaxEnergy = 8;
    int noOfBinsTime = 500;
    int xMinTime = 0;
    int xMaxTime = 50000;

    auto canv = new TCanvas("canv","",1000,600);
    canv->Divide(2,2);

    auto *f= new TFile(filePath.c_str());
    auto *tr1=(TTree*)f->Get(treeName.c_str());
    auto *tr2=(TTree*)f->Get(treeName.c_str());
    auto *tr3=(TTree*)f->Get(treeName.c_str());
    auto *tr4=(TTree*)f->Get(treeName.c_str());

    /* **********************FIRST FIGURE**************************** */
    cout << "Entering first subpad" << endl;
    canv->cd(1);

    auto histEnergyBasic = new TH2D("histEnBasic", titleEnergyBasic.c_str(), noOfBinsEnergy, xMinEnergy,
                                    xMaxEnergy, noOfBinsEnergy, xMinEnergy, xMaxEnergy);
    cout << "Drawing first hist" << endl;
    tr1->Draw((branchNameEnergy +">> histEnBasic").c_str(), selectionCritEnergyBasic.c_str(), option.c_str());

    histEnergyBasic->SetXTitle(xLabelEnergy.c_str());
    histEnergyBasic->SetYTitle(yLabelEnergy.c_str());
    histEnergyBasic->GetYaxis()->SetTitleSize(0.06);
    //histEnergyBasic->GetYaxis()->SetTitleOffset(1);
    histEnergyBasic->GetXaxis()->SetTitleSize(0.06);
    histEnergyBasic->GetYaxis()->SetLabelSize(0.06);
    histEnergyBasic->GetXaxis()->SetLabelSize(0.06);
    histEnergyBasic->GetXaxis()->SetTitleOffset(0.85);
    histEnergyBasic->GetZaxis()->SetLabelSize(0.055);
    histEnergyBasic->GetZaxis()->SetLabelOffset(0.001);

    gPad->SetTicky(); gPad->SetTickx();
    
    histEnergyBasic->SetStats(kFALSE);
    canv->GetPad(1)->SetLogz();

    /* **********************SECOND FIGURE**************************** */
    cout << "Entering second subpad" << endl;
    canv->cd(2);

    auto histTimeBasic = new TH1F("histTBasic", titleTimeBasic.c_str(),
                                  noOfBinsTime, xMinTime, xMaxTime);
    cout << "Drawing second hist" << endl;
    tr2->Draw((branchNameTime +">> histTBasic").c_str(), selectionCritTimeBasic.c_str());

    histTimeBasic->SetXTitle(xLabelTime.c_str());
    histTimeBasic->SetYTitle(yLabelTime.c_str());
    histTimeBasic->GetYaxis()->SetTitleSize(0.06);
    histTimeBasic->GetXaxis()->SetTitleSize(0.06);
    histTimeBasic->GetYaxis()->SetLabelSize(0.06);
    histTimeBasic->GetXaxis()->SetLabelSize(0.06);
    histTimeBasic->GetXaxis()->SetMaxDigits(3);
    histTimeBasic->GetYaxis()->SetTitleOffset(0.75);
    histTimeBasic->GetXaxis()->SetTitleOffset(0.5);

    gPad->SetTicky(); gPad->SetTickx();
    
    histTimeBasic->SetStats(kFALSE);
    canv->GetPad(2)->SetLogy();

    /* **********************THIRD FIGURE**************************** */
    cout << "Entering third subpad" << endl;
    canv->cd(3);

    auto histEnergyAdvanced = new TH2D("histEnAdv", titleEnergyAdvanced.c_str(), noOfBinsEnergy, xMinEnergy,
                                    xMaxEnergy, noOfBinsEnergy, xMinEnergy, xMaxEnergy);
    cout << "Drawing third hist" << endl;
    tr3->Draw((branchNameEnergy +">> histEnAdv").c_str(), selectionCritEnergyAdvanced.c_str(), option.c_str());

    histEnergyAdvanced->SetXTitle(xLabelEnergy.c_str());
    histEnergyAdvanced->SetYTitle(yLabelEnergy.c_str());
    histEnergyAdvanced->GetYaxis()->SetTitleSize(0.06);
    histEnergyAdvanced->GetXaxis()->SetTitleSize(0.06);
    histEnergyAdvanced->GetYaxis()->SetLabelSize(0.06);
    histEnergyAdvanced->GetXaxis()->SetLabelSize(0.06);
    histEnergyAdvanced->GetXaxis()->SetTitleOffset(0.85);
    histEnergyAdvanced->GetZaxis()->SetLabelSize(0.055);
    histEnergyAdvanced->GetZaxis()->SetLabelOffset(0.001);

    gPad->SetTicky(); gPad->SetTickx();
    
    histEnergyAdvanced->SetStats(kFALSE);
    canv->GetPad(3)->SetLogz();

    /* **********************FOURTH FIGURE**************************** */
    cout << "Entering fourth subpad" << endl;
    canv->cd(4);

    auto histTimeAdvanced = new TH1F("histTAdv", titleTimeAdvanced.c_str(),
                                  noOfBinsTime, xMinTime, xMaxTime);
    cout << "Drawing fourth hist" << endl;
    tr4->Draw((branchNameTime +">> histTAdv").c_str(), selectionCritTimeAdvanced.c_str());

    histTimeAdvanced->SetXTitle(xLabelTime.c_str());
    histTimeAdvanced->SetYTitle(yLabelTime.c_str());
    histTimeAdvanced->GetYaxis()->SetTitleSize(0.06);
    histTimeAdvanced->GetXaxis()->SetTitleSize(0.06);
    histTimeAdvanced->GetYaxis()->SetLabelSize(0.06);
    histTimeAdvanced->GetXaxis()->SetLabelSize(0.06);
    histTimeAdvanced->GetXaxis()->SetMaxDigits(3);
    histTimeAdvanced->GetYaxis()->SetTitleOffset(0.75);
    histTimeAdvanced->GetXaxis()->SetTitleOffset(0.5);

    gPad->SetTicky(); gPad->SetTickx();
    
    histTimeAdvanced->SetStats(kFALSE);
    canv->GetPad(4)->SetLogy();

    /* **********************DRAWING AND SAVING FIGURE**************************** */

    canv->Update();
    canv->Draw();

    canv->SaveAs(saveFileName.c_str());
}