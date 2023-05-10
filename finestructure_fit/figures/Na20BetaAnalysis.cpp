//
// Created by Christiane Rahbek on 09-05-2023.
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

void Na20BetaAnalysis() {
    /* INSERT INFORMATION HERE*/

    string filePath = "/mnt/d/IS659/finestructure_fit/analysis1/output/doubleDSSD/Run208mliolio.root";
    string treeName = "a";
    string branchName = "Edep0/1000:Edep1/1000";
    string selectionCritNoBeta =
            "abs(FT0-FT1)<1500 && isParticleBeta==0 && isBetas==0";
    string selectionCrit1Beta =
            "abs(FT0-FT1)<1500 && isParticleBeta==1 && isBetas==0";
    string selectionCrit2Betas =
            "abs(FT0-FT1)<1500 && isParticleBeta==0 && isBetas==1";
    string titleNoBeta = "No electron coincidences";
    string title1Beta = "Maximum 1 electron coincidence";
    string title2Betas = "Maximum 2 electron coincidences";
    string xLabel = "E1 [MeV]";
    string yLabel = "E2 [MeV]";
    string option = "";

    int noOfBins = 1000;
    int xMin = 0;
    int xMax = 6;

    string saveFileDir = "EnergyDistributions/";

    auto canv = new TCanvas("", "", 1200, 550);
    canv->Divide(2,1);

/* ********************** MAKING FIGURE FOR 2 BETAS ************************************************** */
/*
    canv->cd(3);

    auto *f2Betas= new TFile(filePath.c_str());
    auto *tr2Betas=(TTree*)f2Betas->Get(treeName.c_str());
    auto hist2Betas = new TH2D("2Betas", title2Betas.c_str(), noOfBins, xMin, xMax, noOfBins, xMin, xMax);
    tr2Betas->Draw((branchName +">> 2Betas").c_str(), selectionCrit2Betas.c_str(), option.c_str());

    hist2Betas->SetXTitle(xLabel.c_str());
    hist2Betas->SetYTitle(yLabel.c_str());

    hist2Betas->SetStats(kFALSE);
*/

/* ********************** MAKING FIGURE FOR NO BETAS ************************************************** */

    canv->cd(1);

    auto *fNoBeta= new TFile(filePath.c_str());
    auto *trNoBeta=(TTree*)fNoBeta->Get(treeName.c_str());
    auto histNoBeta = new TH2D("noBeta", titleNoBeta.c_str(), noOfBins, xMin, xMax, noOfBins, xMin, xMax);
    trNoBeta->Draw((branchName +">> noBeta").c_str(), selectionCritNoBeta.c_str(), option.c_str());

    histNoBeta->SetXTitle(xLabel.c_str());
    histNoBeta->SetYTitle(yLabel.c_str());

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

    hist1Beta->SetStats(kFALSE);

    gPad->SetTickx();
    gPad->SetTicky();

/* ********************** SAVING FIGURE ************************************************** */
    canv->Update();
    canv->SaveAs((saveFileDir + "20NaBetaAnalysis.pdf").c_str());

/* ********************** KS-TEST ************************************************** */

    std::vector<double> noBeta, oneBeta;
    auto sampleNoBeta = new Double_t[histNoBeta->GetEntries()];
    auto sample1Beta = new Double_t[hist1Beta->GetEntries()];
    //auto sample2Betas = new Double_t[hist2Betas->GetEntries()];

    cout << histNoBeta->GetEntries() << endl;
    cout << hist1Beta->GetEntries() << endl;

    cout << "for loop start" << endl;
    for(int xbin = 1; xbin <= histNoBeta->GetNbinsX(); xbin++) {
        for(int ybin = 1; ybin <= histNoBeta->GetNbinsY(); ybin++) {
            double xy = histNoBeta->GetXaxis()->GetBinCenter(xbin) + histNoBeta->GetYaxis()->GetBinCenter(ybin);
            for(double z = 1; z <= histNoBeta->GetBinContent(xbin,ybin); z++) {
                //cout << "xy = " << xy << endl;
                noBeta.emplace_back(xy);
            }
        }
    }
    cout << "for loop start pt2" << endl;
    for(int xbin = 1; xbin <= hist1Beta->GetNbinsX(); xbin++) {
        for(int ybin = 1; ybin <= hist1Beta->GetNbinsY(); ybin++) {
            double xy = hist1Beta->GetXaxis()->GetBinCenter(xbin) + hist1Beta->GetYaxis()->GetBinCenter(ybin);
            for(double z = 1; z <= hist1Beta->GetBinContent(xbin,ybin); z++) {
                oneBeta.emplace_back(xy);
            }
        }
    }

    for(UInt_t i = 0; i < noBeta.size(); i++){
        //cout << noBeta[i] << endl;
        sampleNoBeta[i] = noBeta[i];
    }
    for(UInt_t i = 0; i < oneBeta.size(); i++){
        sampleNoBeta[i] = oneBeta[i];
    }
/*
    for(UInt_t i = 1; i<= hist2Betas->GetEntries(); i++){
        sample2Betas[i-1] = hist2Betas->GetBinContent(i);
    }
*/

    cout << "did I make it here" << endl;

    auto* testNoBeta1Beta = new ROOT::Math::GoFTest(noBeta.size(), sampleNoBeta,
                                                    oneBeta.size(), sample1Beta);
    //auto* testNoBeta2Betas = new ROOT::Math::GoFTest(histNoBeta->GetEntries(), sampleNoBeta,
    //                                                hist2Betas->GetEntries(), sample2Betas);

    cout << "??" << endl;
    auto pvalNoBeta1Beta = testNoBeta1Beta->KolmogorovSmirnov2SamplesTest();
    //auto pvalNoBeta2Betas = testNoBeta2Betas->KolmogorovSmirnov2SamplesTest();

    cout << "Probability that No beta and 1 beta is from same distribution:" << endl;
    cout << pvalNoBeta1Beta << endl;
    //cout << "Probability that No beta and 2 betas is from same distribution:" << endl;
    //cout << pvalNoBeta2Betas << endl;

    auto canv1 = new TCanvas("", "", 1500, 500);
    canv1->Divide(2,1);
    canv1->cd(1);

    TH1D* hist = new TH1D();

    for(int i = 0; i < noBeta.size(); i++){
        hist->Fill(noBeta[i]);
    }

    hist->Draw();

    canv1->cd(2);

    TH1D* hist1 = new TH1D();

    for(int i = 0; i < oneBeta.size(); i++){
        hist1->Fill(oneBeta[i]);
    }

    hist1->Draw();

}