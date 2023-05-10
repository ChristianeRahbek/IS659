//
// Created by Christiane Rahbek on 09-05-2023.
//

#include <iostream>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>

using namespace std;

void dssdMultiplicityPlot(){
    string filePath = "/mnt/d/IS659/finestructure_fit/analysis1/output/preAnalysis/mergedRuns.root";
    string treeName = "a";
    string branchName = "dssdMul";

    string xLabel = "DSSD multiplicity";
    string yLabel = "N";
    string option = "";

    int noOfBins = 9;
    int xMin = 0;
    int xMax = 9;

    auto canv = new TCanvas("", "", 1200, 550);

    canv->cd(1);

    auto *f = new TFile(filePath.c_str());
    auto *tr = (TTree*)f->Get(treeName.c_str());
    auto hist = new TH1D("hist", "", noOfBins, xMin, xMax);

    hist->GetXaxis()->SetRangeUser(0,8);
    hist->GetXaxis()->SetNdivisions(9, false);
    hist->GetXaxis()->SetBinLabel(1, "0");
    hist->GetXaxis()->SetBinLabel(2, "1");
    hist->GetXaxis()->SetBinLabel(3, "2");
    hist->GetXaxis()->SetBinLabel(4, "3");
    hist->GetXaxis()->SetBinLabel(5, "4");
    hist->GetXaxis()->SetBinLabel(6, "5");
    hist->GetXaxis()->SetBinLabel(7, "6");
    hist->GetXaxis()->SetBinLabel(8, "7");
    hist->GetXaxis()->SetBinLabel(9, "8");
    hist->LabelsOption("h");

    hist->GetXaxis()->SetLabelSize(0.06); //size is in percent of pad size
    hist->GetYaxis()->SetLabelSize();
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetYaxis()->SetTitleSize(0.05);

    tr->Draw((branchName +">> hist").c_str(), "", option.c_str());

    hist->SetXTitle(xLabel.c_str());
    hist->SetYTitle(yLabel.c_str());

    hist->SetStats(kFALSE);

    gPad->SetTickx();
    gPad->SetTicky();

    canv->Update();

    canv->SaveAs("dssdMultiplicityPlot.pdf");
}