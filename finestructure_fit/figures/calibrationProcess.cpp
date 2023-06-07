//
// Made by Christiane on 18/01/23
//

#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TStyle.h"

using namespace std;
//using namespace AUSA;

void calibrationProcess() {
    /* INSERT INFORMATION HERE*/
    string filePath = "/mnt/d/IS659/runs/Run212.root";
    string treeName = "h101";
    string branchName = "U1F_E";
    string selectionCrit = "U1FI == 12";
    string title = "Output from DAQ";
    string xLabel = "Channel number";
    int noOfBins = 3500;
    int xMin = 2000;
    int xMax = 5500;

    string filePath1 = "/mnt/d/IS659/finestructure_fit/calibrationAnalysis/output/Run212mlio.root";
    string treeName1 = "a";
    string branchName1 = "FE/1000";
    string selectionCrit1 = "FI == 12";
    string title1 = "Calibrated Spectrum";
    string xLabel1 = "Energy [MeV]";
    int noOfBins1 = 5000;
    int xMin1 = 2;
    int xMax1 = 7;

    string saveFileName = "CalibrationProcess.pdf";

    /* ********************** MAKING FIGURE ************************************************** */

    auto canv = new TCanvas("", "", 1200, 500);

    canv->Divide(2,1);
    canv->cd(1);

    /* First figure */
    auto *f= new TFile(filePath.c_str());
    auto *tr=(TTree*)f->Get(treeName.c_str());

    auto hist = new TH1F("hist", title.c_str(), noOfBins, xMin, xMax);
    tr->Draw((branchName +">> hist").c_str(), selectionCrit.c_str());

    hist ->SetXTitle(xLabel.c_str());
    hist->SetYTitle("Events/bin [Channel^{-1}]");
    gStyle ->SetOptStat(kFALSE);
    gStyle ->SetPalette(55);
    gPad->SetTickx();
    gPad->SetTicky();
    hist->GetYaxis()->SetMaxDigits(2);
    hist->GetXaxis()->SetMaxDigits(2);
    hist->GetYaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetYaxis()->SetTitleOffset(0.95);
    hist->GetXaxis()->SetTitleOffset(0.8);

    canv->cd(2);

    /* Second figure */

    auto *f1= new TFile(filePath1.c_str());
    auto *tr1=(TTree*)f1->Get(treeName1.c_str());

    auto hist1 = new TH1F("hist1", title1.c_str(), noOfBins1, xMin1, xMax1);
    tr1->Draw((branchName1 +" >> hist1").c_str(), selectionCrit1.c_str());

    hist1 ->SetXTitle(xLabel1.c_str());
    hist1->SetYTitle("Events/bin [MeV^{-1}]");
    gStyle ->SetOptStat(kFALSE);
    gStyle ->SetPalette(55);
    gPad->SetTickx();
    gPad->SetTicky();
    hist1->GetYaxis()->SetMaxDigits(2);
    hist1->GetXaxis()->SetMaxDigits(2);
    hist1->GetYaxis()->SetLabelSize(0.05);
    hist1->GetYaxis()->SetTitleSize(0.05);
    hist1->GetXaxis()->SetLabelSize(0.05);
    hist1->GetXaxis()->SetTitleSize(0.05);
    hist1->GetYaxis()->SetTitleOffset(0.95);
    hist1->GetXaxis()->SetTitleOffset(0.8);


    canv->SaveAs(saveFileName.c_str());
}
