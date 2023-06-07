//
// Created by Christiane Rahbek on 12-05-2023.
//

#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"

using namespace std;

void TPROTONS0deg() {
    string filePath = "/mnt/d/IS659/finestructure_fit/analysis1/output/doubleDSSD/mergedRuns.root";
    string treeName = "a";
    int BG_del = 10;
    int BG_open = 300;

    string title = "";
    string xLabel = "  T_p [ms]";

    string selectionCritId0 = "id0==id1 && id0==0 && abs(FT0-FT1)<1500 && (Edep0<1000 || Edep1<1000)";
    string selectionCritId1 = "id0==id1 && id0==1 && abs(FT0-FT1)<1500 && (Edep0<1000 || Edep1<1000)";
    string selectionCritId2 = "id0==id1 && id0==2 && abs(FT0-FT1)<1500 && (Edep0<1000 || Edep1<1000)";
    string selectionCritId3 = "id0==id1 && id0==3 && abs(FT0-FT1)<1500 && (Edep0<1000 || Edep1<1000)";
    string selectionCrit = "id0==id1 && abs(FT0-FT1)<1500 && (Edep0<1000 || Edep1<1000)";

    int xMin = 0;
    int xMax = 3500;
    int noOfBins = 500;

    string saveFileName = "TPROTONS0deg.pdf";
    string xTitle = "T_{p} [ms]";

    /* ********************** Making fit **************************** */

    auto canv = new TCanvas("", "", 1000, 650);
    canv->Divide(2,2);
    //canv->SetLogy();

    /* Figure */
    auto *f= new TFile(filePath.c_str());
    auto *tr0=(TTree*)f->Get(treeName.c_str());
    auto *tr1=(TTree*)f->Get(treeName.c_str());
    auto *tr2=(TTree*)f->Get(treeName.c_str());
    auto *tr3=(TTree*)f->Get(treeName.c_str());
    auto *tr4=(TTree*)f->Get(treeName.c_str());

    /* ***************************************** */
    canv->cd(1);
    auto hist0 = new TH1F("hist0", "DSSD 1", noOfBins, xMin, xMax);
    tr0->Draw("TPROTONS >> hist0", selectionCritId0.c_str());

    hist0->GetXaxis()->SetRangeUser(xMin,xMax);

    //cout << "Fitting for id0... " << endl;
    //hist0->Fit("expo", "V", "", BG_open+2*BG_del, 1000);

    hist0->SetXTitle(xTitle.c_str());
    hist0->SetYTitle("Entries/bin [ms^{-1}]");
    hist0->GetXaxis()->SetMaxDigits(3);
    hist0->GetYaxis()->SetMaxDigits(3);
    hist0->GetXaxis()->SetTitleOffset(0.77);
    hist0->GetYaxis()->SetTitleOffset(0.75);
    hist0->GetXaxis()->SetLabelOffset(0.015);


    gPad->SetTickx();
    gPad->SetTicky();

    hist0->SetStats(kFALSE);

    /* ***************************************** */
    canv->cd(2);
    auto hist1 = new TH1F("hist1", "DSSD 2", noOfBins, xMin, xMax);
    tr1->Draw("TPROTONS >> hist1", selectionCritId1.c_str());

    hist1->GetXaxis()->SetRangeUser(xMin,xMax);

    //cout << "Fitting for id1... " << endl;
    //hist1->Fit("expo", "V", "", BG_open+2*BG_del, 1000);

    hist1->SetXTitle(xTitle.c_str());
    hist1->SetYTitle("Entries/bin [ms^{-1}]");
    hist1->GetXaxis()->SetMaxDigits(3);
    hist1->GetYaxis()->SetMaxDigits(3);
    hist1->GetYaxis()->SetTitleOffset(0.75);
    hist1->GetXaxis()->SetTitleOffset(0.77);
    hist1->GetXaxis()->SetLabelOffset(0.015);

    gPad->SetTickx();
    gPad->SetTicky();

    hist1->SetStats(kFALSE);

    /* ***************************************** */
    canv->cd(3);
    auto hist2 = new TH1F("hist2", "DSSD 3", noOfBins, xMin, xMax);
    tr2->Draw("TPROTONS >> hist2", selectionCritId2.c_str());

    hist2->GetXaxis()->SetRangeUser(xMin,xMax);

    cout << "Fitting for id2... " << endl;
    hist2->Fit("expo", "V", "", BG_open+2*BG_del, 800);

    cout << "DOF = " << hist2->GetXaxis()->FindBin(800) - hist2->GetXaxis()->FindBin(BG_open+2*BG_del) - 2 << endl;

    hist2->SetXTitle(xTitle.c_str());
    hist2->SetYTitle("Entries/bin [ms^{-1}]");
    hist2->GetXaxis()->SetMaxDigits(3);
    hist2->GetYaxis()->SetMaxDigits(3);
    hist2->GetYaxis()->SetTitleOffset(0.75);
    hist2->GetXaxis()->SetTitleOffset(0.77);
    hist2->GetXaxis()->SetLabelOffset(0.015);

    gPad->SetTickx();
    gPad->SetTicky();

    hist2->SetStats(kFALSE);

    /* ***************************************** */
    canv->cd(4);
    auto hist3 = new TH1F("hist3", "DSSD 4", noOfBins, xMin, xMax);
    tr3->Draw("TPROTONS >> hist3", selectionCritId3.c_str());

    hist3->GetXaxis()->SetRangeUser(xMin,xMax);

    //cout << "Fitting for id3... " << endl;
    //hist3->Fit("expo", "V", "", BG_open+2*BG_del, 1000);

    hist3->SetXTitle(xTitle.c_str());
    hist3->SetYTitle("Entries/bin [ms^{-1}]");
    hist3->GetXaxis()->SetMaxDigits(3);
    hist3->GetYaxis()->SetMaxDigits(3);
    hist3->GetYaxis()->SetTitleOffset(0.75);
    hist3->GetXaxis()->SetTitleOffset(0.77);
    hist3->GetXaxis()->SetLabelOffset(0.015);

    gPad->SetTickx();
    gPad->SetTicky();

    hist3->SetStats(kFALSE);

    /* ***************************************** */
    hist0->SetTitleSize(15);
    hist1->SetTitleSize(0.07);
    hist2->SetTitleSize(0.07);
    hist3->SetTitleSize(0.07);

    hist0->GetYaxis()->SetTitleSize(0.06);
    hist0->GetXaxis()->SetTitleSize(0.06);
    hist0->GetYaxis()->SetLabelSize(0.06);
    hist0->GetXaxis()->SetLabelSize(0.06);

    hist1->GetYaxis()->SetTitleSize(0.06);
    hist1->GetXaxis()->SetTitleSize(0.06);
    hist1->GetYaxis()->SetLabelSize(0.06);
    hist1->GetXaxis()->SetLabelSize(0.06);

    hist2->GetYaxis()->SetTitleSize(0.06);
    hist2->GetXaxis()->SetTitleSize(0.06);
    hist2->GetYaxis()->SetLabelSize(0.06);
    hist2->GetXaxis()->SetLabelSize(0.06);

    hist3->GetYaxis()->SetTitleSize(0.06);
    hist3->GetXaxis()->SetTitleSize(0.06);
    hist3->GetYaxis()->SetLabelSize(0.06);
    hist3->GetXaxis()->SetLabelSize(0.06);


    canv->Update();
    canv->Draw();
    canv->SaveAs(saveFileName.c_str());

    /* *****************NEW CANVAS******************* */
    auto canv1 = new TCanvas("", "", 1000, 650);

    auto hist4 = new TH1F("hist4", "", noOfBins, xMin, xMax);
    tr4->Draw("TPROTONS >> hist4", selectionCrit.c_str());

    hist4->GetXaxis()->SetRangeUser(xMin,xMax);

    cout << "Fitting for the new canvas... " << endl;
    hist4->Fit("expo", "V", "", BG_open+2*BG_del, 800);

    hist4->SetXTitle(xTitle.c_str());
    hist4->SetYTitle("Entries/bin [ms^{-1}]");
    hist4->GetXaxis()->SetMaxDigits(3);
    hist4->GetYaxis()->SetMaxDigits(3);
    hist4->GetXaxis()->SetTitleOffset(0.77);
    hist4->GetYaxis()->SetTitleOffset(0.75);

    hist4->GetYaxis()->SetTitleSize(0.06);
    hist4->GetXaxis()->SetTitleSize(0.06);
    hist4->GetYaxis()->SetLabelSize(0.06);
    hist4->GetXaxis()->SetLabelSize(0.06);

    gPad->SetTickx();
    gPad->SetTicky();

    hist4->SetStats(kFALSE);

    canv1->Update();
    canv1->Draw();
    canv1->SaveAs("TPROTONS0degAllDSSDs.pdf");
}