//
// Created by Christiane Rahbek on 19-05-2023.
//

#include "TH1F.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"

using namespace std;

void distributions180degs(){
    string filePath = "/mnt/d/IS659/finestructure_fit/analysis1/output/doubleDSSD/mergedRuns.root";
    string treeName = "a";
    int BG_del = 10;
    int BG_open = 300;

    string title = "";

    string selectionCrit = "((id0==0 && id1==2) || (id0==1 && id1==3) || (id0==2 && id1==0) || (id0==3 && id1==1)) && abs(FT0-FT1)<1500 && hitAng < 130 && (Edep0<2000||Edep1<2000)";
    string selectionCrit0 = "((id0==0 && id1==2) || (id0==1 && id1==3) || (id0==2 && id1==0) || (id0==3 && id1==1)) && hitAng < 130 && (Edep0<2000||Edep1<2000)";
    string selectionCrit1 = "((id0==0 && id1==2) || (id0==1 && id1==3) || (id0==2 && id1==0) || (id0==3 && id1==1)) && abs(FT0-FT1)<1500 && hitAng < 130";
    string selectionCrit2 = "((id0==0 && id1==2) || (id0==1 && id1==3) || (id0==2 && id1==0) || (id0==3 && id1==1)) && hitAng < 130";
    string selectionCrit_angDistNoTimeCut = "((id0==0 && id1==2) || (id0==1 && id1==3) || (id0==2 && id1==0) || (id0==3 && id1==1))";
    string selectionCrit_angDistTimeCut = "((id0==0 && id1==2) || (id0==1 && id1==3) || (id0==2 && id1==0) || (id0==3 && id1==1)) && abs(FT0-FT1)<1500";
    int xMin = 0;
    int xMax = 3500;
    int noOfBins = 500;

    string saveFileName = "TPROTONS180deg.pdf";
    string xTitle0 = "T_{p} [ms]";
    string xTitle1 = "#DeltaT [ns]";

    /* ************************************************** */

    auto canv = new TCanvas("canv","",1000,300);
    canv->Divide(2,1);

    auto *f= new TFile(filePath.c_str());
    auto *tr0=(TTree*)f->Get(treeName.c_str());
    auto *tr1=(TTree*)f->Get(treeName.c_str());
    auto *tr2=(TTree*)f->Get(treeName.c_str());
    auto *tr3=(TTree*)f->Get(treeName.c_str());
    auto *tr4=(TTree*)f->Get(treeName.c_str());
    auto *tr5=(TTree*)f->Get(treeName.c_str());

    /*************************** Making figures ******************************/
/*
    canv->cd(1);
    auto hist3 = new TH2D("hist3", "Energy distribution, high background", 1000,0,8,1000,0,8);
    tr4->Draw("Edep0/1000:Edep1/1000 >> hist3", selectionCrit1.c_str(), "colz");

    hist3->SetXTitle("E2 [MeV]");
    hist3->SetYTitle("E1 [MeV]");
    hist3->GetYaxis()->SetTitleOffset(1);
    hist3->GetXaxis()->SetTitleOffset(0.9);

    hist3->SetTitleSize(0.06);

    hist3->GetYaxis()->SetTitleSize(0.05);
    hist3->GetXaxis()->SetTitleSize(0.05);
    hist3->GetYaxis()->SetLabelSize(0.05);
    hist3->GetXaxis()->SetLabelSize(0.05);

    gPad->SetTickx();
    gPad->SetTicky();

    hist3->SetStats(kFALSE);

    // ****************Timedist**************************
    canv->cd(2);
    auto hist4 = new TH1F("hist4", "Time between coincidences, high background", noOfBins, 0, 50000);
    tr5->Draw("abs(FT0-FT1) >> hist4", selectionCrit2.c_str());

    canv->GetPad(2)->SetLogy();

    hist4->SetXTitle("#DeltaT [ns]");
    hist4->SetYTitle("Entries/bin");
    hist4->GetXaxis()->SetMaxDigits(3);
    hist4->GetYaxis()->SetMaxDigits(3);
    hist4->GetYaxis()->SetTitleOffset(1);
    hist4->GetXaxis()->SetTitleOffset(0.9);

    hist4->SetTitleSize(0.06);

    hist4->GetYaxis()->SetTitleSize(0.05);
    hist4->GetXaxis()->SetTitleSize(0.05);
    hist4->GetYaxis()->SetLabelSize(0.05);
    hist4->GetXaxis()->SetLabelSize(0.05);

    gPad->SetTickx();
    gPad->SetTicky();

    hist4->SetStats(kFALSE);

*/
    /************************ Energy distibution *****************************/
    canv->cd(1);
    auto hist5 = new TH2D("hist5", "Energy distribution", 1000,0,8,1000,0,8);
    tr3->Draw("Edep0/1000:Edep1/1000 >> hist5", selectionCrit.c_str(), "colz");

    hist5->SetXTitle("E2 [MeV]");
    hist5->SetYTitle("E1 [MeV]");
    hist5->GetYaxis()->SetTitleOffset(0.8);
    hist5->GetXaxis()->SetTitleOffset(0.8);

    hist5->SetTitleSize(0.06);

    hist5->GetYaxis()->SetTitleSize(0.06);
    hist5->GetXaxis()->SetTitleSize(0.06);
    hist5->GetYaxis()->SetLabelSize(0.06);
    hist5->GetXaxis()->SetLabelSize(0.06);
    hist5->GetZaxis()->SetLabelSize(0.06);

    gPad->SetTickx();
    gPad->SetTicky();

    hist5->SetStats(kFALSE);

    /****************timedist***************************/
    canv->cd(2);
    auto hist6 = new TH1F("hist6", "Time between coincidences", noOfBins, 0, 50000);
    tr0->Draw("abs(FT0-FT1) >> hist6", selectionCrit0.c_str());

    canv->GetPad(2)->SetLogy();

    hist6->SetXTitle("#DeltaT [ns]");
    hist6->SetYTitle("Entries/bin [ns^{-1}]");
    hist6->GetXaxis()->SetMaxDigits(3);
    hist6->GetYaxis()->SetMaxDigits(3);
    hist6->GetYaxis()->SetTitleOffset(0.8);
    hist6->GetXaxis()->SetTitleOffset(0.8);
    hist6->SetTitleSize(0.06);

    hist6->GetYaxis()->SetTitleSize(0.06);
    hist6->GetXaxis()->SetTitleSize(0.06);
    hist6->GetYaxis()->SetLabelSize(0.06);
    hist6->GetXaxis()->SetLabelSize(0.06);

    gPad->SetTickx();
    gPad->SetTicky();

    hist6->SetStats(kFALSE);

    /********************** Update and save ****************************/
    canv->Update();
    canv->Draw();
    canv->SaveAs("EnergyTimeDists180Deg.pdf");

    /****************TPROTONS***************************/
    auto canv2 = new TCanvas();
    auto hist0 = new TH1F("hist0", "T_{p} for reduced 180#circ coincidence data", noOfBins, xMin, xMax);
    tr0->Draw("TPROTONS >> hist0", selectionCrit.c_str());

    hist0->GetXaxis()->SetRangeUser(xMin,xMax);

    int fitMin = BG_open+BG_del;
    int fitMax = 800;

    cout << "Fitting TPROTONS... " << endl;
    hist0->Fit("expo", "V", "", fitMin, fitMax);

    hist0->SetXTitle(xTitle0.c_str());
    hist0->SetYTitle("Entries/bin [ms^{-1}]");
    hist0->GetXaxis()->SetMaxDigits(3);
    hist0->GetYaxis()->SetMaxDigits(3);
    hist0->GetYaxis()->SetTitleOffset(1);
    hist0->GetXaxis()->SetTitleOffset(0.9);

    cout << "DOF = " << hist0->GetXaxis()->FindBin(fitMax) - hist0->GetXaxis()->FindBin(fitMin) - 2 << endl;

    hist0->SetTitleSize(0.06);

    hist0->GetYaxis()->SetTitleSize(0.05);
    hist0->GetXaxis()->SetTitleSize(0.05);
    hist0->GetYaxis()->SetLabelSize(0.05);
    hist0->GetXaxis()->SetLabelSize(0.05);

    gPad->SetTickx();
    gPad->SetTicky();

    hist0->SetStats(kFALSE);

    canv2->Update();
    canv2->Draw();
    canv2->SaveAs(saveFileName.c_str());

    /****************** ANGULAR DISTRIBUTION *************************/
    /*
    auto canv1 = new TCanvas("", "", 1000, 750);

    auto hist1 = new TH1F("All coincidences", "", 101, 80, 180);
    tr1->Draw("hitAng >> All coincidences", selectionCrit_angDistNoTimeCut.c_str());
    auto hist2 = new TH1F("Time correlated coincidences", "", 101, 80, 180);
    tr2->Draw("hitAng >> Time correlated coincidences", selectionCrit_angDistTimeCut.c_str(), "SAME");

    hist1->SetLineColor(kBlue);
    hist2->SetLineColor(kRed);

    hist1->SetXTitle("#theta");
    hist1->SetYTitle("Entries/degree");

    canv1->SetLogy();
    canv1->BuildLegend();

    hist1->SetStats(kFALSE);
    hist2->SetStats(kFALSE);

    gPad->SetTickx();
    gPad->SetTicky();

    canv1->Update();
    canv1->Draw();
    //canv1->SaveAs("angularDist180Deg.pdf");
     */

    /* ***************************************** */

}