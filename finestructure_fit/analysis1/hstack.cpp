//
//
//

#include "THStack.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"

void hstack(bool bat) {
    //defining each histogram
    //bool bat = false;
    std::string addition;
    if(bat) addition = "bat/";
    else addition = "";

    auto *f= new TFile(("/mnt/d/IS659/finestructure_fit/analysis1/output/" + addition + "Run167mlio.root").c_str());
    auto *tr=(TTree*)f->Get("a");


    auto h1st = new TH1F("h1st", "h1st", 16, 1, 17);
    tr -> Draw("FI >> h1st");
    h1st ->SetLineColor(kRed);

    auto h2st = new TH1F("h2st", "h2st", 16, 1, 17);
    tr -> Draw("FI >> h2st", "Edep1[0] > 0 && Edep3[1] > 0 && abs(Edep1[0]-Edep3[1])<400 && abs(BT[0]-BT[1]) < 800");
    h2st ->SetLineColor(kBlue);

    auto h3st = new TH1F("h3st", "h3st", 16, 1, 17);
    tr -> Draw("FI >> h3st", "Edep0[0] > 0 && Edep2[1] > 0 && abs(Edep0[0]-Edep2[1])<400 && abs(BT[0]-BT[1]) < 800");
    h3st ->SetLineColor(kGreen);

    auto h4st = new TH1F("h4st", "h4st", 16, 1, 17);
    tr -> Draw("FI >> h4st", "abs(BT[0]-BT[1]) < 800");
    h4st ->SetLineColor(kOrange);

     /*
    auto h4st = new TH1F("h4st", "h4st", 1000, 200 , 8000);
    tr -> Draw("Edep2 >> h4st");
    h4st ->SetLineColor(kMagenta);

    auto h5st = new TH1F("h5st", "h5st", 1000, 200 , 8000);
    tr -> Draw("Edep3 >> h5st");
    h5st ->SetLineColor(kOrange);
    */

    //Plotting every histogram together
    auto hs = new THStack("hs","Strip distribution");
    hs->Add(h1st, "hist");
    hs->Add(h2st, "hist");
    hs->Add(h3st, "hist");
    hs->Add(h4st, "hist");
    //hs->Add(h5st, "hist");

    auto l = new TLegend();

    l -> AddEntry(h1st, "Without gating", "L");
    l -> AddEntry(h2st, "Det 2 & 4, main peak without shoulders", "L");
    l -> AddEntry(h3st, "Det 1 & 3, main peak without shoulders", "L");
    l -> AddEntry(h4st, "Main peak w.o. shoulders, all detectors", "L");
    //l -> AddEntry(h5st, "Edep3", "L");

    hs->Draw("nostack");
    l->Draw();

    hs->SaveAs(("/mnt/d/IS659/finestructure_fit/analysis1/output/" + addition + "Strip_distribution_Back.png").c_str());
/*
 * UNCOMMENT TO STACK HISTOGRAMS AS WELL
    TCanvas *cst = new TCanvas("cst","stacked hists",10,10,700,700);

    //Stacking only Edep
    auto hs1 = new THStack("hs1","Edep");
    auto hs11 = new THStack("hs11", "Compare Edep with Edssd");
    hs1->Add(h2st, "hist");
    hs1->Add(h3st, "hist");
    hs1->Add(h4st, "hist");
    hs1->Add(h5st, "hist");

    hs1->Draw();

    //hs11->Add(hs1->GetHistogram());
    //hs11->Draw("nostak");

    auto l1 = new TLegend();

    l1 -> AddEntry(h2st, "Edep0", "L");
    l1 -> AddEntry(h3st, "Edep0, Edep1", "L");
    l1 -> AddEntry(h4st, "Edep0, Edep1, Edep2", "L");
    l1 -> AddEntry(h5st, "Edep0, Edep1, Edep2, Edep3", "L");

    l1->Draw();

    hs1->SaveAs(("/mnt/d/IS659/finestructure_fit/analysis1/output/" + addition + "CumEdep.png").c_str());
*/
}

