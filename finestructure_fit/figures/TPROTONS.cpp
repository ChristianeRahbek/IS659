//
// Made by Christiane 23/01-2023
//

#include "THStack.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"

void TPROTONS() {
    //defining each histogram

    auto *f= new TFile("/mnt/d/IS659/finestructure_fit/analysis1/output/Run167mlio.root");
    auto *tr=(TTree*)f->Get("a");

    auto h1st = new TH1F("h1st", "h1st", 1000, 0 , 4000);
    tr -> Draw("TPROTONS >> h1st");
    h1st ->SetLineColor(kRed);

    auto h2st = new TH1F("h2st", "h2st", 1000, 0 , 4000);
    tr -> Draw("TPROTONS >> h2st", "Edep1[0] > 0 && Edep3[1] > 0 && abs(Edep1[0]-Edep3[1]) < 400");
    h2st ->SetLineColor(kBlue);

    auto h3st = new TH1F("h3st", "h3st", 1000, 0 , 4000);
    //tr -> Draw("TPROTONS >> h3st", "Edep1[0] > 0 && Edep2[1] > 0 && abs(Edep1[0]-Edep2[1]) > 400");
    tr -> Draw("TPROTONS >> h3st", "Edep1[0] > 0 && Edep2[1] > 0");
    h3st ->SetLineColor(kGreen);

    auto h4st = new TH1F("h4st", "h4st", 1000, 0 , 4000);
    //tr -> Draw("TPROTONS >> h4st", "Edep1[0] > 0 && Edep1[1] > 0 && abs(Edep1[0]-Edep1[1]) > 400");
    tr -> Draw("TPROTONS >> h4st", "Edep1[0] > 0 && Edep1[1] > 0");
    h3st ->SetLineColor(kOrange);

     /*
    auto h4st = new TH1F("h4st", "h4st", 1000, 200 , 8000);
    tr -> Draw("Edep2 >> h4st");
    h4st ->SetLineColor(kMagenta);

    auto h5st = new TH1F("h5st", "h5st", 1000, 200 , 8000);
    tr -> Draw("Edep3 >> h5st");
    h5st ->SetLineColor(kOrange);
    */

    //Plotting every histogram together
    auto hs = new THStack("hs","TPROTONS");
    hs->Add(h1st, "hist");
    hs->Add(h2st, "hist");
    hs->Add(h3st, "hist");
    hs->Add(h4st, "hist");
    //hs->Add(h5st, "hist");

    auto l = new TLegend();

    l -> AddEntry(h1st, "Without gating", "L");
    l -> AddEntry(h2st, "8Li (~180 degrees)", "L");
    l -> AddEntry(h3st, "~90 degrees", "L");
    l -> AddEntry(h4st, "~0 degrees", "L");
    //l -> AddEntry(h5st, "Edep3", "L");

    hs->Draw("nostack");
    l->Draw();

    hs->SaveAs("TPROTONS.png");
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

