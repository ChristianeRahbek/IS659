//
// Made by Christiane 23/01-2023
//

#include "THStack.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"

void DetTimeCompare() {
    //defining each histogram

    std::string runNo = "32";
    std::string index = "1";
    bool frontIndex = true;
    bool noIndex = false; //if false we only look in a single strip either in front or back depending on frontIndex variable above
    std::string selExtra;
    if(frontIndex){
        selExtra = " && FI == " + index;
    }
    else selExtra = " && BI == " + index;

    if(noIndex) selExtra = "";

    auto canv = new TCanvas("", "", 700, 700);
    canv->cd()->SetLogy();

    auto *f= new TFile(("/mnt/d/IS659/finestructure_fit/analysis1/output/Run" + runNo + "mlio.root").c_str());
    auto *tr=(TTree*)f->Get("a");

    auto h1st = new TH1F("h1st", "h1st", 500, -10 , 50000);
    tr -> Draw("abs(FT[0]-BT[0]) >> h1st", ("id == 0" + selExtra).c_str());
    h1st ->SetLineColor(kRed);

    auto h2st = new TH1F("h2st", "h2st", 500, -10 , 50000);
    tr -> Draw("abs(FT[0]-BT[0]) >> h2st", ("id == 1" + selExtra).c_str());
    h2st ->SetLineColor(kBlue);

    auto h3st = new TH1F("h3st", "h3st", 500, -10 , 50000);
    tr -> Draw("abs(FT[0]-BT[0]) >> h3st", ("id == 2" + selExtra).c_str());
    h3st ->SetLineColor(kGreen);

    auto h4st = new TH1F("h4st", "h4st", 500, -10 , 50000);
    tr -> Draw("abs(FT[0]-BT[0]) >> h4st", ("id == 3" + selExtra).c_str());
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
    auto hs = new THStack("hs",("Detector front/back time comparison, run " + runNo).c_str());
    hs->Add(h1st, "hist");
    hs->Add(h2st, "hist");
    hs->Add(h3st, "hist");
    hs->Add(h4st, "hist");
    //hs->Add(h5st, "hist");

    auto l = new TLegend();

    l -> AddEntry(h1st, "Id = 0", "L");
    l -> AddEntry(h2st, "Id = 1", "L");
    l -> AddEntry(h3st, "Id = 2", "L");
    l -> AddEntry(h4st, "Id = 3", "L");
    //l -> AddEntry(h5st, "Edep3", "L");

    hs->Draw("nostack");
    l->Draw();

    if(selExtra == "") {
        canv->SaveAs(("DetTimeCompare/DetTimeCompareRun" + runNo + ".png").c_str());
    }
    else if (frontIndex) {
        canv->SaveAs(("DetTimeCompare/DetTimeCompareRun" + runNo + "FI" + index + ".png").c_str());
    }
    else canv->SaveAs(("DetTimeCompare/DetTimeCompareRun" + runNo + "BI" + index + ".png").c_str());

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

