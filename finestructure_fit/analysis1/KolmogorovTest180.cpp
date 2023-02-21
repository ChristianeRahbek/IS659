//
// Created by Chris on 13-02-2023.
//

#include <TFile.h>
#include <iostream>
#include <TH1F.h>
#include <TTree.h>
#include <TCanvas.h>
#include "Math/GoFTest.h"

using namespace std;

void KolmogorovTest180(){
    string filePath = "../../TPROTONS_id8_980.root";
    string myFilePath = "output/Run167mlio.root";

    auto canv = new TCanvas();
    canv->cd();

    auto *f= new TFile(filePath.c_str());
    auto *hist=(TH1F*)f->Get("h");

    auto madsHist = new TH1F("madshist", "madshist", 7000,0,3500);

    for(UInt_t i = 1; i <= hist->GetEntries(); i++) { //cut histogram to only the usable part
        if(hist->GetXaxis()->GetBinCenter(i) < 3500) {
            auto data = hist->GetBinContent(i);
            madsHist->SetBinContent(i, hist->GetBinContent(i));
        }
    }

    madsHist->Draw();

    auto *myfile = new TFile(myFilePath.c_str());
    auto *myTree=(TTree*)myfile->Get("a");

    auto myHist = new TH1F("myhist", "myhist", 7000,0,3500);
    myTree->Draw("TPROTONS >> myhist","180/3.14159*acos((dir[0].X()*dir[1].X() + dir[0].Y()*dir[1].Y()+dir[0].Z()*dir[1].Z())/(pow(pow(dir[0].X(),2)+pow(dir[0].Y(),2)+pow(dir[0].Z(),2), 1/2)*pow(pow(dir[1].X(),2)+pow(dir[1].Y(),2)+pow(dir[1].Z(),2),1/2)))<130 && ((Edep0[0]>0 && Edep2[1]>0) || (Edep2[0]>0 && Edep0[1]>0) || (Edep1[0]>0 && Edep3[1]>0) || (Edep3[0]>0 && Edep1[1]>0)) && abs(FT[1]-FT[0])<1000 && abs(Edssd[0]-Edssd[1])>500 && Edssd[0]>400 && Edssd[1]>400", "same");

    Double_t* sampleMads = new Double_t[madsHist->GetEntries()];
    Double_t* sampleMyAng = new Double_t[myHist->GetEntries()];

    for(UInt_t i = 1; i <= madsHist->GetEntries(); i++) {
        sampleMads[i-1] = madsHist->GetBinContent(i);
    }
    for(UInt_t i = 1; i <= myHist->GetEntries(); i++) {
        sampleMyAng[i-1] = myHist->GetBinContent(i);
    }

    ROOT::Math::GoFTest* testMine = new ROOT::Math::GoFTest(madsHist->GetEntries(), sampleMads, myHist->GetEntries(), sampleMyAng);

    auto pVal = testMine->KolmogorovSmirnov2SamplesTest();
    auto maxDist = testMine->KolmogorovSmirnov2SamplesTest("t");

    cout << "pval: " << to_string(pVal) << endl;
    cout << "max distance: " << to_string(maxDist) << endl;
    
    auto canvCum = new TCanvas();

    auto *cumMads = madsHist->DrawNormalized()->GetCumulative();
    auto *cumMine = myHist->DrawNormalized()->GetCumulative();

    cumMads->Draw();
    cumMine->Draw("same");

    auto canvNormHists = new TCanvas();

    myHist->Scale(1/myHist->GetMaximum(), "nosw2");
    madsHist->Scale(1/madsHist->GetMaximum(), "nosw2");
    myHist->Draw();
    madsHist->Draw("same");

}