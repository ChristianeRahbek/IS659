//
// Created by Chris on 13-02-2023.
//

#include <TFile.h>
#include <iostream>
#include <TH1F.h>
#include <TTree.h>

using namespace std;

Double_t KolmogorovTest(){
    string filePath = "../../TPROTONS_id8_980.root";
    string myFilePath = "output/Run167mlio.root";

    auto *f= new TFile(filePath.c_str());
    auto *hist=(TH1F*)f->Get("h");

    Double_t* sample = new Double_t[hist->GetEntries()];

    auto madsHist = new TH1F("madshist", "madshist", 1000,0,3500);
    /*
    for(UInt_t i = 1; i <= hist->GetEntries(); i++) { //cut histogram to only the usable part
        if(hist->GetXaxis()->GetBinCenter(i) < 3500) {
            auto data = hist->GetBinContent(i);
            sample[i-1] = data;
            madsHist->Fill(data);
        }
    }
     madsHist->Draw();
     */
    for(UInt_t i = 1; i <= hist->GetEntries(); i++) { //cut histogram to only the usable part
        if(hist->GetXaxis()->GetBinCenter(i) >= 3500) {
            hist->SetBinContent(i, 0);
        }
    }

    hist->Draw();


    auto *myfile = new TFile(myFilePath.c_str());
    auto *myTree=(TTree*)myfile->Get("a");

    auto myHist = new TH1F("myhist", "myhist", 1000,0,3500);
    myTree->Draw("TPROTONS >> myhist", "180/3.14159*acos((dir[0].X()*dir[1].X() + dir[0].Y()*dir[1].Y()+dir[0].Z()*dir[1].Z())/(pow(pow(dir[0].X(),2)+pow(dir[0].Y(),2)+pow(dir[0].Z(),2), 1/2)*pow(pow(dir[1].X(),2)+pow(dir[1].Y(),2)+pow(dir[1].Z(),2),1/2)))<=10 && ((Edep0[0]>0 && Edep0[1]>0) || (Edep1[0]>0 && Edep1[1]>0) || (Edep2[0]>0 && Edep2[1]>0) ||(Edep3[0]>0 && Edep3[1]>0)) && (Edssd[0]<800 || Edssd[1]<800) && abs(FT[1]-FT[0])<1500", "same");

    auto madsHistNormalized = hist->DrawNormalized();
    auto myHistNormalized = myHist->DrawNormalized();

    auto maxKolmogorovDistance = myHistNormalized->KolmogorovTest(madsHistNormalized, "N");

    return maxKolmogorovDistance;
}