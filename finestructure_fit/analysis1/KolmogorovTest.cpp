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

void KolmogorovTest(){
    string filePath = "../../TPROTONS_id8_980.root";
    string myFilePath = "output/2coincidence/Run120mlio.root";
    int noOfBins = 1000;//7000;
    double xlow = 0;
    double xup = 500;//3500;

    auto canv = new TCanvas();
    canv->cd();

    auto *f= new TFile(filePath.c_str());
    auto *hist=(TH1F*)f->Get("h");

    auto madsHist = new TH1F("madshist", "madshist", noOfBins, xlow, xup);

    int madsEntries = 0;
    for(UInt_t i = 1; i <= hist->GetEntries(); i++) { //cut histogram to only the usable part
        if(hist->GetXaxis()->GetBinCenter(i) < xup) {
            auto data = hist->GetBinContent(i);
            madsHist->SetBinContent(i, hist->GetBinContent(i));
            madsEntries+=hist->GetBinContent(i);
        }
    }

    madsHist->Draw();

    auto *myfile = new TFile(myFilePath.c_str());
    auto *myTree=(TTree*)myfile->Get("a");

    auto myHist = new TH1F("myhist", "myhist", noOfBins, xlow, xup);
    //myTree->Draw("TPROTONS >> myhist", "id0==id1 && abs(FT0-FT1)<1500 && (Edep0<800 || Edep1<800)", "same");
    myTree->Draw("TPROTONS >> myhist", "((id0==0 && id1==2) || (id0==1 && id1==3) ||(id0==2 && id1==0) || (id0==3 && id1==1)) && abs(FT0-FT1)<1500 && hitAng<130 && abs(Edep0-Edep1)>500", "same");
/*
    auto pValLowAngle = madsHist->KolmogorovTest(myHist);

    cout << "p-val for <=10 degrees: " << pValLowAngle << endl;
*/
    auto canvAlt = new TCanvas();
    canvAlt->cd();

    madsHist->Draw();

    //auto myHistAltnative = new TH1F("myhistAlt", "myhistAlt", 7000,0,3500);
    //myTree->Draw("TPROTONS >> myhistAlt", "180/3.14159*acos((dir[0].X()*dir[1].X() + dir[0].Y()*dir[1].Y()+dir[0].Z()*dir[1].Z())/(pow(pow(dir[0].X(),2)+pow(dir[0].Y(),2)+pow(dir[0].Z(),2), 1/2)*pow(pow(dir[1].X(),2)+pow(dir[1].Y(),2)+pow(dir[1].Z(),2),1/2)))>10 && ((Edep0[0]>0 && Edep0[1]>0) || (Edep1[0]>0 && Edep1[1]>0) || (Edep2[0]>0 && Edep2[1]>0) ||(Edep3[0]>0 && Edep3[1]>0)) && (Edssd[0]<800 || Edssd[1]<800) && abs(FT[1]-FT[0])<1500", "same");
/*
    auto pValBigAngle = madsHist->KolmogorovTest(myHistAltnative);
    cout << "p-val for >10 degrees: " << pValBigAngle << endl;
*/


    Double_t* sampleMads = new Double_t[madsHist->GetEntries()];
    Double_t* sampleLowAng = new Double_t[myHist->GetEntries()];
    //Double_t* sampleBigAng = new Double_t[myHistAltnative->GetEntries()];

    for(UInt_t i = 1; i <= madsHist->GetEntries(); i++) {
        sampleMads[i-1] = madsHist->GetBinContent(i);
    }
    for(UInt_t i = 1; i <= myHist->GetEntries(); i++) {
        sampleLowAng[i-1] = myHist->GetBinContent(i);
    }
    /*
    for(UInt_t i = 1; i <= myHistAltnative->GetEntries(); i++){
        sampleBigAng[i-1] = myHistAltnative->GetBinContent(i);
    }
    */

    ROOT::Math::GoFTest* testLowAng = new ROOT::Math::GoFTest(madsEntries, sampleMads, myHist->GetEntries(), sampleLowAng);
    //ROOT::Math::GoFTest* testLowAng = new ROOT::Math::GoFTest(madsHist->GetEntries(), sampleMads, myHist->GetEntries(), sampleLowAng);
    //ROOT::Math::GoFTest* testBigAng = new ROOT::Math::GoFTest(madsHist->GetEntries(), sampleMads, myHistAltnative->GetEntries(), sampleBigAng);

    auto pValLow = testLowAng->KolmogorovSmirnov2SamplesTest();
    //auto pValBig = testBigAng->KolmogorovSmirnov2SamplesTest();

    cout << "pval low angles: " << pValLow << endl;
    //cout << "pval big angles: " << pValBig << endl;

    auto canvCumLow = new TCanvas();

    auto *cumMads = madsHist->DrawNormalized()->GetCumulative();
    auto *cumLow = myHist->DrawNormalized()->GetCumulative();
    //auto *cumBig = myHistAltnative->DrawNormalized()->GetCumulative();

    cumMads->Draw();
    cumLow->Draw("same");

    auto canvCum = new TCanvas();

    auto *cumMine = myHist->DrawNormalized()->GetCumulative();

    cumMads->Draw();
    cumMine->Draw("same");

    auto canvNormHists = new TCanvas();

    myHist->Scale(1/myHist->GetMaximum(), "nosw2");
    madsHist->Scale(1/madsHist->GetMaximum(), "nosw2");
    madsHist->Draw();
    myHist->Draw("same");

    //auto canvCumBig = new TCanvas();

    //cumMads->Draw();
    //cumBig->Draw("same");

}