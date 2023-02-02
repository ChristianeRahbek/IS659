//
// Created by Christiane Rahbek on 02-02-2023.
//

#include <iostream>
#include <string>
#include <TROOT.h>
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include <ausa/util/DynamicBranchVector.h>
#include "TCanvas.h"
#include "TH1F.h"
//#include "projectutil.h"

using namespace std;
using namespace AUSA;

void AngleAnalysis() {
    cout << "Opening document and making tree..." << endl;
    //string filePath = EUtil::getProjectRoot() + "output/Run167mlio.root";
    string filePath = "output/Run167mlio.root";
    string treeName = "a";

    auto *f= new TFile(filePath.c_str());
    auto *tr=(TTree*)f->Get(treeName.c_str());

    cout << "Done, now initialising branches..." << endl;

    int det_a = 0;
    int det_b = 1;

    UInt_t max_hits = 10;
    //Branches I would like to access
    UInt_t mul, TPROTONS;
    TVector3 dir[max_hits];
    //TClonesArray dir;
    Short_t id[max_hits]; //list with length max_hits
    Double_t FT[max_hits];

    cout << "Setting branch adresses..." << endl;

    //Referencing branches, so I can access in them in this code
    tr->SetBranchAddress("mul", &mul);
    tr->SetBranchAddress("TPROTONS", &TPROTONS);
    tr->SetBranchAddress("dir", &dir);
    tr->SetBranchAddress("id", id);
    tr->SetBranchAddress("FT", FT);

    //looping over all entries

    cout << "Making new tree" << endl;
    auto tr1 = new TTree("a1", "a1");
    //auto ang1 = make_unique<DynamicBranchVector<double>>(*tr1, "ang1", "mul");
    cout << "Making new array" << endl;
    cout << tr->GetEntries() << endl;
    //Double_t ang1[tr->GetEntries()];
    Double_t ang1[max_hits];

    cout << "Done, now looping through entries" << endl;
    for(UInt_t i = 0; i < tr->GetEntries(); i++) {
        tr->GetEntry(i);

        //looping over every instance in an entry
        for(UInt_t j = 0; j < mul; j++){
            TVector3 dir0, dir1;
            if(mul < 2) break; //only do the following if we have at least 2 instances
            if(j == 0 && id[j] == 0){
                dir0 = dir[j];
            }
            else if(j==1 && id[j]==1){
                dir1 = dir[j];
            }
            else continue; //only looking at instance 0 in det 1 and instance 1 in det 2

            //calculating angle between the vectors
            auto ang = dir1.Angle(dir0);

            //deciding what angles should be saved
            if (TPROTONS > 1000 && abs(FT[0]-FT[1])<2500) {
                ang1[i] = ang;
            }
        }
    }
    cout << "Done, making canvas..." << endl;
    /* ********************** MAKING FIGURE ************************************************** */
    auto canv = new TCanvas("", "", 1200, 500);


    auto hist = new TH1F("hist", "", 100, 0, 360);
    for(int i = 0; i < tr->GetEntries(); i++){
        hist->Fill(ang1[i]);
    }
    hist->Draw();


    /*
    for(auto i = tr->GetMinimum("id"); i < tr->GetMaximum("id")+1; i++) {
        for(auto j = i; j < tr->GetMaximum("id")+1; j++){

        }
    }
    */
}