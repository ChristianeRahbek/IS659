//
// Created by Christiane Rahbek on 02-03-2023
// This script is made to prove the concept of finding incidences in the same detector
//

#include <TFile.h>
#include <iostream>
#include <TTree.h>
#include <TVector3.h>
#include <TClonesArray.h>

using namespace std;
using namespace ROOT;

void proofOfConcept() {
    string filePath = "output/Run167mlio.root";

    auto *f = new TFile(filePath.c_str());
    auto *a = (TTree *) f->Get("a");

    UInt_t max_hits = 20;
    Int_t num;
    UInt_t mul, TPROTONS;
    UShort_t id[max_hits];
    Double_t FT[max_hits], Edep[max_hits], posX[max_hits], posY[max_hits], posZ[max_hits];
    //TVector3 pos[max_hits];
    TClonesArray pos("TVector3", max_hits);


    a->SetBranchAddress("num", &num);
    a->SetBranchAddress("id", &id);
    a->SetBranchAddress("mul", &mul);
    a->SetBranchAddress("FT", FT);
    a->SetBranchAddress("posX", posX);
    a->SetBranchAddress("posY", posY);
    a->SetBranchAddress("posY", posY);
    //a->SetBranchAddress("pos", pos);
    a->SetBranchAddress("Edep", Edep);

    int instancesLeft = 0;
    int totalInstances = 0;
    double sumAng = 0;

    auto entries = a->GetEntries();

    for (int ei = 0; ei < entries; ei++) {
        //cout << "Entry number is " << ei << endl;
        a->GetEntry(ei);
        if (mul < 2) continue;

        for (int i = 0; i < 2; i++) { //looping through every instance in an entry
            auto FT0 = FT[i];
            auto id0 = id[i];

            if(id0 > 3) continue;

            for (int j = i + 1; j < 2; j++) {
                totalInstances++;
                auto FT1 = FT[j];
                auto id1 = id[j];


                bool crit1 = !(id1 < 4 && id0 == id1 && abs(FT0 - FT1) < 1500);

                if (crit1) continue;

                TVector3 posa, posb;
                double_t Edepa, Edepb;

                Edepa = Edep[i];
                Edepb = Edep[j];
                posa = TVector3(posX[i], posY[i], posZ[i]);
                posb = TVector3(posX[j], posY[j], posZ[j]);

                if(Edepa >= 800 && Edepb >= 800) continue; //energies of less than 800 keV in one of the detectors is needed

                auto ang = posa.Angle(posb)*180/3.14159;//TMath::RadToDeg();
                if(ang > 10) continue; //we only want hits at 10 degrees or less.

                sumAng += ang;

                instancesLeft++;
            }
        }
    }


    cout << "There are " << instancesLeft << " valid coincidences left out of " << totalInstances << " possible coincidences" << endl;
    cout << "Mean angle of the coincidences left was " << sumAng/instancesLeft << " degrees" << endl;
}