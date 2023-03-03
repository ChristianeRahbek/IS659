//
// Created by Christiane Rahbek on 27-02-2023.
//

#include <TFile.h>
#include <iostream>
#include <TTree.h>
#include <TVector3.h>
#include <TClonesArray.h>

using namespace std;
using namespace ROOT;

//void FindingCoincidentialHits() {
void BetaCoincidenceAnalysis() {
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
        if (mul < 3) continue; //we need at least three hits to have a coincidences with plastics as well

        for (int i = 0; i < mul; i++) { //looping through every instance in an entry
            auto FT0 = FT[i];
            auto id0 = id[i];

            for (int j = i + 1; j < mul; j++) {
                auto FT1 = FT[j];
                auto id1 = id[j];

                for (int k = j + 1; k < mul; k++) {
                    totalInstances++;
                    //cout << "ijk = " << i << j << k << endl;
                    auto FT2 = FT[k];
                    auto id2 = id[k];

                    //we want two of the instances in the same detector, and one in a pad


                    bool crit1 = !(id0 < 4 && id0 == id1 && abs(FT0 - FT1) < 1500 && id2 > 3);
                    bool crit2 = !(id0 < 4 && id0 == id2 && abs(FT0 - FT2) < 1500 && id1 > 3);
                    bool crit3 = !(id1 < 4 && id1 == id2 && abs(FT1 - FT2) < 1500 && id0 > 3);
                    /*
                    bool crit1 = !(id0 == id1);
                    bool crit2 = !(id0 == id2);
                    bool crit3 = !(id1 == id2);
                    */


                    if (crit1 && crit2 && crit3) continue;

                    TVector3 posa, posb;
                    double_t Edepa, Edepb;
                    if(id0==id1) {
                        Edepa = Edep[i];
                        Edepb = Edep[j];
                        posa = TVector3(posX[i], posY[i], posZ[i]);
                        posb = TVector3(posX[j], posY[j], posZ[j]);
                    }
                    else if(id0==id2) {
                        Edepa = Edep[i];
                        Edepb = Edep[k];
                        posa = TVector3(posX[i], posY[i], posZ[i]);
                        posb = TVector3(posX[k], posY[k], posZ[k]);
                    }
                    else {
                        Edepa = Edep[k];
                        Edepb = Edep[j];
                        posa = TVector3(posX[k], posY[k], posZ[k]);
                        posb = TVector3(posX[j], posY[j], posZ[j]);
                    }

                    if(Edepa >= 800 && Edepb >= 800) continue; //energies of less than 800 keV in one of the detectors is needed

                    auto ang = posa.Angle(posb)*TMath::RadToDeg();
                    if(ang > 10) continue; //we only want hits at 10 degrees or less.

                    sumAng += ang;

                    instancesLeft++;
                }
            }
        }
    }


    cout << "There are " << instancesLeft << " valid coincidences left out of " << totalInstances << " possible coincidences" << endl;
    cout << "Mean angle of the coincidences left was " << sumAng/instancesLeft << " degrees" << endl;
}