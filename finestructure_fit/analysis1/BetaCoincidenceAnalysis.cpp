//
// Created by Christiane Rahbek on 27-02-2023.
//

#include <TFile.h>
#include <iostream>
#include <TTree.h>

using namespace std;

//void FindingCoincidentialHits() {
void BetaCoincidenceAnalysis() {
    string filePath = "output/Run167mlio.root";

    auto *f = new TFile(filePath.c_str());
    auto *a = (TTree *) f->Get("a");

    UInt_t max_hits = 20;
    Int_t num;
    UInt_t mul, TPROTONS;
    UShort_t id[max_hits];
    Double_t FT[max_hits];

    a->SetBranchAddress("num", &num);
    a->SetBranchAddress("id", &id);
    a->SetBranchAddress("mul", &mul);
    a->SetBranchAddress("FT", &FT);

    int instancesLeft = 0;
    int totalInstances = 0;

    auto entries = a->GetEntries();

    for (int ei = 0; ei < entries; ei++) {
        a->GetEntry(ei);
        if (mul < 3) continue; //we need at least three hits to have a coincidences with plastics as well

        for (int i = 0; i < mul; i++) {
            auto n0 = num;
            auto FT0 = FT[i];
            auto id0 = id[i];

            for (int j = i + 1; j < mul; j++) {
                auto n1 = num;
                auto FT1 = FT[j];
                auto id1 = id[j];

                for (int k = j + 1; k < mul; k++) {
                    totalInstances++;
                    auto n2 = num;
                    auto FT2 = FT[k];
                    auto id2 = id[k];

                    //we want two of the instances in the same detector, and one in a pad


                    bool crit1 = !(id0 < 4 && id0 == id1 && abs(FT0 - FT1) < 1000 && id2 > 3);
                    bool crit2 = !(id0 < 4 && id0 == id2 && abs(FT0 - FT2) < 1000 && id1 > 3);
                    bool crit3 = !(id1 < 4 && id1 == id2 && abs(FT1 - FT2) < 1000 && id0 > 3);
                    /*
                    bool crit1 = !(id0 == id1);
                    bool crit2 = !(id0 == id2);
                    bool crit3 = !(id1 == id2);
                    */
                    if (crit1 && crit2 && crit3) continue;


                    //if (id0 > 3) continue;

                    instancesLeft++;
                }
            }
        }
    }


    cout << "There are " << instancesLeft << " instances left out of " << totalInstances << endl;
}