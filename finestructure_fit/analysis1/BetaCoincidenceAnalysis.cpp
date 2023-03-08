//
// Created by Christiane Rahbek on 27-02-2023.
//

#include <TFile.h>
#include <iostream>
#include <TTree.h>
#include <TVector3.h>
#include <TClonesArray.h>
#include <ausa/util/DynamicBranchVector.h>
#include "projectutil.h"

using namespace std;
using namespace ROOT;
using namespace AUSA;
//using namespace libconfig;

//void FindingCoincidentialHits() {

class BetaCoincidenceAnalysis {
public:
    BetaCoincidenceAnalysis() {
        string filePath = EUtil::getProjectRoot() + "/output/Run167mlio.root";

        hasBeta = false;

        auto *f = new TFile(filePath.c_str());
        a = (TTree *) f->Get("a");

        a->SetBranchAddress("num", &num);
        a->SetBranchAddress("id", &id);
        a->SetBranchAddress("mul", &mul);
        a->SetBranchAddress("FT", FT);
        a->SetBranchAddress("posX", posX);
        a->SetBranchAddress("posY", posY);
        a->SetBranchAddress("posY", posY);
        //a->SetBranchAddress("pos", pos);
        a->SetBranchAddress("Edep", Edep);

        instancesLeft = 0;
        totalInstances = 0;

        // Defing Ttree and branches for the output tree
        string fileName = "90DegDets_10Deg.root";
        out = new TFile((EUtil::getProjectRoot() + "/output/betaAnalysis/" + fileName).c_str(), "RECREATE");
        tree = new TTree("a", "a");

        tree->Branch("mul", &multi);
        tree->Branch("hasBeta", &hasBeta);
        Ea = make_unique<DynamicBranchVector<double>>(*tree, "Ea", "mul");
        Eb = make_unique<DynamicBranchVector<double>>(*tree, "Eb", "mul");
        EPlastic = make_unique<DynamicBranchVector<double>>(*tree, "EPlastic", "mul");
        timeA = make_unique<DynamicBranchVector<double>>(*tree, "timeA", "mul");
        timeB = make_unique<DynamicBranchVector<double>>(*tree, "timeB", "mul");
        timeC = make_unique<DynamicBranchVector<double>>(*tree, "timeB", "mul");
        //betaBool = make_unique<DynamicBranchVector<Bool_t>>(*tree, "hasBeta", "mul");

        entries = a->GetEntries();

    }

    void Analyze() {
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
                        multi = 0;
                        hasBeta = false;
                        AUSA::clear(*Ea, *Eb, *EPlastic,*timeA, *timeB, *timeC);

                        totalInstances++;
                        //cout << "ijk = " << i << j << k << endl;
                        auto FT2 = FT[k];
                        auto id2 = id[k];

                        if(FT0 > 3e9 || FT1 > 3e9 || FT2 > 3e9) continue;
                        //we want two of the instances in the same detector, and one in a pad
                        bool crit1 = !(id0 < 4 && id0 == id1 && abs(FT0 - FT1) < 1500 && id2 > 3);
                        bool crit2 = !(id0 < 4 && id0 == id2 && abs(FT0 - FT2) < 1500 && id1 > 3);
                        bool crit3 = !(id1 < 4 && id1 == id2 && abs(FT1 - FT2) < 1500 && id0 > 3);

                        if (crit1 && crit2 && crit3) hasBeta = true;
                        //now making sure we at least have two coincidences in the same DSSD
                        else if (!((id0 < 4 && id0 == id1) || (id0 < 4 && id0 == id2) || (id1 < 4 && id1 == id2))) continue;
                        totalInstances++;

                        if(id0==id1) {
                            SetValues(i,j,k);
                        }
                        else if(id0==id2) {
                            SetValues(i,k,j);
                        }
                        else {
                            SetValues(k,j, i);
                        }

                        if(Edepa >= 800 && Edepb >= 800) continue; //energies of less than 800 keV in one of the detectors is needed

                        auto ang = posa.Angle(posb)*TMath::RadToDeg();
                        if(ang > 10) continue; //we only want hits at 10 degrees or less.
                        //if(abs(tA-tC)>1500 || abs(tB-tC)>1500) continue; // only want coincidential hits...

                        Ea->add(Edepa);
                        Eb->add(Edepb);
                        EPlastic->add(EdepPlastic);
                        timeA->add(tA);
                        timeB->add(tB);
                        timeC->add(tC);

                        multi++;
                        tree->Fill();
                        instancesLeft++;
                    }
                }
            }
        }

        tree->Write();
        out->Close();

        cout << "There are " << instancesLeft << " valid coincidences left out of " << totalInstances << " possible coincidences" << endl;

    }

    void SetValues(int i, int j, int k) {
        Edepa = Edep[i];
        Edepb = Edep[j];
        EdepPlastic = NAN;
        posa = TVector3(posX[i], posY[i], posZ[i]);
        posb = TVector3(posX[j], posY[j], posZ[j]);
        tA = FT[i];
        tB = FT[j];
        tC = FT[k];

    }

    UInt_t mul, TPROTONS;
    Int_t num;
    UShort_t id[20];
    Double_t FT[20], Edep[20], posX[20], posY[20], posZ[20];

    unique_ptr<DynamicBranchVector<double>> Ea, Eb, EPlastic, timeA, timeB, timeC;
    //unique_ptr<DynamicBranchVector<Bool_t>> betaBool;
    TFile *out;
    TTree *a, *tree;
    int instancesLeft, totalInstances, multi; //multi is multiplicity variable in the output tree

    TVector3 posa, posb;
    double_t Edepa, Edepb, EdepPlastic, tA, tB, tC;
    Bool_t hasBeta;

    Long64_t entries;
};

int main() {
    auto analysis = new BetaCoincidenceAnalysis();

    analysis->Analyze();

    return EXIT_SUCCESS;
}