//
// Created by erik on 6/10/22.
// Split branches shared between all detectors from output of bat.cpp into separate detector branches.
//

#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <ausa/util/DynamicBranchVector.h>
#include <iomanip>
#include "projectutil.h"

using namespace std;
using namespace AUSA;

void bat(int run_num) {
  ostringstream run_num_padded;
  run_num_padded << setw(3) << setfill('0') << run_num;
  string file_name = "Run" + run_num_padded.str() + "mlio.root";
  auto c = new TChain("a");
  // possible to give file (char* path) as relative or absolute path, but it MUST be located in data/bpa
  c->Add((EUtil::getProjectRoot() + "/data/bat/" + EUtil::getBasename(file_name)).c_str());

  UInt_t max_hits = 10;
  Int_t num;
  UInt_t mul, TPROTONS;
  UShort_t id[max_hits];
  Double_t Edep[max_hits];

  c->SetBranchAddress("num", &num);
  c->SetBranchAddress("mul", &mul);
  c->SetBranchAddress("TPROTONS", &TPROTONS);
  c->SetBranchAddress("id", id);
  c->SetBranchAddress("Edep", Edep);

  system(("mkdir -p " + EUtil::getProjectRoot() + "/data/bat_split/").c_str());
  auto out = new TFile((EUtil::getProjectRoot() + "/data/bat_split/" + EUtil::getBasename(file_name)).c_str(), "RECREATE");
  auto tree = new TTree("a", "a");

  tree->Branch("num", &num);
  tree->Branch("mul", &mul);
  tree->Branch("TPROTONS", &TPROTONS);
  auto idout = make_unique<DynamicBranchVector<unsigned short>>(*tree, "id", "mul");
  auto Edep0 = make_unique<DynamicBranchVector<double>>(*tree, "Edep0", "mul");
  auto Edep1 = make_unique<DynamicBranchVector<double>>(*tree, "Edep1", "mul");
  auto Edep2 = make_unique<DynamicBranchVector<double>>(*tree, "Edep2", "mul");
  auto Edep3 = make_unique<DynamicBranchVector<double>>(*tree, "Edep3", "mul");

  for (UInt_t i = 0; i < c->GetEntries(); i++) {
    AUSA::clear(*idout, *Edep0, *Edep1, *Edep2, *Edep3);
    c->GetEntry(i);

    for (UInt_t j = 0; j < mul; j++) {
      idout->add(id[j]);
      if (id[j] == 0) {
        Edep0->add(Edep[j]);
      } else if (id[j] == 1) {
        Edep1->add(Edep[j]);
      } else if (id[j] == 2) {
        Edep2->add(Edep[j]);
      } else {
        Edep3->add(Edep[j]);
      }
    }
    tree->Fill();
  }

  tree->Write();
  out->Close();
}

int main(int argc, char* argv[]) {
  if (argc == 2) bat(stoi(argv[1]));
  return EXIT_SUCCESS;
}