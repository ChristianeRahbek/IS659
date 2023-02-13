//
// Created by Chris on 13-02-2023.
//

#include <TFile.h>
#include <TChain.h>
#include <iostream>
#include <TH1F.h>

using namespace std;

int KolmogorovTest(){
    string filePath = "../../TPROTONS_id8_980.root";

    auto *f= new TFile(filePath.c_str());
    auto *hist=(TH1F*)f->Get("h");

    Double_t* sample = new Double_t[hist->GetEntries()];

    for(UInt_t i = 0; i < hist->GetEntries(); i++) {
        sample[i] = hist->GetBinContent(i);
    }

    return 0;
}