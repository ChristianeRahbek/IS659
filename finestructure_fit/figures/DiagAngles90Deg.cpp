//
// Created by Christiane Rahbek on 16-05-2023.
//

#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"

using namespace std;

void DiagAngles90Deg() {
    string filePath = "/mnt/d/IS659/finestructure_fit/analysis1/output/doubleDSSD/mergedRuns.root";
    string treeName = "a";
    string saveFilename = "diagAngles90Deg.pdf";

    string selectionCrit = "((id0 == 0 && (id1 == 3 || id1 == 1)) || (id0 == 1 && (id1 == 0 || id1 == 2)) || (id0 == 2 && (id1 == 1 || id1 == 3)) || (id0 == 3 && (id1 == 2 || id1 == 0))) && abs(FT0-FT1)<2500 && Edep0>2500 && Edep1>2500";

    auto canv = new TCanvas();
    auto *f= new TFile(filePath.c_str());
    auto *tr=(TTree*)f->Get(treeName.c_str());
    auto hist = new TH1F("hist", "hist", 181, 0, 180);
    tr->Draw("hitAng >> hist", selectionCrit.c_str());

    hist->SetXTitle("#theta");
    hist->SetYTitle("Entries");
    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetYaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetTitleSize(0.05);
    hist->SetTitle("");
    hist->SetTitleSize(0.06);
    hist->SetStats(kFALSE);
    hist->GetXaxis()->SetTitleOffset(0.75);
    hist->GetYaxis()->SetTitleOffset(0.8);

    gPad->SetTickx();
    gPad->SetTicky();

    canv->Update();
    canv->Draw();
    canv->SaveAs(saveFilename.c_str());
}