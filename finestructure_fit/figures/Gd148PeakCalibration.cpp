//
// Created by Chris on 18-01-2023.
//

#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPad.h"
#include "TArrow.h"
#include <TROOT.h>
#include "TColor.h"
#include "TBox.h"

using namespace std;

void Gd148PeakCalibration() {
/* INSERT INFORMATION HERE*/
string filePath = "/mnt/d/IS659/finestructure_fit/calibrationAnalysis/output/Run212mlio.root";
string treeName = "a";
string branchName = "E";
string selectionCrit = "";
string title = "Run 212";
string xLabel = "Energy [keV]";
int noOfBins = 5500;
int xMin = 2500;
int xMax = 8000;

string filePath1 = "/mnt/d/IS659/finestructure_fit/calibrationAnalysis/output/Run32mlio.root";
string treeName1 = "a";
string branchName1 = "E";
string selectionCrit1 = "";
string title1 = "Run 32";
string xLabel1 = "Energy [keV]";
int noOfBins1 = 5500;
int xMin1 = 2500;
int xMax1 = 8000;

string saveFileName = "148GdPeakCalibration.pdf";

/* ********************** MAKING FIGURE ************************************************** */

auto canv = new TCanvas("", "", 1200, 500);

canv->Divide(2,1);
canv->cd(1);

/* First figure */
auto *f= new TFile(filePath.c_str());
auto *tr=(TTree*)f->Get(treeName.c_str());
auto hist = new TH1F("hist", title.c_str(), noOfBins, xMin, xMax);
tr->Draw((branchName +">> hist").c_str(), selectionCrit.c_str());

hist ->SetXTitle(xLabel.c_str());
hist->SetYTitle("Events");
canv->SetLeftMargin(0.1);
hist->GetYaxis()->SetTitleOffset(1);
hist->GetYaxis()->SetMaxDigits(3); //Setting numbers on axis to x*10^y
gStyle ->SetOptStat(kFALSE);
gStyle ->SetPalette(55);

//TArrow *ar = new TArrow(0.3,0.9, 0.95, 0.95, 0.05, "|>");
//ar->Draw();

//Making zoom on the graph in the corner
TPad *p = new TPad("p", "p", .58, .52, 0.95, 0.97); //Where the overlaying canvas (pad) is placed (numbers between 0 and 1)
p->Draw();
p->cd();
p->DrawFrame(3150,0,3220,40000);

auto hist01 = new TH1F("hist01", "hist01", noOfBins, xMin, xMax); //Making this so fit only shows up in the corner figure.
tr->Draw((branchName +">> hist01").c_str(), selectionCrit.c_str(), "Same");
hist01->Fit("gaus", "", "", 3175, 3195); //fitting
hist01->GetYaxis()->SetMaxDigits(3);

//hist01->GetYaxis()->ImportAttributes(hist->GetYaxis());
//p->RedrawAxis();


canv->cd(1);

TBox *box = new TBox(500, 500, 1000, 1000);
box->SetFillColorAlpha(0, 0.);
box->SetFillStyle(0);
box->SetLineWidth(2);
box->SetLineColor(1);
box->Draw();

//TPad *p01 = new TPad("p01", "p01", .15, .1, .20, .88);
//p01->SetFillColorAlpha(kWhite,0);
//p01->Draw();
//p01->cd();
//p01->PaintBox(.3, .1, .35, .8);
//p01->PaintPadFrame(0, 0, 1, 1);

canv->cd(2);

/* Second figure */

auto *f1= new TFile(filePath1.c_str());
auto *tr1=(TTree*)f1->Get(treeName1.c_str());

auto hist1 = new TH1F("hist1", title1.c_str(), noOfBins1, xMin1, xMax1);
tr1->Draw((branchName1 +">> hist1").c_str(), selectionCrit1.c_str());

hist1 ->SetXTitle(xLabel1.c_str());
hist1->SetYTitle("Events");
canv->SetLeftMargin(0.1);
hist1->GetYaxis()->SetTitleOffset(1);
hist1->GetYaxis()->SetMaxDigits(3); //Setting numbers on axis to x*10^y
gStyle ->SetOptStat(kFALSE);
gStyle ->SetPalette(55);

//Making zoom on the graph in the corner
TPad *p1 = new TPad("p1", "p1", .58, .52, 0.95, 0.97);
p1->Draw();
p1->cd();
p1->DrawFrame(3150,0,3220,18000);
auto hist11 = new TH1F("hist11", "hist11", noOfBins, xMin, xMax);
tr1->Draw((branchName +">> hist11").c_str(), selectionCrit.c_str(), "Same");
hist11->Fit("gaus", "", "", 3175, 3195);
hist11->GetYaxis()->SetMaxDigits(3);

canv->SaveAs(saveFileName.c_str());
}