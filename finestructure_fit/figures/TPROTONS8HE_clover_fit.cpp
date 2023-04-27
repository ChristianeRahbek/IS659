//
// Created by Christiane Rahbek on 21-04-2023.
// Finds the norm for each clover from the fitting done in TPROTONS8He_fit.cpp
//

#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TGraphErrors.h"

using namespace std;

void TPROTONS8HE_clover_fit() {
    string filePath, histName;
    double BG_del, BG_open;
    filePath = "/mnt/d/IS659/christiane_eff/spectra_christiane_R121.root";
    histName = "C";
    BG_del = 10;
    BG_open = 300;

    string title = "TPROTONS fit";
    string xLabel = "T [ms]";
    int xMin = 0;
    int xMax = 3500;

    string saveFileName = "TPROTONS_fit.pdf";

    /* ********************** Making fit **************************** */

    auto canv = new TCanvas("", "", 1000, 800);

    canv->Divide(2,2);


    /* Figure */
    auto *f = new TFile(filePath.c_str());
    auto *hist1 = (TH1F *) f->Get((histName+to_string(1)).c_str());
    auto *hist2 = (TH1F *) f->Get((histName+to_string(2)).c_str());
    auto *hist3 = (TH1F *) f->Get((histName+to_string(3)).c_str());
    auto *hist4 = (TH1F *) f->Get((histName+to_string(4)).c_str());

    hist1->Rebin(3);
    hist2->Rebin(3);
    hist3->Rebin(3);
    hist4->Rebin(4);

    canv->cd(1); hist1->Draw();
    canv->cd(2); hist2->Draw();
    canv->cd(3); hist3->Draw();
    canv->cd(4); hist4->Draw();

    //hist1->GetXaxis()->SetRangeUser(xMin, xMax);

    //double lambda8He = 1 / 119.1; // 1/halflife of 8He
    //double lambda8Li = 1 / 839.9; // 1/halflife of 8Li

    string formula;
    formula = "[0]*exp(-log(2)*(x+[3])/119.1)*(1/[1]*(1-exp(-[1]*(x+[3]))) + 1/([2]+[1])*(exp(-([2]+[1])*(x+[3]))-1)) + [4]*x + [5]";

    auto fitFunc = new TF1("fitFunc", formula.c_str(), xMin, xMax);
    fitFunc->SetParameter(0, 100);
    fitFunc->SetParName(0, "N");
    fitFunc->SetParNames("N", "lambda_r", "lambda_f", "t'", "a", "b");
    fitFunc->FixParameter(1,1.01127e-02);
    fitFunc->FixParameter(2, 7.44872e-02);
    fitFunc->FixParameter(3,-1.4410e+01);
    fitFunc->FixParameter(4, 4.47505e-02);
    fitFunc->FixParameter(5, 1.01994e+00);


    //auto effFunc = new TF1("effFunc", "[3] * (TMath::Exp( [0]*TMath::Log(x/[4]) + [1]*TMath::Power(TMath::Log(x/[4]), 2) - [2]/TMath::Power(x,3) ) )");
    auto effFunc = new TF1("effFunc", "[3] * (exp( [0]*log(x/[4]) + [1]*pow(log(x/[4]), 2) - [2]/pow(x,3) ) )");
    auto eff_d0 = new TF1("d0", "log(x/1000) * effFunc");
    auto eff_d1 = new TF1("d1", "pow(log(x/1000), 2) * effFunc");
    auto eff_d2 = new TF1("d2", "-1/pow(x,3) * effFunc");
    auto eff_d3 = new TF1("d3", "1/[0] * effFunc");
    auto effFuncErr = new TF1("effFuncErr", "pow([3]*d3,2) + pow([2]*d2,2) + pow([1]*d1,2) + pow([0]*d0,2) + [4]*d0*d1 + [5]*d0*d2 + [6]*d0*d3 + [7]*d1*d2 + [8]*d1*d3 + [9]*d2*d3");
    auto effErr = new TF1("effErr", "sqrt(pow([0]*effFunc,2)*(pow(effFuncErr/effFunc,2) + pow([0]/[1],2)))");
            //REMEMBER TO SET PARAMS FOR effFuncErr

    cout << "Fitting C1:" << endl;
    hist1->Fit("fitFunc", "", "", xMin + BG_del, BG_open + BG_del);
    effFunc->SetParameters(-0.656915, 0.0453174, 761991, 0.000172134, 1000); //setting parameters given from Mads' calibration
    effFuncErr->SetParameters(0.040494, 0.0423398, 209288  , 1.79534e-06, 0.001564, 6196, 1.113e-08, 8230, 2.732e-08, 0.1435 );
    eff_d3->SetParameter(0, effFunc->GetParameter(3));
    effErr->SetParameters(fitFunc->GetParameter(0),fitFunc->GetParError(0));
    auto eff1 = fitFunc->GetParameter(0)*effFunc->Eval(980);
    auto effErr1 = effErr->Eval(980);
    cout << "N1_eff = " << eff1 << " +- " << effErr1 << endl;

    cout << "Fitting C2:" << endl;
    hist2->Fit("fitFunc", "", "", xMin + BG_del, BG_open + BG_del);
    effFunc->SetParameters(-0.693445, 0.0303209, 726016, 0.000172092, 1000); //setting parameters given from Mads' calibration
    effFuncErr->SetParameters(0.0418139 , 0.0435551, 213236, 1.83373e-06,0.001672, 6651, 1.113e-08, 8657, 2.806e-08, 0.148);
    eff_d3->SetParameter(0, effFunc->GetParameter(3));
    effErr->SetParameters(fitFunc->GetParameter(0),fitFunc->GetParError(0));
    auto eff2 = fitFunc->GetParameter(0)*effFunc->Eval(980);
    auto effErr2 = effErr->Eval(980);
    cout << "N2_eff = " << eff2 << " +- " << effErr2 << endl;

    cout << "Fitting C3:" << endl;
    hist3->Fit("fitFunc", "", "", xMin + BG_del, BG_open + BG_del);
    effFunc->SetParameters(-0.675181, 0.0356154, 539443, 0.000168379, 1000); //setting parameters given from Mads' calibration
    effFuncErr->SetParameters(0.0408049, 0.0430538, 213215, 1.81564e-06,0.001602, 6390, 1.22e-08, 8548, 2.93e-08, 0.1538);
    eff_d3->SetParameter(0, effFunc->GetParameter(3));
    effErr->SetParameters(fitFunc->GetParameter(0),fitFunc->GetParError(0));
    auto eff3 = fitFunc->GetParameter(0)*effFunc->Eval(980);
    auto effErr3 = effErr->Eval(980);
    cout << "N3_eff = " << eff3 << " +- " << effErr3 << endl;

    /*
     * MADE IT TO HERE IN THE CORRECTION!!!
     * */

    cout << "Fitting C4:" << endl;
    hist4->Fit("fitFunc", "", "", xMin + BG_del, BG_open + BG_del);
    effFunc->SetParameters(-0.65969, 0.0459889, 595718, 0.000174134, 1000); //setting parameters given from Mads' calibration
    effFuncErr->SetParameters(0.0403567, 0.0427729, 212226, 1.82149e-06,0.001577, 6326,1.261e-08, 8460,2.94e-08, 0.154);
    eff_d3->SetParameter(0, effFunc->GetParameter(3));
    effErr->SetParameters(fitFunc->GetParameter(0),fitFunc->GetParError(0));
    auto eff4 = fitFunc->GetParameter(0)*effFunc->Eval(980);
    auto effErr4 = effErr->Eval(980);
    cout << "N4_eff = " << eff4 << " +- " << effErr4 << endl;

    cout << "number of bins fitted = " << hist1->GetXaxis()->FindBin(BG_open+BG_del) - hist1->GetXaxis()->FindBin(xMin+BG_del) << endl;

    double effs[4] = {eff1,eff2,eff3,eff4};
    double effErrs[4] = {effErr1,effErr2,effErr3,effErr4};
    double clovers[4] = {1,2,3,4};

    auto effCanv = new TCanvas("", "", 1000, 800);

    //TGraph *effGraph = new TGraph(4,clovers, effs);
    TGraph *effGraph = new TGraphErrors(4,clovers, effs, 0, effErrs);
    effGraph->SetTitle("Efficiency for each clover");
    effGraph->SetMarkerStyle(3);
    effGraph->Draw("ap");

    double effMean = TMath::Mean(4, effs);
    double meanXs[2] = {0,5};
    double means[2] = {effMean,effMean};
    auto *meanGraph = new TGraph(2, meanXs, means);
    meanGraph->SetLineStyle(2);
    meanGraph->SetTitle("Mean");

    meanGraph->Draw("Same");

    effCanv->BuildLegend();
    effCanv->SetTitle("Efficiency corregated value of N");

    effCanv->Update(); effCanv->Draw();
}