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

    auto effFunc = new TF1("effFunc", "[3] * (exp( [0]*log(x/[4]) + [1]*pow(log(x/[4]), 2) - [2]/pow(x,3) ) )",1,3000);
    effFunc->SetParameters(-0.656915, 0.0453174, 761991, 0.000172134, 1000); //setting parameters given from Mads' calibration
    /*
    auto eff_d0 = new TF1("d0", "log(x/1000) * effFunc",1,3000);
    auto eff_d1 = new TF1("d1", "pow(log(x/1000), 2) * effFunc",1,3000);
    auto eff_d2 = new TF1("d2", "-1/pow(x,3) * effFunc",1,3000);
    auto eff_d3 = new TF1("d3", "1/[0] * effFunc",1,3000);
    eff_d3->SetParameter(0, effFunc->GetParameter(3));
    cout << eff_d3->GetParameter(0) << " & " << eff_d3->Eval(980)<< endl;
    auto effFuncErr = new TF1("effFuncErr", "pow([3]*d3,2) + pow([2]*d2,2) + pow([1]*d1,2) + pow([0]*d0,2) + [4]*d0*d1 + [5]*d0*d2 + [6]*d0*d3 + [7]*d1*d2 + [8]*d1*d3 + [9]*d2*d3",1,3000);
    effFuncErr->SetParameters(0.040494, 0.0423398, 209288  , 1.79534e-06, 0.001564, 6196, 1.113e-08, 8230, 2.732e-08, 0.1435);
    */
    auto effErr = new TF1("effErr", "sqrt(pow([0]*effFunc,2)*(pow(effFuncErr/effFunc,2) + pow([0]/[1],2)))",1,3000);
    effErr->SetParameters(fitFunc->GetParameter(0),fitFunc->GetParError(0));
    string effFunc_U1 = "0.000172134 * (exp(-0.656915*log(x/1000) + 0.0453174*pow(log(x/1000), 2) - 761991/pow(x,3)))";
    string eff_d0_U1 = "log(x/1000) * (" + effFunc_U1 + ")";
    string eff_d1_U1 = "pow(log(x/1000), 2) * (" + effFunc_U1 + ")";
    string eff_d2_U1 = "-1/pow(x,3) * (" + effFunc_U1 + ")";
    string eff_d3_U1 = "1/0.000172134 * (" + effFunc_U1 + ")";
    string effFuncErr_U1 = "sqrt(pow((1.79534e-06)*(" + eff_d3_U1 + "),2) + pow(209288*(" + eff_d2_U1 + "),2) + pow(0.0423398*(" + eff_d1_U1 + "),2) + pow(0.040494*(" + eff_d0_U1 + "),2) + 0.001564*(" + eff_d0_U1 + ")*(" + eff_d1_U1 + ") + 6196*(" + eff_d0_U1 + ")*(" + eff_d2_U1 + ") - 1.113e-08*(" + eff_d0_U1 + ")*(" + eff_d3_U1 + ") + 8230*(" + eff_d1_U1 + ")*(" + eff_d2_U1 + ") - 2.732e-08*(" + eff_d1_U1 + ")*(" + eff_d3_U1 + ") - 0.1435*(" + eff_d2_U1 + ")*(" + eff_d3_U1 + "))";
    auto effFuncErrU1 = new TF1("effFuncErrU1", effFuncErr_U1.c_str());
    auto effFuncErr980U1 = effFuncErrU1->Eval(980);
    auto effFunc980U1 = effFunc->Eval(980);
    string effErr_U1 = "sqrt(pow((9.08513e-01)/(" + effFunc_U1 + "),2)*(pow((" + effFuncErr_U1 + ")/(" + effFunc_U1 + "),2) + pow((9.08513e-01)/(2.48847e-02),2)))";
    auto effErrU1 = new TF1("effErrU1", effErr_U1.c_str());

    cout << "Fitting C1:" << endl;
    hist1->Fit("fitFunc", "", "", xMin + BG_del, BG_open + BG_del);
    auto NU1 = fitFunc->GetParameter(0);
    auto NErrU1 = fitFunc->GetParError(0);
    auto eff1 = fitFunc->GetParameter(0)/effFunc->Eval(980);
    auto effErr1 = abs(NU1/effFunc980U1)*sqrt(pow(effFuncErr980U1/effFunc980U1,2)+pow(NErrU1/NU1,2));
    cout << "N1_eff = " << eff1 << " +- " << effErr1 << endl;

    // NEXT DETECTOR
    string effFunc_U2 = "0.000172092 * (exp(-0.693445*log(x/1000) + 0.0303209*pow(log(x/1000), 2) - 726016/pow(x,3)))";
    string eff_d0_U2 = "log(x/1000) * (" + effFunc_U2 + ")";
    string eff_d1_U2 = "pow(log(x/1000), 2) * (" + effFunc_U2 + ")";
    string eff_d2_U2 = "-1/pow(x,3) * (" + effFunc_U2 + ")";
    string eff_d3_U2 = "1/0.000172092 * (" + effFunc_U2 + ")";
    string effFuncErr_U2 = "pow(1.83373e-06*(" + eff_d3_U2 + "),2) + pow(213236*(" + eff_d2_U2 + "),2) + pow(0.0435551*(" + eff_d1_U2 + "),2) + pow(0.0418139*(" + eff_d0_U2 + "),2) + 0.001672*(" + eff_d0_U2 + ")*(" + eff_d1_U2 + ") + 6651*(" + eff_d0_U2 + ")*(" + eff_d2_U2 + ") - 1.113e-08*(" + eff_d0_U2 + ")*(" + eff_d3_U2 + ") + 8657*(" + eff_d1_U2 + ")*(" + eff_d2_U2 + ") - 2.806e-08*(" + eff_d1_U2 + ")*(" + eff_d3_U2 + ") - 0.148*(" + eff_d2_U2 + ")*(" + eff_d3_U2 + ")";
    string effErr_U2 = "sqrt(pow((5.93017e-01)/(" + effFunc_U2 + "),2)*(pow((" + effFuncErr_U2 + ")/(" + effFunc_U2 + "),2) + pow(5.93017e-01/(2.15680e-02),2)))";
    auto effErrU2 = new TF1("effErrU2", effErr_U2.c_str());
    cout << "Fitting C2:" << endl;
    hist2->Fit("fitFunc", "", "", xMin + BG_del, BG_open + BG_del);
    effFunc->SetParameters(-0.693445, 0.0303209, 726016, 0.000172092, 1000); //setting parameters given from Mads' calibration
    auto effFuncErrU2 = new TF1("effFuncErrU2", effFuncErr_U2.c_str());
    auto effFuncErr980U2 = effFuncErrU2->Eval(980);
    auto effFunc980U2 = effFunc->Eval(980);
    auto NU2 = fitFunc->GetParameter(0);
    auto NErrU2 = fitFunc->GetParError(0);
    /*
    effFuncErr->SetParameters(0.0418139 , 0.0435551, 213236, 1.83373e-06,0.001672, 6651, 1.113e-08, 8657, 2.806e-08, 0.148);
    eff_d3->SetParameter(0, effFunc->GetParameter(3));
    */
    //effErr->SetParameters(fitFunc->GetParameter(0),fitFunc->GetParError(0));
    auto eff2 = fitFunc->GetParameter(0)/effFunc->Eval(980);
    auto effErr2 = abs(NU2/effFunc980U2)*sqrt(pow(effFuncErr980U2/effFunc980U2,2)+pow(NErrU2/NU2,2));
    cout << "N2_eff = " << eff2 << " +- " << effErr2 << endl;


    // NEXT DETECTOR

    string effFunc_U3 = "0.000168379 * (exp(-0.675181*log(x/1000) + 0.0356154*pow(log(x/1000), 2) - 539443/pow(x,3)))";
    string eff_d0_U3 = "log(x/1000) * (" + effFunc_U3 + ")";
    string eff_d1_U3 = "pow(log(x/1000), 2) * (" + effFunc_U3 + ")";
    string eff_d2_U3 = "-1/pow(x,3) * (" + effFunc_U3 + ")";
    string eff_d3_U3 = "1/0.000168379 * (" + effFunc_U3 + ")";
    string effFuncErr_U3 = "pow(1.81564e-06*(" + eff_d3_U3 + "),2) + pow(213215*(" + eff_d2_U3 + "),2) + pow(0.0430538*(" + eff_d1_U3 + "),2) + pow(0.0408049*(" + eff_d0_U3 + "),2) + 0.001602*(" + eff_d0_U3 + ")*(" + eff_d1_U3 + ") + 6390*(" + eff_d0_U3 + ")*(" + eff_d2_U3 + ") - 1.22e-08*(" + eff_d0_U3 + ")*(" + eff_d3_U3 + ") + 8548*(" + eff_d1_U3 + ")*(" + eff_d2_U3 + ") - 2.93e-08*(" + eff_d1_U3 + ")*(" + eff_d3_U3 + ") - 0.1538*(" + eff_d2_U3 + ")*(" + eff_d3_U3 + ")";
    string effErr_U3 = "sqrt(pow((1.04967e+00)/(" + effFunc_U3 + "),2)*(pow((" + effFuncErr_U3 + ")/(" + effFunc_U3 + "),2) + pow(1.04967e+00/(3.07743e-02),2)))";
    auto effErrU3 = new TF1("effErrU3", effErr_U3.c_str());
    cout << "Fitting C3:" << endl;
    hist3->Fit("fitFunc", "", "", xMin + BG_del, BG_open + BG_del);
    effFunc->SetParameters(-0.675181, 0.0356154, 539443, 0.000168379, 1000); //setting parameters given from Mads' calibration
    auto effFuncErrU3 = new TF1("effFuncErrU3", effFuncErr_U3.c_str());
    auto effFuncErr980U3 = effFuncErrU3->Eval(980);
    auto effFunc980U3 = effFunc->Eval(980);
    auto NU3 = fitFunc->GetParameter(0);
    auto NErrU3 = fitFunc->GetParError(0);
    /*
    effFuncErr->SetParameters(0.0408049, 0.0430538, 213215, 1.81564e-06,0.001602, 6390, 1.22e-08, 8548, 2.93e-08, 0.1538);
    eff_d3->SetParameter(0, effFunc->GetParameter(3));
    effErr->SetParameters(fitFunc->GetParameter(0),fitFunc->GetParError(0));
    */
    auto eff3 = fitFunc->GetParameter(0)/effFunc->Eval(980);
    auto effErr3 = abs(NU3/effFunc980U3)*sqrt(pow(effFuncErr980U3/effFunc980U3,2)+pow(NErrU3/NU3,2));
    cout << "N3_eff = " << eff3 << " +- " << effErr3 << endl;

    // NEXT DETECTOR

    string effFunc_U4 = "0.000174134 * (exp(-0.65969*log(x/1000) + 0.0459889*pow(log(x/1000), 2) - 595718/pow(x,3)))";
    string eff_d0_U4 = "log(x/1000) * (" + effFunc_U4 + ")";
    string eff_d1_U4 = "pow(log(x/1000), 2) * (" + effFunc_U4 + ")";
    string eff_d2_U4 = "-1/pow(x,3) * (" + effFunc_U4 + ")";
    string eff_d3_U4 = "1/0.000174134 * (" + effFunc_U4 + ")";
    string effFuncErr_U4 = "pow(1.82149e-06*(" + eff_d3_U4 + "),2) + pow(212226*(" + eff_d2_U4 + "),2) + pow(0.0427729*(" + eff_d1_U4 + "),2) + pow(0.0403567*(" + eff_d0_U4 + "),2) + 0.001577*(" + eff_d0_U4 + ")*(" + eff_d1_U4 + ") + 6326*(" + eff_d0_U4 + ")*(" + eff_d2_U4 + ") - 1.261e-08*(" + eff_d0_U4 + ")*(" + eff_d3_U4 + ") + 8460*(" + eff_d1_U4 + ")*(" + eff_d2_U4 + ") - 2.94e-08*(" + eff_d1_U4 + ")*(" + eff_d3_U4 + ") - 0.154*(" + eff_d2_U4 + ")*(" + eff_d3_U4 + ")";
    string effErr_U4 = "sqrt(pow((1.10658e+00)/(" + effFunc_U4 + "),2)*(pow((" + effFuncErr_U4 + ")/(" + effFunc_U4 + "),2) + pow(1.10658e+00/(3.07743e-02),2)))";
    auto effErrU4 = new TF1("effErrU4", effErr_U4.c_str());
    cout << "Fitting C4:" << endl;
    hist4->Fit("fitFunc", "", "", xMin + BG_del, BG_open + BG_del);
    effFunc->SetParameters(-0.65969, 0.0459889, 595718, 0.000174134, 1000); //setting parameters given from Mads' calibration
    auto effFuncErrU4 = new TF1("effFuncErrU4", effFuncErr_U4.c_str());
    auto effFuncErr980U4 = effFuncErrU4->Eval(980);
    auto effFunc980U4 = effFunc->Eval(980);
    auto NU4 = fitFunc->GetParameter(0);
    auto NErrU4 = fitFunc->GetParError(0);
    /*
    effFuncErr->SetParameters(0.0403567, 0.0427729, 212226, 1.82149e-06,0.001577, 6326,1.261e-08, 8460,2.94e-08, 0.154);
    eff_d3->SetParameter(0, effFunc->GetParameter(3));
    effErr->SetParameters(fitFunc->GetParameter(0),fitFunc->GetParError(0));
    */
    auto eff4 = fitFunc->GetParameter(0)/effFunc->Eval(980);
    auto effErr4 = abs(NU4/effFunc980U4)*sqrt(pow(effFuncErr980U4/effFunc980U4,2)+pow(NErrU4/NU4,2));

    cout << "N4_eff = " << eff4 << " +- " << effErr4 << endl;

    cout << "number of bins fitted = " << hist1->GetXaxis()->FindBin(BG_open+BG_del) - hist1->GetXaxis()->FindBin(xMin+BG_del) << endl;

    double effs[4] = {eff1,eff2,eff3,eff4};
    double effErrs[4] = {effErr1,effErr2,effErr3,effErr4};
    double clovers[4] = {1,2,3,4};

    auto effCanv = new TCanvas("", "", 1000, 800);

    //TGraph *effGraph = new TGraph(4,clovers, effs);
    TGraph *effGraph = new TGraphErrors(4,clovers, effs, 0, effErrs);
    effGraph->SetTitle("Efficiency corregated value of N");
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