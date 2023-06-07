//
// Created by Christiane Rahbek on 04-05-2023.
// Calculates solid angle detector efficiency and its uncertainty for non-uniform decays.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <TMath.h>
#include <ausa/sort/analyzer/AbstractSortedAnalyzer.h>
#include <ausa/json/IO.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <ausa/util/FileUtil.h>
#include <libconfig.h++>
#include "projectutil.h"
#include <TH2D.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>

using namespace std;
using namespace AUSA::Sort;
using namespace AUSA;
using namespace libconfig;

class SolidAngleClac_alpha_t {
public:
    string filename;
    std::vector<std::vector<double>> U1_SAM, U2_SAM, U3_SAM, U4_SAM, U1_un, U2_un, U3_un, U4_un;
    std::vector<double> SAs, SAerrs, hitAngles;
    double SA_found, SA_found_err;

    int matrixSize;

    TH1F* effHist;

    std::vector<double> hitAngBins, solidAngsBins, solidAngsBinsAlt;

    shared_ptr<Setup> setup, errSetup;

    shared_ptr<DoubleSidedDetector> det1, det2, det3, det4;

    double saU11, saU12, saU13, saU14, saU22, saU23, saU24, saU33, saU34, saU44;
    double saU11err, saU12err, saU13err, saU14err, saU22err, saU23err, saU24err, saU33err, saU34err, saU44err;


    SolidAngleClac_alpha_t(string filename, string setupFilename, string errFilename) {
        this->filename = filename;

        setup = JSON::readSetupFromJSON(setupFilename);

        effHist = new TH1F("effHist", "effHist", 37, 0, 180);

        det1 = setup->getDSSD(0);
        det2 = setup->getDSSD(1);
        det3 = setup->getDSSD(2);
        det4 = setup->getDSSD(3);

        matrixSize = 16;

        U1_SAM = readMatrixFromFile(filename, 9, matrixSize);
        U2_SAM = readMatrixFromFile(filename, 27, matrixSize);
        U3_SAM = readMatrixFromFile(filename, 45, matrixSize);
        U4_SAM = readMatrixFromFile(filename, 63, matrixSize);
        U1_un = readMatrixFromFile(errFilename, 9, matrixSize);
        U2_un = readMatrixFromFile(errFilename, 27, matrixSize);
        U3_un = readMatrixFromFile(errFilename, 45, matrixSize);
        U4_un = readMatrixFromFile(errFilename, 63, matrixSize);

        saU11 = 0; saU12 = 0; saU13 = 0; saU14 = 0; saU22 = 0; saU23 = 0; saU24 = 0; saU33 = 0; saU34 = 0; saU44 = 0;
        saU11err = 0; saU12err = 0; saU13err = 0; saU14err = 0; saU22err = 0; saU23err = 0; saU24err = 0; saU33err = 0;
        saU34err = 0; saU44err = 0;

        int a = 5;
        for(int n = 0; n <=180/a; n++) {
            hitAngBins.push_back(a*n);
        }
    }

    void calculateSAs() {
        for (int i1 = 0; i1 < matrixSize; i1++) {
            for (int j1 = 0; j1 < matrixSize; j1++) { //now I can find i'th, j'th entry in first matrix
                for (int i2 = 0; i2 < matrixSize; i2++) {
                    for (int j2 = 0; j2 < matrixSize; j2++) { //now I can find i'th, j'th entry in second matrix

                        saU11 = calculateSA(U1_SAM, U1_SAM, i1, j1, i2, j2);
                        saU12 = calculateSA(U1_SAM, U2_SAM, i1, j1, i2, j2);
                        saU13 = calculateSA(U1_SAM, U3_SAM, i1, j1, i2, j2);
                        saU14 = calculateSA(U1_SAM, U4_SAM, i1, j1, i2, j2);
                        saU22 = calculateSA(U2_SAM, U2_SAM, i1, j1, i2, j2);
                        saU23 = calculateSA(U2_SAM, U3_SAM, i1, j1, i2, j2);
                        saU24 = calculateSA(U2_SAM, U4_SAM, i1, j1, i2, j2);
                        saU33 = calculateSA(U3_SAM, U3_SAM, i1, j1, i2, j2);
                        saU34 = calculateSA(U3_SAM, U4_SAM, i1, j1, i2, j2);
                        saU44 = calculateSA(U4_SAM, U4_SAM, i1, j1, i2, j2);

                        saU11err = calculateSAerrSq(U1_SAM, U1_SAM,U1_un, U1_un, i1, j1, i2, j2);
                        saU12err = calculateSAerrSq(U1_SAM, U2_SAM,U1_un, U2_un, i1, j1, i2, j2);
                        saU13err = calculateSAerrSq(U1_SAM, U3_SAM,U1_un, U3_un, i1, j1, i2, j2);
                        saU14err = calculateSAerrSq(U1_SAM, U4_SAM,U1_un, U4_un, i1, j1, i2, j2);
                        saU22err = calculateSAerrSq(U2_SAM, U2_SAM,U2_un, U2_un, i1, j1, i2, j2);
                        saU23err = calculateSAerrSq(U2_SAM, U3_SAM,U2_un, U3_un, i1, j1, i2, j2);
                        saU24err = calculateSAerrSq(U2_SAM, U4_SAM,U2_un, U4_un, i1, j1, i2, j2);
                        saU33err = calculateSAerrSq(U3_SAM, U3_SAM,U3_un, U3_un, i1, j1, i2, j2);
                        saU34err = calculateSAerrSq(U3_SAM, U4_SAM,U3_un, U4_un, i1, j1, i2, j2);
                        saU44err = calculateSAerrSq(U4_SAM, U4_SAM,U4_un, U4_un, i1, j1, i2, j2);

                        auto& p1U1 = det1->getPixelPosition(i1+1, j1+1);
                        auto& p1U2 = det2->getPixelPosition(i1+1, j1+1);
                        auto& p1U3 = det3->getPixelPosition(i1+1, j1+1);
                        auto& p1U4 = det4->getPixelPosition(i1+1, j1+1);
                        auto& p2U1 = det1->getPixelPosition(i2+1, j2+1);
                        auto& p2U2 = det2->getPixelPosition(i2+1, j2+1);
                        auto& p2U3 = det3->getPixelPosition(i2+1, j2+1);
                        auto& p2U4 = det4->getPixelPosition(i2+1, j2+1);

                        SAs.push_back(saU11);
                        SAs.push_back(saU12);
                        SAs.push_back(saU13);
                        SAs.push_back(saU14);
                        SAs.push_back(saU22);
                        SAs.push_back(saU23);
                        SAs.push_back(saU24);
                        SAs.push_back(saU33);
                        SAs.push_back(saU34);
                        SAs.push_back(saU44);

                        SAerrs.push_back(saU11err);
                        SAerrs.push_back(saU12err);
                        SAerrs.push_back(saU13err);
                        SAerrs.push_back(saU14err);
                        SAerrs.push_back(saU22err);
                        SAerrs.push_back(saU23err);
                        SAerrs.push_back(saU24err);
                        SAerrs.push_back(saU33err);
                        SAerrs.push_back(saU34err);
                        SAerrs.push_back(saU44err);

                        hitAngles.push_back(p1U1.Angle(p2U1) * TMath::RadToDeg());
                        hitAngles.push_back(p1U1.Angle(p2U2) * TMath::RadToDeg());
                        hitAngles.push_back(p1U1.Angle(p2U3) * TMath::RadToDeg());
                        hitAngles.push_back(p1U1.Angle(p2U4) * TMath::RadToDeg());
                        hitAngles.push_back(p1U2.Angle(p2U2) * TMath::RadToDeg());
                        hitAngles.push_back(p1U2.Angle(p2U3) * TMath::RadToDeg());
                        hitAngles.push_back(p1U2.Angle(p2U4) * TMath::RadToDeg());
                        hitAngles.push_back(p1U3.Angle(p2U3) * TMath::RadToDeg());
                        hitAngles.push_back(p1U3.Angle(p2U4) * TMath::RadToDeg());
                        hitAngles.push_back(p1U4.Angle(p2U4) * TMath::RadToDeg());


                    }
                }
            }
        }
    }

    void calcSABins() {
        for(int i = 0; i<hitAngBins.size(); i++) {
            double solidAngs = 0;
            double solidAngErr = 0;
            for(int j = 0; j < hitAngles.size(); j++) {
                auto ang = hitAngles[j];
                if((ang >= hitAngBins[i]) && (ang < hitAngBins[i+1])) {
                    solidAngs += SAs[j];
                    solidAngErr += SAerrs[j]; //FIX UNCERTAINTIES, THIS IS NOT RIGHT
                }
            }
            solidAngsBins.push_back(solidAngs);
            solidAngsBinsAlt.push_back(solidAngs/TMath::Power(4*TMath::Pi(),2));
        }
    }

    double calcEffectiveNoOfDetections(string dataFile, string selectionCrit, int degs){
        //return the effective number of detections in file "dataFile" with selection criteria "selectionCrit"

        cout << "hitAngBins.size() = " << hitAngBins.size() << endl;
        cout << "solidAngsBins.size() = " << solidAngsBins.size() << endl;

        auto canv = new TCanvas();
        double totalEff = 0;
        auto *f= new TFile(dataFile.c_str());
        auto *tr=(TTree*)f->Get("a");
        auto hist = new TH1F("hist", "hist", hitAngBins.size(), 0, 180);
        tr->Draw("hitAng >> hist", selectionCrit.c_str());
        canv->Update();
        canv->Draw();
        canv->SaveAs((EUtil::getProjectRoot()+"/DataHitAng" + to_string(degs) + "degs.pdf").c_str());

        auto canv1 = new TCanvas();
        auto effNHist = new TH1F("effNHist", "effNHist", hitAngBins.size(), 0, 180);

        for(int i = 0; i < hist->GetNbinsX(); i++) {
            auto entry = hist->GetBinContent(i+1);
            auto eff = solidAngsBins[i];
            double effEntry;
            if(isnan(eff) || isnan(entry)) {
                cout << "i = " << i << endl;
                cout << "bincontent = " << entry << endl;
                cout << "eff = " << eff << endl;
            }
            if(eff==0) effEntry = 0;
            else effEntry = (entry/eff);//TMath::Power(4*TMath::Pi(),2);

            totalEff += effEntry;

            effNHist->SetBinContent(i+1, effEntry);
            //effHist->SetBinContent(i,effEntry); //tried to return this, but I changed my mind...
        }

        effNHist->Draw();
        canv1->Update();
        canv1->Draw();
        canv1->SaveAs((EUtil::getProjectRoot()+"/EffNHist" + to_string(degs) + "degCoincidences.pdf").c_str());

        cout << "total eff = " << totalEff << endl;
        cout << "" << endl;
        return totalEff;
    }

    void plotSA() {

        auto *canv = new TCanvas("", "", 1500, 700);

        canv->cd();

        int size = SAs.size();

        auto *graph = new TGraph(size, &hitAngles[0], &SAs[0]);


        graph->Draw("");

        graph->DrawGraph(size,&hitAngles[0], &SAs[0], "B");

        graph->GetXaxis()->SetTitle("Hit angle");
        graph->GetYaxis()->SetTitle("Solid angle pr bin");
        graph->SetTitle("Solid angle");

        canv->Update();

        canv->SaveAs((EUtil::getProjectRoot() + "/SolidAngleTritonAlpha.png").c_str());


    }

    void plotEff(bool isAlt) {
        auto canv = new TCanvas();
        TGraph *graph;
        if(isAlt) {
            graph = new TGraph(solidAngsBins.size(), &hitAngBins[0], &solidAngsBinsAlt[0]);
        }
        else {
            graph = new TGraph(solidAngsBins.size(), &hitAngBins[0], &solidAngsBins[0]);
        }
        graph->Draw();
        canv->Update();
        canv->Draw();
        if(isAlt) {
            canv->SaveAs((EUtil::getProjectRoot()+"/EffPlotAlt.pdf").c_str());
        } else {
            canv->SaveAs((EUtil::getProjectRoot()+"/EffPlot.pdf").c_str());
        }
    }

    std::vector<double> findingTotalSAinRange(double theta_min, double theta_max){
        std::vector<double> returnVec;
        SA_found = 0;
        SA_found_err = 0;
        for(int i = 0; i < hitAngles.size(); i++){
            double hitAng = hitAngles[i];
            double SA = SAs[i];
            double SAerr = SAerrs[i];
            if(hitAng>=theta_min && hitAng<=theta_max) {
                SA_found += SA;
                SA_found_err += SAerr;
            }
        }
        returnVec.push_back(SA_found);
        returnVec.push_back(sqrt(SA_found_err));
        return returnVec;
    }

    void printDataToTxt() {
        string filename = EUtil::getProjectRoot() + "/solidAnglePairData.txt";
        ofstream outfile(filename);
        if (!outfile.is_open()) {
            cerr << "Error: Unable to open file \n";
        }
        outfile << "hitAngles" << " " << "Solid angles" << endl;
        for(int i = 0; i < SAs.size(); i++) {
            if(SAs[i]>0){
                outfile << hitAngles[i] << " " << SAs[i] << endl;
            }
        }
        outfile.close();
    }

private:
    std::vector<std::vector<double>> readMatrixFromFile(string file, int from, int matrixSize) {
        std::ifstream input(file);
        std::vector<std::vector<double>> Matrix;
        std::string line;
        int lineNumber = 0;

        while (std::getline(input, line)) {
            lineNumber++;
            if (lineNumber <= from) {
                continue;
            }
            if (lineNumber > from + matrixSize) {
                break;
            }
            std::istringstream iss(line);
            std::vector<double> row;
            double value;
            for (int i = 0; i < matrixSize; i++) {
                iss >> value;
                row.push_back(value);
            }
            Matrix.push_back(row);
        }
        return Matrix;
    }

    double calculateSA(std::vector<std::vector<double>> m1, std::vector<std::vector<double>> m2, int i1, int j1, int i2,
                       int j2) {
        double omega1 = m1[i1][j1]/(4*TMath::Pi());
        double omega2 = m2[i2][j2]/(4*TMath::Pi());

        if(i1 == i2 && j1 == j2 && m1 == m2) return 0;
        else return omega1 * omega2;
        /*
        if (i1 != i2 && j1 != j2 && m1 != m2) {
            return omega1 * omega2;
        } else return 0;
         */
    }

    double calculateSAerrSq(std::vector<std::vector<double>> m1, std::vector<std::vector<double>> m2,
                            std::vector<std::vector<double>> m1err, std::vector<std::vector<double>> m2err,
                            int i1, int j1, int i2, int j2) {
        double omega1 = m1[i1][j1]/(4*TMath::Pi());
        double omega2 = m2[i2][j2]/(4*TMath::Pi());
        double omega1err = TMath::Abs(m1err[i1][j1]/(4*TMath::Pi())-omega1);
        double omega2err = TMath::Abs(m2err[i1][j1]/(4*TMath::Pi())-omega2);

        if(omega1==0 || omega2==0) return 0;

        if (i1 != i2 && j1 != j2 && m1 != m2) {
            double f = 2 * omega1 * omega2;
            return TMath::Power(f,2) * (TMath::Power(omega1err/omega1,2)
                   + TMath::Power(omega2err/omega2,2));
        } else {
            double f1 = pow(omega1, 2);
            double f2 = pow(omega2, 2);
            return TMath::Power(f1,2) * (TMath::Power(2*omega1err/omega1,2))
                   + TMath::Power(f2,2) * (TMath::Power(2*omega2err/omega2,2));
        }
    }
};

int main() {
    string filename = EUtil::getProjectRoot() + "/../SAM8He.txt";
    string setupFilename = EUtil::getProjectRoot() + "/../../setup/setup.json";
    string errFilename = EUtil::getProjectRoot() + "/../SAM8He_mUS.txt";
    string dataFileName = EUtil::getProjectRoot() + "/../analysis1/output/doubleDSSD/mergedRuns.root";

    cout << "defining SAcalc" << endl;
    auto SAcalc = new SolidAngleClac_alpha_t(filename, setupFilename, errFilename);

    cout << "calculating solid angles" << endl;
    SAcalc->calculateSAs();
    cout << "plotting" << endl;
    SAcalc->plotSA();

    cout << "saving to txt" << endl;
    SAcalc->printDataToTxt();
    SAcalc->calcSABins();

    double totalSA = 0;
    for(int i = 0; i < SAcalc->SAs.size(); i++) {
        totalSA += (SAcalc->SAs)[i];
    }

    cout << "total efficiency of all detectors " << totalSA << "radians/(4pi)^2" << endl;

    string outfileName = EUtil::getProjectRoot() + "/efficiencyOfCounts.txt";
    ofstream outfile(outfileName);
    if (!outfile.is_open()) {
        cerr << "Error: Unable to open file \n";
    } else {
        outfile << "Alpha-t coincidences" << endl;
        auto effReactions0Deg = SAcalc->calcEffectiveNoOfDetections(
                dataFileName, "id0==id1 && abs(FT0-FT1)<1500 && (Edep0<1000 || Edep1<1000)", 0);
        outfile << "Total number of decays in 0deg range is " << effReactions0Deg << endl;

        auto effReactions180Deg = SAcalc->calcEffectiveNoOfDetections(
                dataFileName,
                "((id0==0 && id1==2) || (id0==1 && id1==3) || (id0==2 && id1==0) || (id0==3 && id1==1)) && abs(FT0-FT1)<1500 && hitAng < 130", 180);
        outfile << "Total number of decays in 180deg range is " << effReactions180Deg << endl;

        auto eff8Li = SAcalc->calcEffectiveNoOfDetections(
                dataFileName,
                "((id0==0 && id1==2) || (id0==1 && id1==3) || (id0==2 && id1==0) || (id0==3 && id1==1)) && abs(FT0-FT1)<1500 && hitAng > 130 && abs(Edep0-Edep1)<500",
                130);
        outfile << "Total number of decays in 8Li is " << eff8Li << endl;

        SAcalc->plotEff(true);
        SAcalc->plotEff(false);
        outfile.close();
    }

    return EXIT_SUCCESS;
}