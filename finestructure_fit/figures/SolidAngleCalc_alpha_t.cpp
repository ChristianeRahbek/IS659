//
// Created by Christiane Rahbek on 04-05-2023.
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

using namespace std;
using namespace AUSA::Sort;
using namespace AUSA;
using namespace libconfig;

class SolidAngleClac_alpha_t {
public:
    string filename;
    std::vector<std::vector<double>> U1_SAM, U2_SAM, U3_SAM, U4_SAM, U1_un, U2_un, U3_un, U4_un;
    std::vector<double> SAs, hitAngles;

    int matrixSize;

    shared_ptr<Setup> setup;

    shared_ptr<DoubleSidedDetector> det1, det2, det3, det4;

    SolidAngleClac_alpha_t(string filename, string setupFilename) {
        this->filename = filename;

        setup = JSON::readSetupFromJSON(setupFilename);

        det1 = setup->getDSSD(0);
        det2 = setup->getDSSD(1);
        det3 = setup->getDSSD(2);
        det4 = setup->getDSSD(3);

        matrixSize = 16;

        U1_SAM = readMatrixFromFile(filename, 9, matrixSize);
        U2_SAM = readMatrixFromFile(filename, 27, matrixSize);
        U3_SAM = readMatrixFromFile(filename, 45, matrixSize);
        U4_SAM = readMatrixFromFile(filename, 63, matrixSize);

    }

    void calculateSAs() {
        for (int i1 = 0; i1 < matrixSize; i1++) {
            for (int j1 = 0; j1 < matrixSize; j1++) { //now I can find i'th, j'th entry in first matrix
                for (int i2 = 0; i2 < matrixSize; i2++) {
                    for (int j2 = 0; j2 < matrixSize; j2++) { //now I can find i'th, j'th entry in second matrix
                        double saU11 = calculateSA(U1_SAM, U1_SAM, i1, j1, i2, j2);
                        double saU12 = calculateSA(U1_SAM, U2_SAM, i1, j1, i2, j2);
                        double saU13 = calculateSA(U1_SAM, U3_SAM, i1, j1, i2, j2);
                        double saU14 = calculateSA(U1_SAM, U4_SAM, i1, j1, i2, j2);
                        double saU22 = calculateSA(U2_SAM, U2_SAM, i1, j1, i2, j2);
                        double saU23 = calculateSA(U2_SAM, U3_SAM, i1, j1, i2, j2);
                        double saU24 = calculateSA(U2_SAM, U4_SAM, i1, j1, i2, j2);
                        double saU33 = calculateSA(U3_SAM, U3_SAM, i1, j1, i2, j2);
                        double saU34 = calculateSA(U3_SAM, U4_SAM, i1, j1, i2, j2);
                        double saU44 = calculateSA(U4_SAM, U4_SAM, i1, j1, i2, j2);

                        auto p1U1 = det1->getPixelPosition(i1, j1);
                        auto p1U2 = det2->getPixelPosition(i1, j1);
                        auto p1U3 = det3->getPixelPosition(i1, j1);
                        auto p1U4 = det4->getPixelPosition(i1, j1);
                        auto p2U1 = det1->getPixelPosition(i2, j2);
                        auto p2U2 = det2->getPixelPosition(i2, j2);
                        auto p2U3 = det2->getPixelPosition(i2, j2);
                        auto p2U4 = det4->getPixelPosition(i2, j2);

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

    double calculateSA(std::vector<std::vector<double>> m1, std::vector<std::vector<double>> m2, int i1, int j1, int i2,
                       int j2) {
        double omega1 = m1[i1][j1];
        double omega2 = m2[i2][j2];
        if (i1 != i2 && j1 != j2) {
            return TMath::Power(omega1, 2) + TMath::Power(omega2, 2) + 2 * omega1 * omega2;
        } else return TMath::Power(omega1, 2) + TMath::Power(omega2, 2);
    }

    void plotSA() {

        auto *canv = new TCanvas("", "", 1000, 800);

        int size = SAs.size();

        auto *graph = new TGraph(size);

        for (int i = 0; i < size; i++) {
            graph->SetPoint(i, hitAngles[i], SAs[i]);
        }

        cout << graph->GetN() << endl;

        graph->Draw();

        graph->GetXaxis()->SetTitle("Hit angle");
        graph->GetYaxis()->SetTitle("Solid angle pr bin");
        graph->SetTitle("Solid angle");

        canv->SaveAs("SolidAngleTritonAlpha.pdf");
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
};

int main() {
    string filename = EUtil::getProjectRoot() + "/../SAM8He.txt";
    string setupFilename = EUtil::getProjectRoot() + "/../../setup/setup.json";
    string errFilename = EUtil::getProjectRoot() + "/../SAM8He.txt_mUS";
    string errSetupFilename = EUtil::getProjectRoot() + "/../../setup/setup_mUS.json";

    auto SAcalc = new SolidAngleClac_alpha_t(filename, setupFilename);
    auto SAcalcErr = new SolidAngleClac_alpha_t(errFilename, errSetupFilename);

    SAcalc->calculateSAs();
    SAcalc->plotSA();

    return EXIT_SUCCESS;
}