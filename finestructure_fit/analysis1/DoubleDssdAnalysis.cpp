#include <iostream>
#include <string>
#include <ausa/json/IO.h>
#include <ausa/sort/analyzer/AbstractSortedAnalyzer.h>
#include <ausa/sort/SortedReader.h>
#include <ausa/util/StringUtil.h>
#include <ausa/util/FileUtil.h>
#include <ausa/setup/DoubleSidedSiliconDetector.h>
#include <ausa/util/DynamicBranchVector.h>
#include <ausa/eloss/Ion.h>
#include <ausa/eloss/Default.h>
#include <ausa/constants/Mass.h>
#include <ausa/output/OutputConvenience.h>
#include "Hit.h"
#include <Math/Vector3D.h>
#include <TROOT.h>
#include <ctime>
#include <libconfig.h++>
#include "projectutil.h"
#include "IS659Detector.h"


using namespace std;
using namespace AUSA;
using namespace ROOT::Math;
using namespace TMath;
using namespace AUSA::Sort;
using namespace AUSA::EnergyLoss;
using namespace AUSA::Constants;
using namespace libconfig;

#define NAN_TVECTOR3 TVector3(NAN, NAN, NAN)


class SingleDssdAnalysis {
public:
    explicit SingleDssdAnalysis(TString outfilename) {
        NUM = 0;

        //this->filePath = filePath;

        isBetas = 0;
        isParticleBeta = 0;

        output = new TFile(outfilename, "RECREATE");
        t = new TTree("a", "a");
        t->Branch("mul", &mul);
        t->Branch("num", &NUM);
        t->Branch("dssdMul", &dssdMul);
        t->Branch("plasticMul", &plasticMul);
        t->Branch("TPROTONS", &TPROTONS);

        v_dir0 = make_unique<DynamicBranchVector<TVector3>>(*t, "dir0");
        v_pos0 = make_unique<DynamicBranchVector<TVector3>>(*t, "pos0");
        v_dir1 = make_unique<DynamicBranchVector<TVector3>>(*t, "dir1");
        v_pos1 = make_unique<DynamicBranchVector<TVector3>>(*t, "pos1");

        v_theta0 = make_unique<DynamicBranchVector<double>>(*t, "theta0", "mul");
        v_ang0 = make_unique<DynamicBranchVector<double>>(*t, "angle0", "mul");
        v_theta1 = make_unique<DynamicBranchVector<double>>(*t, "theta1", "mul");
        v_ang1 = make_unique<DynamicBranchVector<double>>(*t, "angle1", "mul");

        v_Edep0 = make_unique<DynamicBranchVector<double>>(*t, "Edep0", "mul");
        v_Ea0 = make_unique<DynamicBranchVector<double>>(*t, "Ea0", "mul");
        v_BE0 = make_unique<DynamicBranchVector<double>>(*t, "BE0", "mul");
        v_FE0 = make_unique<DynamicBranchVector<double>>(*t, "FE0", "mul");
        v_Et0 = make_unique<DynamicBranchVector<double>>(*t, "Et0", "mul");
        v_Edep1 = make_unique<DynamicBranchVector<double>>(*t, "Edep1", "mul");
        v_Ea1 = make_unique<DynamicBranchVector<double>>(*t, "Ea1", "mul");
        v_BE1 = make_unique<DynamicBranchVector<double>>(*t, "BE1", "mul");
        v_FE1 = make_unique<DynamicBranchVector<double>>(*t, "FE1", "mul");
        v_Et1 = make_unique<DynamicBranchVector<double>>(*t, "Et1", "mul");

        v_FT0 = make_unique<DynamicBranchVector<double>>(*t, "FT0", "mul");
        v_BT0 = make_unique<DynamicBranchVector<double>>(*t, "BT0", "mul");
        v_FT1 = make_unique<DynamicBranchVector<double>>(*t, "FT1", "mul");
        v_BT1 = make_unique<DynamicBranchVector<double>>(*t, "BT1", "mul");

        v_dE0 = make_unique<DynamicBranchVector<double>>(*t, "dE0", "mul");
        v_Ecm0 = make_unique<DynamicBranchVector<double>>(*t, "Ecm0", "mul");
        v_dE1 = make_unique<DynamicBranchVector<double>>(*t, "dE1", "mul");
        v_Ecm1 = make_unique<DynamicBranchVector<double>>(*t, "Ecm1", "mul");

        v_i0 = make_unique<DynamicBranchVector<short>>(*t, "id0", "mul");
        v_i1 = make_unique<DynamicBranchVector<short>>(*t, "id1", "mul");

        v_isBetas = make_unique<DynamicBranchVector<int>>(*t, "isBetas", "mul");
        v_isParticleBeta = make_unique<DynamicBranchVector<int>>(*t, "isParticleBeta", "mul");

        v_F0 = make_unique<DynamicBranchVector<short>>(*t, "FI0", "mul");
        v_B0 = make_unique<DynamicBranchVector<short>>(*t, "BI0", "mul");
        v_F1 = make_unique<DynamicBranchVector<short>>(*t, "FI1", "mul");
        v_B1 = make_unique<DynamicBranchVector<short>>(*t, "BI1", "mul");

        v_hitAng = make_unique<DynamicBranchVector<double>>(*t, "hitAng", "mul");
    }

    void findSingleDssdHits(string filePath) {
        cout << "making tree" << endl;
        auto *f = new TFile(filePath.c_str());
        auto *tr = (TTree *) f->Get("a");

        UInt_t maxHits = 20;
        Int_t numOrig;
        UInt_t mulOrig, TPROTONSOrig, plasticMulOrig, dssdMulOrig;
        Short_t idOrig[maxHits], FIOrig[maxHits], BIOrig[maxHits];
        Double_t EdepOrig[maxHits], EaOrig[maxHits], EtOrig[maxHits], BEOrig[maxHits], FEOrig[maxHits],
                thetaOrig[maxHits], angOrig[maxHits], FTOrig[maxHits], BTOrig[maxHits], dEOrig[maxHits],
                EcmOrig[maxHits], posX[maxHits], posY[maxHits], posZ[maxHits], dirX[maxHits],
                dirY[maxHits], dirZ[maxHits];

        cout << "setting branches in tree" << flush << endl;
        tr->SetBranchAddress("mul", &mulOrig);
        tr->SetBranchAddress("TPROTONS", &TPROTONSOrig);
        tr->SetBranchAddress("num", &numOrig);
        tr->SetBranchAddress("plasticMul", &plasticMulOrig);
        tr->SetBranchAddress("dssdMul", &dssdMulOrig);
        tr->SetBranchAddress("id", idOrig);
        tr->SetBranchAddress("FI", FIOrig);
        tr->SetBranchAddress("BI", BIOrig);
        tr->SetBranchAddress("Edep", EdepOrig);
        tr->SetBranchAddress("Ea", EaOrig);
        tr->SetBranchAddress("Et", EtOrig);
        tr->SetBranchAddress("BE", BEOrig);
        tr->SetBranchAddress("FE", FEOrig);
        tr->SetBranchAddress("theta", thetaOrig);
        tr->SetBranchAddress("angle", angOrig);
        tr->SetBranchAddress("FT", FTOrig);
        tr->SetBranchAddress("BT", BTOrig);
        tr->SetBranchAddress("dE", dEOrig);
        tr->SetBranchAddress("Ecm", EcmOrig);
        tr->SetBranchAddress("posX", posX);
        tr->SetBranchAddress("posY", posY);
        tr->SetBranchAddress("posZ", posZ);
        tr->SetBranchAddress("dirX", dirX);
        tr->SetBranchAddress("dirY", dirY);
        tr->SetBranchAddress("dirZ", dirZ);


        output->cd();

        cout << "looping through entries" << flush << endl;
        for (UInt_t i = 0; i < tr->GetEntries(); i++) { //looping through each entry
            clear();
            tr->GetEntry(i);

            if (dssdMulOrig != 2) continue; //this is a double DSSD analysis

            for (UInt_t j = 0; j < mulOrig; j++) { //looping through each hit
                if (idOrig[j] > 3) continue; //we only want to save DSSD data

                for (UInt_t k = j + 1; k < mulOrig; k++) {//we want to find pairs of hits
                    if (idOrig[k] > 3) continue; //we only want to save DSSD data

                    if (plasticMulOrig > 0) {

                        isBetas = 0;
                        isParticleBeta = 0;

                        for (UInt_t l = 0; l < mulOrig; l++) {
                            if (idOrig[l] < 4) continue; //now searching for plastics

                            //is the plastic hit in one of the DSSDs that has been hit
                            if (!(short(idOrig[j] + 4) == idOrig[l] || short(idOrig[k] + 4) == idOrig[l])) continue;

                            double timeCutMin = 12000;
                            double timeCutMax = 21000;
                            bool timePB = (short(idOrig[j] + 4) == idOrig[l] && abs(FTOrig[j]-FTOrig[l]) < timeCutMax
                                                                             && abs(FTOrig[j]-FTOrig[l]) > timeCutMin)
                                    || (short(idOrig[k] + 4) == idOrig[l] && abs(FTOrig[j]-FTOrig[l]) < timeCutMax
                                                                             && abs(FTOrig[j]-FTOrig[l]) > timeCutMin);

                            //look for beta-particle here!
                            if (timePB) {
                                isParticleBeta = 1;
                                isBetas = 0;
                            }

                            for (UInt_t m = l + 1; m < mulOrig; m++) { //searching for beta beta
                                if (idOrig[l] < 4) continue; //searching for plastics

                                //is the plastic hit in one of the DSSDs that has been hit
                                bool hitBool = short(idOrig[j] + 4) == idOrig[m] || short(idOrig[k] + 4) == idOrig[m];
                                //beta beta time criteria
                                auto timeCutBBMin = timeCutMin;
                                auto timeCutBBMax = timeCutMax;
                                bool timeBB = (short(idOrig[j] + 4) == idOrig[l] && short(idOrig[k] + 4) == idOrig[m]
                                        && abs(FTOrig[j]-FTOrig[l]) < timeCutBBMax && abs(FTOrig[k]-FTOrig[m]) < timeCutBBMax
                                        && abs(FTOrig[j]-FTOrig[l]) > timeCutBBMin && abs(FTOrig[k]-FTOrig[m]) > timeCutBBMin)
                                        || (short(idOrig[k] + 4) == idOrig[l] && short(idOrig[j] + 4) == idOrig[m]
                                        && abs(FTOrig[k]-FTOrig[l]) < timeCutBBMax && abs(FTOrig[j]-FTOrig[m]) < timeCutBBMax
                                        && abs(FTOrig[k]-FTOrig[l]) > timeCutBBMin && abs(FTOrig[j]-FTOrig[m]) > timeCutBBMin);

                                //if we hit 2 different dssds
                                if (hitBool && idOrig[m] != idOrig[l] && timeBB) {
                                    isBetas = 1;
                                    isParticleBeta = 0;
                                    break; //we do not want to search more
                                }
                                //if we  hit the same dssd
                                if (hitBool && idOrig[m] == idOrig[l] && idOrig[j] == idOrig[k] && timeBB) {
                                    isBetas = 1;
                                    isParticleBeta = 0;
                                    break; //we do not want to search more
                                }
                            }
                        }
                    } else {
                        isBetas = 0; //false
                        isParticleBeta = 0; //false
                    }

                    v_theta0->add(thetaOrig[j]);
                    v_ang0->add(angOrig[j]);
                    v_Et0->add(EtOrig[j]);
                    v_BE0->add(BEOrig[j]);
                    v_FE0->add(FEOrig[j]);
                    v_FT0->add(FTOrig[j]);
                    v_BT0->add(BTOrig[j]);
                    v_dE0->add(dEOrig[j]);
                    v_Ecm0->add(EcmOrig[j]);
                    v_Ea0->add(EaOrig[j]);
                    v_Edep0->add(EdepOrig[j]);
                    v_i0->add(idOrig[j]);
                    v_F0->add(FIOrig[j]);
                    v_B0->add(BIOrig[j]);

                    v_theta1->add(thetaOrig[k]);
                    v_ang1->add(angOrig[k]);
                    v_Et1->add(EtOrig[k]);
                    v_BE1->add(BEOrig[k]);
                    v_FE1->add(FEOrig[k]);
                    v_FT1->add(FTOrig[k]);
                    v_BT1->add(BTOrig[k]);
                    v_dE1->add(dEOrig[k]);
                    v_Ecm1->add(EcmOrig[k]);
                    v_Ea1->add(EaOrig[k]);
                    v_Edep1->add(EdepOrig[k]);
                    v_i1->add(idOrig[k]);
                    v_F1->add(FIOrig[k]);
                    v_B1->add(BIOrig[k]);

                    auto pos0 = *new TVector3(posX[j], posY[j], posZ[j]);
                    auto pos1 = *new TVector3(posX[k], posY[k], posZ[k]);
                    auto dir0 = *new TVector3(dirX[j], dirY[j], dirZ[j]);
                    auto dir1 = *new TVector3(dirX[k], dirY[k], dirZ[k]);

                    v_pos0->add(pos0);
                    v_pos1->add(pos1);
                    v_dir0->add(dir0);
                    v_dir1->add(dir1);

                    v_hitAng->add(pos0.Angle(pos1) * TMath::RadToDeg());

                    v_isBetas->add(isBetas);
                    v_isParticleBeta->add(isParticleBeta);

                    TPROTONS = TPROTONSOrig;
                    plasticMul = plasticMulOrig;
                    dssdMul = dssdMulOrig;
                    NUM = numOrig;

                    mul++;
                }
            }

            //cout << "filling tree" << endl;
            t->Fill();
        }

        //cout << "returning tree" << flush << endl;
        t->Write();

        output->Close();
    }

    void clear() {
        mul = 0;
        dssdMul = 0;
        plasticMul = 0;
        AUSA::clear(
                *v_Et0, *v_Ea0, *v_theta0, *v_Edep0, *v_i0, *v_FE0, *v_BE0,
                *v_F0, *v_B0, *v_Ecm0, *v_ang0, *v_pos0, *v_dir0, *v_ang1,
                *v_dE0, *v_FT0, *v_BT0, *v_isBetas, *v_isParticleBeta, *v_hitAng,
                *v_Et1, *v_Ea1, *v_theta1, *v_Edep1, *v_i1, *v_FE1, *v_BE1,
                *v_F1, *v_B1, *v_Ecm1, *v_pos1, *v_dir1, *v_dE1, *v_FT1, *v_B1
        );
    }

    const map<unsigned short, IS659Detector *> detector_map = {
            {0, new IS659Detector(0, "U1", IS659Type::SquareDSSD)},
            {1, new IS659Detector(1, "U2", IS659Type::SquareDSSD)},
            {2, new IS659Detector(2, "U3", IS659Type::SquareDSSD)},
            {3, new IS659Detector(3, "U4", IS659Type::SquareDSSD)},
            {4, new IS659Detector(4, "P1", IS659Type::Plastic)},
            {5, new IS659Detector(5, "P2", IS659Type::Plastic)},
            {6, new IS659Detector(6, "P3", IS659Type::Plastic)},
            {7, new IS659Detector(7, "P4", IS659Type::Plastic)}
    };

    int NUM;
    TTree *t;
    unique_ptr<DynamicBranchVector<TVector3>> v_dir0, v_pos0, v_dir1, v_pos1;
    unique_ptr<DynamicBranchVector<double>> v_Edep0, v_Ea0, v_Et0, v_BE0, v_FE0, v_theta0, v_dE0, v_Ecm0;
    unique_ptr<DynamicBranchVector<double>> v_Edep1, v_Ea1, v_Et1, v_BE1, v_FE1, v_theta1, v_dE1, v_Ecm1;
    unique_ptr<DynamicBranchVector<short>> v_i0, v_i1;
    unique_ptr<DynamicBranchVector<int>> v_isBetas, v_isParticleBeta;
    unique_ptr<DynamicBranchVector<short>> v_F0, v_B0, v_F1, v_B1;
    unique_ptr<DynamicBranchVector<double>> v_ang0, v_ang1, v_hitAng;
    unique_ptr<DynamicBranchVector<double>> v_FT0, v_BT0, v_FT1, v_BT1;


    TFile *output;
    UInt_t mul{}, TPROTONS{}, dssdMul{}, plasticMul{};
//string filePath;
    int isBetas, isParticleBeta;
};

string setup_path, target_path, input_path, output_dir;

void prepareFileIO(const string &configfile) {
    Config cfg;
    cfg.readFile(configfile.c_str());

    if (cfg.lookup("paths_relative_to_project_root")) {
        setup_path = EUtil::getProjectRoot() + "/" + cfg.lookup("setup_file").c_str();
        target_path = EUtil::getProjectRoot() + "/" + cfg.lookup("target_file").c_str();
        input_path = EUtil::getProjectRoot() + "/" + cfg.lookup("data_input_singleDSSD_analysis").c_str();
        output_dir = EUtil::getProjectRoot() + "/" + cfg.lookup("data_output_dir_DoubleDSSD").c_str();
    } else {
        setup_path = cfg.lookup("setup_file").c_str();
        target_path = cfg.lookup("target_file").c_str();
        input_path = cfg.lookup("data_input_singleDSSD_analysis").c_str();
        output_dir = cfg.lookup("data_output_dir_DoubleDSSD").c_str();
    }

    if (cfg.lookup("verbose")) {
        cout << "------------------------ IO configuration ------------------------" << endl;
        cout << "Setup:  " << setup_path << endl;
        cout << "Target: " << target_path << endl;
        cout << "Input:  " << input_path << endl;
        cout << "Output: " << output_dir << endl;
        cout << "------------------------------------------------------------------" << endl << endl;
    }
}


int main(int argc, char *argv[]) {
    prepareFileIO(EUtil::getProjectRoot() + "/Analysis.cfg");

    system(("mkdir -p " + output_dir).c_str());

    vector<string> input;
    int run;

    for (int i = 1; i < argc; i++) {
        run = stoi(argv[i]);
        findFilesMatchingWildcard(Form(input_path.c_str(), run), input);
    }

    for (auto &in: input) {
        clock_t start = clock();

        string stem = EUtil::getStem(in);
        TString outfile = (output_dir + "/" + stem + "lio.root").c_str();

        cout << "Reading from: " << in << endl;
        cout << "Printing to:  " << outfile << endl;

        TFile output(outfile, "RECREATE");
        cout << "starting analysis" << flush << endl;
        auto analysis = new SingleDssdAnalysis(outfile);
        cout << "making out tree" << flush << endl;
        analysis->findSingleDssdHits(in);
        cout << "closing file" << flush << endl;
        //output.Close();

        clock_t stop = clock();
        double elapsed = (double) (stop - start) / CLOCKS_PER_SEC;
        printf("\nFile: %s.\tTime elapsed: %.5f\n", (stem + "lio.root").c_str(), elapsed);
    }

    return 0;
}
