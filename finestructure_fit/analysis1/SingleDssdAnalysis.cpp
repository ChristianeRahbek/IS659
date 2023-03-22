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
    SingleDssdAnalysis(TString outfilename) {
        NUM = 0;

        //this->filePath = filePath;

        isBeta = 0;

        output = new TFile(outfilename, "RECREATE");
        t = new TTree("a", "a");
        t->Branch("mul", &mul);
        t->Branch("num", &NUM);
        t->Branch("dssdMul", &dssdMul);
        t->Branch("plasticMul", &plasticMul);

        v_dir = make_unique<DynamicBranchVector<TVector3>>(*t, "dir");
        v_pos = make_unique<DynamicBranchVector<TVector3>>(*t, "pos");

        v_theta = make_unique<DynamicBranchVector<double>>(*t, "theta", "mul");
        v_ang = make_unique<DynamicBranchVector<double>>(*t, "angle", "mul");

        v_Edep = make_unique<DynamicBranchVector<double>>(*t, "Edep", "mul");
        v_Ea = make_unique<DynamicBranchVector<double>>(*t, "Ea", "mul");
        v_BE = make_unique<DynamicBranchVector<double>>(*t, "BE", "mul");
        v_FE = make_unique<DynamicBranchVector<double>>(*t, "FE", "mul");
        v_Et = make_unique<DynamicBranchVector<double>>(*t, "Et", "mul");

        v_FT = make_unique<DynamicBranchVector<double>>(*t, "FT", "mul");
        v_BT = make_unique<DynamicBranchVector<double>>(*t, "BT", "mul");

        v_dE = make_unique<DynamicBranchVector<double>>(*t, "dE", "mul");
        v_Ecm = make_unique<DynamicBranchVector<double>>(*t, "Ecm", "mul");

        v_i = make_unique<DynamicBranchVector<short>>(*t, "id", "mul");

        v_isBeta = make_unique<DynamicBranchVector<int>>(*t, "isBeta", "mul");

        v_F = make_unique<DynamicBranchVector<short>>(*t, "FI", "mul");
        v_B = make_unique<DynamicBranchVector<short>>(*t, "BI", "mul");

        t->Branch("TPROTONS", &TPROTONS);
    }

    void findSingleDssdHits(string filePath) {
        cout << "making tree" << endl;
        auto *f= new TFile(filePath.c_str());
        auto *tr=(TTree*)f->Get("a");

        UInt_t maxHits = 20;
        Int_t numOrig;
        UInt_t mulOrig, TPROTONSOrig, plasticMulOrig, dssdMulOrig;
        Short_t idOrig[maxHits], FIOrig[maxHits], BIOrig[maxHits];
        Double_t EdepOrig[maxHits], EaOrig[maxHits], EtOrig[maxHits], BEOrig[maxHits], FEOrig[maxHits],
                    thetaOrig[maxHits], angOrig[maxHits], FTOrig[maxHits], BTOrig[maxHits], dEOrig[maxHits],
                    EcmOrig[maxHits];

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

        output->cd();

        cout << "looping through entries" << flush << endl;
        for(UInt_t i = 0; i < tr->GetEntries(); i++) { //looping through each entry
            clear();
            tr->GetEntry(i);

            if(dssdMulOrig != 1) continue; //this is a single DSSD analysis

            for(UInt_t j = 0; j < mulOrig; j++) {
                //cout << "i, j = " << i << ", " << j << flush << endl;

                double tol = 0; //FIX THIS
                double deltaT = 0; //FIX THIS

                if(plasticMulOrig > 0 && deltaT < tol) { //should plasticMul>1??
                    isBeta = 1; // true
                }
                else isBeta = 0; // false

                //cout << "adding to tree" << flush << endl;
                v_theta->add(thetaOrig[j]);
                v_ang->add(angOrig[j]);
                //cout << "angOrig added" << flush << endl;
                v_Et->add(EtOrig[j]);
                v_BE->add(BEOrig[j]);
                v_FE->add(FEOrig[j]);
                v_FT->add(FTOrig[j]);
                v_BT->add(BTOrig[j]);
                v_dE->add(dEOrig[j]);
                v_Ecm->add(EcmOrig[j]);
                v_Ea->add(EaOrig[j]);
                //cout << "EaOrig added" << flush << endl;
                v_Edep->add(EdepOrig[j]);
                //cout << "EdepOrig added" << flush << endl;
                v_i->add(idOrig[j]);
                v_isBeta->add(isBeta);
                v_F->add(FIOrig[j]);
                v_B->add(BIOrig[j]);


                //cout << "finished adding" << flush << endl;
                TPROTONS = TPROTONSOrig;
                plasticMul = plasticMulOrig;
                dssdMul = dssdMulOrig;
                NUM = numOrig;

                mul++;
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
                *v_Et, *v_Ea, *v_theta, *v_Edep,*v_i, *v_FE, *v_BE,
                *v_F, *v_B, *v_Ecm,*v_ang, *v_pos, *v_dir,
                *v_dE, *v_FT, *v_BT, *v_isBeta
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
    unique_ptr<DynamicBranchVector<TVector3>> v_dir, v_pos;
    //unique_ptr<DynamicBranchVector<double>> v_Edep;
    unique_ptr<DynamicBranchVector<double>> v_Edep, v_Ea, v_Et, v_BE, v_FE, v_theta, v_dE, v_Ecm;
    unique_ptr<DynamicBranchVector<short>> v_i;
    unique_ptr<DynamicBranchVector<int>> v_isBeta;
    unique_ptr<DynamicBranchVector<short>> v_F, v_B;
    unique_ptr<DynamicBranchVector<double>> v_ang;
    unique_ptr<DynamicBranchVector<double>> v_FT, v_BT;


    TFile *output;
    UInt_t mul{}, TPROTONS{}, dssdMul{}, plasticMul{};
    //string filePath;
    int isBeta;
};

string setup_path, target_path, input_path, output_dir;

void prepareFileIO(const string &configfile) {
    Config cfg;
    cfg.readFile(configfile.c_str());

    if (cfg.lookup("paths_relative_to_project_root")) {
        setup_path = EUtil::getProjectRoot() + "/" + cfg.lookup("setup_file").c_str();
        target_path = EUtil::getProjectRoot() + "/" + cfg.lookup("target_file").c_str();
        input_path = EUtil::getProjectRoot() + "/" + cfg.lookup("data_input_singleDSSD_analysis").c_str();
        output_dir = EUtil::getProjectRoot() + "/" + cfg.lookup("data_output_dir_singleDSSD").c_str();
    } else {
        setup_path = cfg.lookup("setup_file").c_str();
        target_path = cfg.lookup("target_file").c_str();
        input_path = cfg.lookup("data_input_singleDSSD_analysis").c_str();
        output_dir = cfg.lookup("data_output_dir_singleDSSD").c_str();
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
        cout << "starting analysis" << flush <<  endl;
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
