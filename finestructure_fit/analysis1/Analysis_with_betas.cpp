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

class MyAnalysis : public AbstractSortedAnalyzer {
public:
    MyAnalysis(Target &target, TFile *output) : target(target) {
        NUM = 0;

        t = new TTree("a", "a");
        t->Branch("mul", &mul);
        t->Branch("num", &NUM);

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

        v_posX = make_unique<DynamicBranchVector<double>>(*t, "posX", "mul");
        v_posY = make_unique<DynamicBranchVector<double>>(*t, "posY", "mul");
        v_posZ = make_unique<DynamicBranchVector<double>>(*t, "posZ", "mul");

        v_i = make_unique<DynamicBranchVector<short>>(*t, "id", "mul");

        v_F = make_unique<DynamicBranchVector<short>>(*t, "FI", "mul");
        v_B = make_unique<DynamicBranchVector<short>>(*t, "BI", "mul");

        t->Branch("TPROTONS", &TPROTONS);

        aSiCalc = defaultRangeInverter("a", "Silicon"); //for alphas through deadlayer
        tSiCalc = defaultRangeInverter("t", "Silicon"); //for tritons through deadlayer

        for (auto &layer: target.getLayers()) {
            aTargetCalcs.push_back(
                    defaultRangeInverter(Ion::predefined("a"), layer.getMaterial())); //for alphas through target
            tTargetCalcs.push_back(
                    defaultRangeInverter(Ion::predefined("t"), layer.getMaterial())); //for tritons through target
        }
    }

    void setup(const SortedSetupOutput &output) override {
        AbstractSortedAnalyzer::setup(output);

        double dl, dlB;
        for(const auto &entry : detector_map) {
            auto id = entry.first; //getting inputs of detector
            auto det = entry.second;

            if(det->getDetType() == IS659Type::SquareDSSD) {
                dl = getFrontDeadLayer(output.getDssdOutput(id).detector());
                dlB = getBackDeadLayer(output.getDssdOutput(id).detector());
                deadlayerF.push_back(dl);
                deadlayerB.push_back(dlB);
            }
            else{
                //no deadlayer corrections in Plastics
            }
        }
    }


    void analyze() override {
        clear();
        TPROTONS = output.getScalerOutput("TPROTONS").getValue();
        findHits();
        doAnalysis();
        if (mul > 0) { t->Fill(); } //looking for coincidences
        hits.clear();
        NUM++;
    }

    void findHits() {
        for(const auto &entry : detector_map) {
            auto det = entry.second;
            auto detType = det->getDetType();
            if(detType == IS659Type::SquareDSSD) findDSSDHit(det);
            else findPlasticHit(det);
        }
    }

    void findDSSDHit(IS659Detector* detector) {
        auto name = detector->getName();
        auto &o = output.getDssdOutput(name);
        auto &d = o.detector();
        auto m = AUSA::mul(o);

        for (UInt_t j = 0; j < m; j++) {
            Hit hit;

            auto id = detector->getId();

            auto dE = fEnergy(o, j) - bEnergy(o, j);
            hit.dE = dE;

            auto Edep = energy(o, j);
            auto eFDssd = fEnergy(o, j);
            auto eBDssd = bEnergy(o, j);
            //auto ePlastic = p.energy(0);

            auto BI = bSeg(o, j);
            auto FI = fSeg(o, j);
            hit.fseg = short(FI);
            hit.bseg = short(BI);

            auto position = o.detector().getUniformPixelPosition(FI, BI);
            auto origin = target.getCenter();
            hit.position = position;
            auto direction = (position - origin).Unit();
            hit.direction = direction;
            hit.theta = hit.direction.Theta();

            if (!simulation) {
                hit.TF = fTime(o, j);
                hit.TB = bTime(o, j);
                //hit.TPlastic = p.time(0);
            } else {
                hit.TF = 42;
                hit.TB = 42;
                //hit.TPlastic = 42;
            }

            auto angle = hit.direction.Angle(-d.getNormal());
            hit.angle = angle;

            auto tF = deadlayerF[id] / abs(cos(angle));
            auto tB = deadlayerB[id] / abs(cos(angle));

            double Ea = 0.0;
            double Et = 0.0;
            double FE = 0.0;
            double BE = 0.0;

            /* Energy corrections in deadlayer */
            Ea += Edep;
            Ea += aSiCalc->getTotalEnergyCorrection(Ea, tF);
            Et += Edep;
            Et += tSiCalc->getTotalEnergyCorrection(Et, tF);

            FE += eFDssd; //only energy correction on Ea and Et.
            BE += eBDssd;

            /* stop_length has to be given in mm, so it is also 0.1221 um
             * It is calculated at eloss.kern.phys.au.dk for E_beam = 30 keV*/
            double stop_length = 0.2868 * pow(10, -3); //how far the beam goes to be stopped

            /* Energy corrections in target */
            auto &from = position;
            for (auto &intersection: target.getIntersections(from, target.getCenter() /*NOT IN CENTER!*/)) {
                auto &calca = aTargetCalcs[intersection.index];
                auto &calct = tTargetCalcs[intersection.index];
                auto traveled = target.getThickness() - stop_length;
                if (id == 1 || id == 2) { //downstream detectorer
                    Ea += calca->getTotalEnergyCorrection(Ea, traveled / abs(cos(from.Angle(target.getCenter()))));
                    Et += calct->getTotalEnergyCorrection(Et, traveled / abs(cos(from.Angle(target.getCenter()))));
                } else {
                    Ea += calca->getTotalEnergyCorrection(Ea,
                                                          stop_length / abs(cos(from.Angle(target.getCenter()))));
                    Et += calct->getTotalEnergyCorrection(Et,
                                                          stop_length / abs(cos(from.Angle(target.getCenter()))));
                }
            }

            hit.Ea = Ea;
            hit.Et = Et;
            hit.Edep = Edep; //deposited energy in detector
            hit.BE = BE;
            hit.FE = FE;

            hit.index = id;

            hits.emplace_back(move(hit));
        }
    }

    void findPlasticHit(IS659Detector* detector) {
        auto name = detector->getName();
        auto &o = output.getSingleOutput(name);
        auto &d = o.detector();
        auto m = AUSA::mul(o);

        for (UInt_t j = 0; j < m; j++) {
            Hit hit;

            auto id = detector->getId();

            auto Edep = o.energy(j);

            hit.fseg = short(o.segment(j)); //there are two front IDs
            hit.bseg = short(NAN);

            auto origin = target.getCenter();
            hit.position = NAN_TVECTOR3;
            hit.direction = NAN_TVECTOR3;
            hit.theta = NAN;

            if (!simulation) {
                hit.TF = o.time(j);
                hit.TB = NAN;
            } else {
                hit.TF = 42;
                hit.TB = 42;
            }

            hit.angle = NAN;

            hit.Ea = NAN;
            hit.Et = NAN;
            hit.Edep = Edep; //deposited energy in detector
            hit.BE = NAN;
            hit.FE = NAN;
            hit.dE = NAN;

            hit.index = id;

            hits.emplace_back(move(hit));
        }
    }

    void doAnalysis() {
        auto mult = hits.size();

        if(hits.empty()) return;
        //if (mult < 2) return;

        for (size_t i = 0; i < mult; i++) {
            auto h = hits[i];

            auto posi = h.position;
            v_pos->add(posi);
            v_posX->add(posi.X());
            v_posY->add(posi.Y());
            v_posZ->add(posi.Z());
            v_dir->add(h.direction);
            v_theta->add(h.theta * TMath::RadToDeg());

            v_Ea->add(h.Ea);
            v_Et->add(h.Et);
            v_Edep->add(h.Edep);
            v_BE->add(h.BE);
            v_FE->add(h.FE);

            v_dE->add(h.dE);
            //v_Ecm->add(Ecm);
            v_ang->add(h.angle * TMath::RadToDeg());


            v_i->add(static_cast<short>(h.index));
            v_F->add(h.fseg);
            v_B->add(h.bseg);
            v_FT->add(h.TF);
            v_BT->add(h.TB);

            mul++;
        }
    }

    void terminate() override {
        AbstractSortedAnalyzer::terminate();
        gDirectory->WriteTObject(t);
    }

    void clear() {
        mul = 0;
        AUSA::clear(
                *v_dir, *v_pos, *v_Edep, *v_Ea, *v_Et, *v_BE,
                *v_FE, *v_theta, *v_dE, *v_Ecm, *v_i, *v_F, *v_B,
                *v_ang, *v_FT, *v_BT, *v_posX, *v_posY, *v_posZ
        );
    }

    const map<unsigned short, IS659Detector*> detector_map = {
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
    unique_ptr<DynamicBranchVector<double>> v_Edep;
    unique_ptr<DynamicBranchVector<double>> v_posX, v_posY, v_posZ;
    unique_ptr<DynamicBranchVector<double>> v_Ea, v_Et, v_BE, v_FE, v_theta, v_dE, v_Ecm;;
    unique_ptr<DynamicBranchVector<short>> v_i, v_F, v_B;
    unique_ptr<DynamicBranchVector<double>> v_ang;
    unique_ptr<DynamicBranchVector<double>> v_FT, v_BT;


    UInt_t mul{}, TPROTONS{};
    vector<Hit> hits;

    unique_ptr<EnergyLossRangeInverter> aSiCalc, tSiCalc;
    vector<unique_ptr<EnergyLossRangeInverter>> aTargetCalcs, tTargetCalcs;
    vector<double> deadlayerF, deadlayerB, deadlayerP;
    Target &target;
    ParticleType *TRITON;
    TLorentzVector beamVector;
    double beamEnergy;
    TVector3 cmBoost;
    bool simulation = false;
};

string setup_path, target_path, input_path, output_dir;

void prepareFileIO(const string &configfile) {
    Config cfg;
    cfg.readFile(configfile.c_str());

    if (cfg.lookup("paths_relative_to_project_root")) {
        setup_path = EUtil::getProjectRoot() + "/" + cfg.lookup("setup_file").c_str();
        target_path = EUtil::getProjectRoot() + "/" + cfg.lookup("target_file").c_str();
        input_path = EUtil::getProjectRoot() + "/" + cfg.lookup("data_input").c_str();
        output_dir = EUtil::getProjectRoot() + "/" + cfg.lookup("data_output_dir").c_str();
    } else {
        setup_path = cfg.lookup("setup_file").c_str();
        target_path = cfg.lookup("target_file").c_str();
        input_path = cfg.lookup("data_input").c_str();
        output_dir = cfg.lookup("data_output_dir").c_str();
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
    auto setup = JSON::readSetupFromJSON(setup_path);
    auto target = JSON::readTargetFromJSON(target_path);

    system(("mkdir -p " + output_dir).c_str());

    vector<string> input;
    int run;

    for (int i = 1; i < argc; i++) {
        run = stoi(argv[i]);
        findFilesMatchingWildcard(Form(input_path.c_str(), run), input);
    }

    for (auto &in: input) {
        clock_t start = clock();

        SortedReader reader{*setup};
        reader.add(in);
        reader.setVerbose(true);

        string stem = EUtil::getStem(in);
        TString outfile = (output_dir + "/" + stem + "lio.root").c_str();

        cout << "Reading from: " << in << endl;
        cout << "Printing to:  " << outfile << endl;

        TFile output(outfile, "RECREATE");
        auto analysis = make_shared<MyAnalysis>(target, &output);
        reader.attach(analysis);
        reader.run();

        clock_t stop = clock();
        double elapsed = (double) (stop - start) / CLOCKS_PER_SEC;
        printf("\nFile: %s.\tTime elapsed: %.5f\n", (stem + "lio.root").c_str(), elapsed);
    }

    return 0;
}
