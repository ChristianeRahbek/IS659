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


using namespace std;
using namespace AUSA;
using namespace ROOT::Math;
using namespace TMath;
using namespace AUSA::Sort;
using namespace AUSA::EnergyLoss;
using namespace AUSA::Constants;
using namespace libconfig;


class MyAnalysis : public AbstractSortedAnalyzer {
public:
    MyAnalysis(Target &target, TFile *output) : target(target) {
        NUM = 0;

        t = new TTree("a", "a");
        t->Branch("mul", &mul);
        t->Branch("num", &NUM);

        v_dir = make_unique < DynamicBranchVector < TVector3 >> (*t, "dir");
        v_pos = make_unique < DynamicBranchVector < TVector3 >> (*t, "pos");

        v_theta = make_unique < DynamicBranchVector < double >> (*t, "theta", "mul");
        v_ang = make_unique < DynamicBranchVector < double >> (*t, "angle", "mul");

        v_E = make_unique < DynamicBranchVector < double >> (*t, "E", "mul");
        v_Ex = make_unique < DynamicBranchVector < double >> (*t, "Ex", "mul");
        v_BE = make_unique < DynamicBranchVector < double >> (*t, "BE", "mul");
        v_FE = make_unique < DynamicBranchVector < double >> (*t, "FE", "mul");

        v_FT = make_unique < DynamicBranchVector < double >> (*t, "FT", "mul");
        v_BT = make_unique < DynamicBranchVector < double >> (*t, "BT", "mul");

        v_dE = make_unique < DynamicBranchVector < double >> (*t, "dE", "mul");
        v_Ecm = make_unique < DynamicBranchVector < double >> (*t, "Ecm", "mul");

        v_i = make_unique < DynamicBranchVector < short >> (*t, "id", "mul");

        v_F = make_unique < DynamicBranchVector < short >> (*t, "FI", "mul");
        v_B = make_unique < DynamicBranchVector < short >> (*t, "BI", "mul");

        t->Branch("TPROTONS", &TPROTONS);


        string product = "He4"; //insert decay product here
        SiCalc = defaultRangeInverter(product, "Silicon"); //used to calculate energy loss in deadlayer

        for (auto &layer: target.getLayers()) {
            //this is used to calculate energy loss in target
            targetCalcs.push_back(defaultRangeInverter(Ion::predefined(product), layer.getMaterial()));
        }
    }

    void setup(const SortedSetupOutput &output) override {
        AbstractSortedAnalyzer::setup(output);
        for (size_t i = 0; i < output.dssdCount(); ++i) {
            auto dl = getFrontDeadLayer(output.getDssdOutput(i).detector());
            auto dlB = getBackDeadLayer(output.getDssdOutput(i).detector());
            deadlayerF.push_back(dl);
            deadlayerB.push_back(dlB);
        }
    }


    void analyze() override {
        clear();
        TPROTONS = output.getScalerOutput("TPROTONS").getValue();
        findHits();
        doAnalysis();
        if (mul > 0) { t->Fill(); }
        hits.clear();
        NUM++;
    }


    void findHits() {
        for (size_t i = 0; i < output.dssdCount(); i++) { //running through detectors
            auto &o = output.getDssdOutput(i);
            auto &p = output.getSingleOutput(i);
            auto &d = o.detector();
            auto m = AUSA::mul(o);

            for (UInt_t j = 0; j < m; j++) { //running through every hit in each detector
                Hit hit;

                auto dE = fEnergy(o, j) - bEnergy(o, j);
                hit.dE = dE;

                auto eDssd = energy(o, j);
                auto eFDssd = fEnergy(o, j);
                auto eBDssd = bEnergy(o, j);
                auto ePad = p.energy(0);


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
                } else {
                    hit.TF = 42;
                    hit.TB = 42;
                }

                auto angle = hit.direction.Angle(-d.getNormal());
                hit.angle = angle;

                auto tF = deadlayerF[i] / abs(cos(angle)); //actual thickness that the particle travels through, when it comes in at an angle.
                auto tB = deadlayerB[i] / abs(cos(angle));
                //auto tP = deadlayerP[i] / abs(cos(angle));

                double E = 0.0;
                double FE = 0.0;
                double BE = 0.0;

                E += eDssd;
                E += SiCalc->getTotalEnergyCorrection(E, tF);
                FE += eFDssd; //only energy correction on E.
                BE += eBDssd;

                hit.E = E;
                hit.BE = BE;
                hit.FE = FE;

                auto &t_c = target.getCenter();
                auto &from = position; //coordinates of where the particle hit the detector

                auto dir = t_c.Unit();
                auto t_thick = target.getThickness(); //given in unit of mm.


                /* stop_length has to be given in mm, so it is also 0.1221 um
                 * It is calculated at eloss.kern.phys.au.dk for E_beam = 30 keV*/
                double stop_length = 0.1221* pow(10,-3); //how far the beam goes to be stopped
                //double stop_length = 7*t_thick/8; //how far the beam goes to be stopped

                /*
                 * Med dette kode som det er lige nu, så er upstreamdetektorerne ens og downstream detektorerne ens.
                 * Downstream detektorerne er overkorregeret en smule, og MÅSKE er upstream detektorerne også...
                 * */

                double transversed_extra = t_thick/2 - stop_length;
                auto stop_coord = t_c + direction*transversed_extra; //coordinate where the beam is stopped
                //auto stop_coord = t_c - TVector3(0,0,transversed_extra); //coordinate where the beam is stopped


/*                auto &calc = targetCalcs[0];
                if(i == 2 || i == 3){ //if downstream
                    auto traveled = t_thick-stop_length;
                    hit.E += calc->getTotalEnergyCorrection(hit.E, traveled/abs(cos(from.Angle(t_c))));
                }
                else hit.E += calc->getTotalEnergyCorrection(hit.E, stop_length/abs(cos(from.Angle(t_c))));
*/

                for (auto &intersection: target.getIntersections(from, stop_coord)) {
                    auto &calc = targetCalcs[intersection.index];
                    //hit.E += calc->getTotalEnergyCorrection(hit.E, intersection.transversed);

                    if(i == 1 || i == 2){ //if downstream
                        auto traveled = t_thick-stop_length;
                        hit.E += calc->getTotalEnergyCorrection(hit.E, traveled/abs(cos(from.Angle(t_c))));
                    }
                    else hit.E += calc->getTotalEnergyCorrection(hit.E, stop_length/abs(cos(from.Angle(t_c))));
                }

                hit.index = i;
                hit.lVector = {sqrt(2 * hit.E * ALPHA_MASS) * hit.direction, hit.E + ALPHA_MASS};
                hits.emplace_back(move(hit));
            }
        }
    }

    void doAnalysis() {
        if (hits.empty()) return;
        mul = hits.size();

        double alpha_sep = 4730; //keV from

        for (auto &hit: hits) {
            v_pos->add(hit.position);
            v_dir->add(hit.direction);
            v_theta->add(hit.theta * TMath::RadToDeg());

            auto Ex = hit.E * 5 / 4 + alpha_sep; //finding excitation energy

            v_E->add(hit.E);
            v_Ex->add(Ex);
            v_BE->add(hit.BE);
            v_FE->add(hit.FE);

            v_dE->add(hit.dE);
            //v_Ecm->add(E0cm);
            v_ang->add(hit.angle * TMath::RadToDeg());

            v_i->add(static_cast<short>(hit.index));
            v_F->add(hit.fseg);
            v_B->add(hit.bseg);
            v_FT->add(hit.TF);
            v_BT->add(hit.TB);
        }
    }

    void terminate() override {
        AbstractSortedAnalyzer::terminate();
        gDirectory->WriteTObject(t);
    }

    void clear() {
        mul = 0;
        AUSA::clear(
                *v_E, *v_theta, *v_Ex,
                *v_i, *v_FE, *v_BE,
                *v_F, *v_B, *v_Ecm,
                *v_ang, *v_pos, *v_dir,
                *v_dE, *v_FT, *v_BT
        );
    }

    int NUM;
    TTree *t;
    unique_ptr<DynamicBranchVector<TVector3>> v_dir, v_pos;
    unique_ptr<DynamicBranchVector<double>> v_E, v_BE, v_FE, v_theta, v_dE, v_Ecm, v_Ex;
    unique_ptr<DynamicBranchVector<short>> v_i;
    unique_ptr<DynamicBranchVector<short>> v_F, v_B;
    unique_ptr<DynamicBranchVector<double>> v_ang;
    unique_ptr<DynamicBranchVector<double>> v_FT, v_BT;

    UInt_t mul{}, TPROTONS{}; //TPATTERN{}, EGPS{};
    vector<Hit> hits;

    unique_ptr<EnergyLossRangeInverter> SiCalc;
    vector<unique_ptr<EnergyLossRangeInverter>> targetCalcs;
    vector<double> deadlayerF, deadlayerB, deadlayerP;
    Target &target;
    bool simulation = false;
};

string setup_path, target_path, input_path, output_dir;
void prepareFileIO(const string &configfile) {
    Config cfg;
    cfg.readFile(configfile.c_str());

    if (cfg.lookup("paths_relative_to_project_root")) {
        setup_path = ANALYSIS::getProjectRoot() + "/" + cfg.lookup("setup_file").c_str();
        target_path = ANALYSIS::getProjectRoot() + "/" + cfg.lookup("target_file").c_str();
        input_path = ANALYSIS::getProjectRoot() + "/" + cfg.lookup("data_input").c_str();
        output_dir = ANALYSIS::getProjectRoot() + "/" + cfg.lookup("data_output_dir").c_str();
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
    prepareFileIO(ANALYSIS::getProjectRoot() + "/Analysis.cfg");
    auto setup = JSON::readSetupFromJSON(setup_path);
    auto target = JSON::readTargetFromJSON(target_path);

    system(("mkdir -p " + output_dir).c_str());

    vector <string> input;
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

        string stem = ANALYSIS::getStem(in);
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