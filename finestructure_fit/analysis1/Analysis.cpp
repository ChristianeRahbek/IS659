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

    //beamVector =

    t = new TTree("a", "a");
    t->Branch("mul", &mul);
    t->Branch("num", &NUM);

    v_dir0 = make_unique<DynamicBranchVector<TVector3>>(*t, "dir0");
    v_pos0 = make_unique<DynamicBranchVector<TVector3>>(*t, "pos0");
    v_dir1 = make_unique<DynamicBranchVector<TVector3>>(*t, "dir1");
    v_pos1 = make_unique<DynamicBranchVector<TVector3>>(*t, "pos1");

    v_theta0 = make_unique<DynamicBranchVector<double>>(*t, "theta0", "mul");
    v_ang0 = make_unique<DynamicBranchVector<double>>(*t, "angle0", "mul");
    v_theta1 = make_unique<DynamicBranchVector<double>>(*t, "theta1", "mul");
    v_ang1 = make_unique<DynamicBranchVector<double>>(*t, "angle1", "mul");

    v_E0 = make_unique<DynamicBranchVector<double>>(*t, "E0", "mul");
    v_BE0 = make_unique<DynamicBranchVector<double>>(*t, "BE0", "mul");
    v_FE0 = make_unique<DynamicBranchVector<double>>(*t, "FE0", "mul");
    v_E1 = make_unique<DynamicBranchVector<double>>(*t, "E1", "mul");
    v_BE1 = make_unique<DynamicBranchVector<double>>(*t, "BE1", "mul");
    v_FE1 = make_unique<DynamicBranchVector<double>>(*t, "FE1", "mul");

    v_FT0 = make_unique<DynamicBranchVector<double>>(*t, "FT0", "mul");
    v_BT0 = make_unique<DynamicBranchVector<double>>(*t, "BT0", "mul");
    v_FT1 = make_unique<DynamicBranchVector<double>>(*t, "FT1", "mul");
    v_BT1 = make_unique<DynamicBranchVector<double>>(*t, "BT1", "mul");

    v_dE0 = make_unique<DynamicBranchVector<double>>(*t, "dE0", "mul");
    v_Ecm0 = make_unique<DynamicBranchVector<double>>(*t, "E0cm", "mul");
    v_dE1 = make_unique<DynamicBranchVector<double>>(*t, "dE1", "mul");
    v_Ecm1 = make_unique<DynamicBranchVector<double>>(*t, "E1cm", "mul");

    v_i0 = make_unique<DynamicBranchVector<short>>(*t, "id0", "mul");
    v_i1 = make_unique<DynamicBranchVector<short>>(*t, "id1", "mul");

    v_F0 = make_unique<DynamicBranchVector<short>>(*t, "FI0", "mul");
    v_B0 = make_unique<DynamicBranchVector<short>>(*t, "BI0", "mul");
    v_F1 = make_unique<DynamicBranchVector<short>>(*t, "FI1", "mul");
    v_B1 = make_unique<DynamicBranchVector<short>>(*t, "BI1", "mul");

    //t->Branch("TPATTERN", &TPATTERN);
    t->Branch("TPROTONS", &TPROTONS);
    //t->Branch("EGPS", &EGPS);

    //SiCalc = defaultRangeInverter("p", "Silicon");
    SiCalc = defaultRangeInverter("8He", "Silicon");

    for (auto &layer: target.getLayers()) {
      //targetCalcs.push_back(defaultRangeInverter(Ion::predefined("p"), layer.getMaterial()));
        targetCalcs.push_back(defaultRangeInverter(Ion::predefined("8He"), layer.getMaterial()));
    }
  }

  void setup(const SortedSetupOutput &output) override {
    AbstractSortedAnalyzer::setup(output);
    for (size_t i = 0; i < output.dssdCount(); ++i) {
      auto dl = getFrontDeadLayer(output.getDssdOutput(i).detector());
      auto dlB = getBackDeadLayer(output.getDssdOutput(i).detector());
      //auto dlP = getFrontDeadLayer(output.getSingleOutput(i).detector());
      deadlayerF.push_back(dl);
      deadlayerB.push_back(dlB);
      //deadlayerP.push_back(dlP);
    }
  }


  void analyze() override {
    clear();
    //TPATTERN = output.getScalerOutput("TPATTERN").getValue();
    TPROTONS = output.getScalerOutput("TPROTONS").getValue();
    //EGPS = output.getScalerOutput("EGPS").getValue();
    findHits();
    doAnalysis();
    if (mul > 0) { t->Fill(); }
    hits.clear();
    NUM++;
  }


  void findHits() {
    for (size_t i = 0; i < output.dssdCount(); i++) {
      auto &o = output.getDssdOutput(i);
      auto &p = output.getSingleOutput(i);
      auto &d = o.detector();
      auto m = AUSA::mul(o);

      for (UInt_t j = 0; j < m; j++) {
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
          //hit.TPad = p.time(0);
        } else {
          hit.TF = 42;
          hit.TB = 42;
          //hit.TPad = 42;
        }

        auto angle = hit.direction.Angle(-d.getNormal());
        hit.angle = angle;

        auto tF = deadlayerF[i] / abs(cos(angle));
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


        auto &from = position;
        for (auto &intersection: target.getIntersections(from, target.getCenter())) {
          auto &calc = targetCalcs[intersection.index];
          hit.E += calc->getTotalEnergyCorrection(hit.E, intersection.transversed);
        }

        hit.index = i;
        hit.lVector = {sqrt(2 * hit.E * ALPHA_MASS) * hit.direction, hit.E + ALPHA_MASS};
        hits.emplace_back(move(hit));
      }
    }
  }

  void doAnalysis() {
    if (hits.empty()) return;
    mul = hits.size(); //mul0, mul1 or mul??? It started as mul.

    for(auto &hit: hits) {
        v_pos0->add(hit.position);
        v_dir0->add(hit.direction);
        v_theta0->add(hit.theta * TMath::RadToDeg());

        v_E0->add(hit.E);
        v_BE0->add(hit.BE);
        v_FE0->add(hit.FE);

        v_dE0->add(hit.dE);
            //v_E0cm->add(E0cm);
        v_ang0->add(hit.angle * TMath::RadToDeg());

        v_i0->add(static_cast<short>(hit.index));
        v_F0->add(hit.fseg);
        v_B0->add(hit.bseg);
        v_FT0->add(hit.TF);
        v_BT0->add(hit.TB);
    }
  }

  void terminate() override {
    AbstractSortedAnalyzer::terminate();
    gDirectory->WriteTObject(t);
  }

  void clear() {
    mul = 0;
    AUSA::clear(
        *v_E1, *v_theta1, *v_E0, *v_theta0,
        *v_i1, *v_FE1, *v_BE1, *v_i0, *v_FE0, *v_BE0,
        *v_F1, *v_B1, *v_Ecm1, *v_F0, *v_B0, *v_Ecm0,
        *v_ang1, *v_pos1, *v_dir1, *v_ang0, *v_pos0, *v_dir0,
        *v_dE1, *v_FT1, *v_BT1, *v_dE0, *v_FT0, *v_BT0
    );
  }

  int NUM;
  TTree *t;
  unique_ptr<DynamicBranchVector<TVector3>> v_dir1, v_pos1, v_dir0, v_pos0;
  unique_ptr<DynamicBranchVector<double>> v_E0, v_BE0, v_FE0, v_theta0, v_dE0, v_Ecm0;
  unique_ptr<DynamicBranchVector<double>> v_E1, v_BE1, v_FE1, v_theta1, v_dE1, v_Ecm1;
  unique_ptr<DynamicBranchVector<short>> v_i1, v_i0;
  unique_ptr<DynamicBranchVector<short>> v_F1, v_B1, v_F0, v_B0;
  unique_ptr<DynamicBranchVector<double>> v_ang1, v_ang0;
  unique_ptr<DynamicBranchVector<double>> v_FT1, v_BT1, v_FT0, v_BT0;

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