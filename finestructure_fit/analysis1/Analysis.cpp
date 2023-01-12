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

    v_dir = make_unique<DynamicBranchVector<TVector3>>(*t, "dir");
    v_pos = make_unique<DynamicBranchVector<TVector3>>(*t, "pos");

    v_theta = make_unique<DynamicBranchVector<double>>(*t, "theta", "mul");
    v_ang = make_unique<DynamicBranchVector<double>>(*t, "angle", "mul");

    v_Edssd = make_unique<DynamicBranchVector<double>>(*t, "Edssd", "mul");
    v_Ea = make_unique<DynamicBranchVector<double>>(*t, "Ea", "mul");
    v_BE = make_unique<DynamicBranchVector<double>>(*t, "BE", "mul");
    v_FE = make_unique<DynamicBranchVector<double>>(*t, "FE", "mul");
    v_Et = make_unique<DynamicBranchVector<double>>(*t, "Et", "mul");

    v_Edep_alphas = make_unique<DynamicBranchVector<double>>(*t, "Edep_alphas", "mul");

    v_Edep0 = make_unique<DynamicBranchVector<double>>(*t, "Edep0", "mul");
    v_Edep1 = make_unique<DynamicBranchVector<double>>(*t, "Edep1", "mul");
    v_Edep2 = make_unique<DynamicBranchVector<double>>(*t, "Edep2", "mul");
    v_Edep3 = make_unique<DynamicBranchVector<double>>(*t, "Edep3", "mul");

    /*
    v_Ea0 = make_unique<DynamicBranchVector<double>>(*t, "Ea0", "mul");
    v_Ea1 = make_unique<DynamicBranchVector<double>>(*t, "Ea1", "mul");
    v_Ea2 = make_unique<DynamicBranchVector<double>>(*t, "Ea2", "mul");
    v_Ea3 = make_unique<DynamicBranchVector<double>>(*t, "Ea3", "mul");

    v_Et0 = make_unique<DynamicBranchVector<double>>(*t, "Et0", "mul");
    v_Et1 = make_unique<DynamicBranchVector<double>>(*t, "Et1", "mul");
    v_Et2 = make_unique<DynamicBranchVector<double>>(*t, "Et2", "mul");
    v_Et3 = make_unique<DynamicBranchVector<double>>(*t, "Et3", "mul");
    */

    v_FT = make_unique<DynamicBranchVector<double>>(*t, "FT", "mul");
    v_BT = make_unique<DynamicBranchVector<double>>(*t, "BT", "mul");

    v_dE = make_unique<DynamicBranchVector<double>>(*t, "dE", "mul");
    v_Ecm = make_unique<DynamicBranchVector<double>>(*t, "Ecm", "mul");

    v_i = make_unique<DynamicBranchVector<short>>(*t, "id", "mul");

    v_F = make_unique<DynamicBranchVector<short>>(*t, "FI", "mul");
    v_B = make_unique<DynamicBranchVector<short>>(*t, "BI", "mul");

    t->Branch("TPROTONS", &TPROTONS);

    aSiCalc = defaultRangeInverter("a", "Silicon"); //for alphas through deadlayer
    tSiCalc = defaultRangeInverter("t", "Silicon"); //for tritons through deadlayer

    for (auto &layer: target.getLayers()) {
        aTargetCalcs.push_back(defaultRangeInverter(Ion::predefined("a"), layer.getMaterial())); //for alphas through target
        tTargetCalcs.push_back(defaultRangeInverter(Ion::predefined("t"), layer.getMaterial())); //for tritons through target
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
    if (mul > 0) { t->Fill(); } //looking for coincidences
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

        double Ea = 0.0;
        double Et = 0.0;
        double FE = 0.0;
        double BE = 0.0;

        /* Energy corrections in deadlayer */
        Ea += eDssd;
        Ea += aSiCalc->getTotalEnergyCorrection(Ea, tF);
        Et += eDssd;
        Et += tSiCalc->getTotalEnergyCorrection(Et, tF);

        FE += eFDssd; //only energy correction on Ea and Et.
        BE += eBDssd;


        /* stop_length has to be given in mm, so it is also 0.1221 um
         * It is calculated at eloss.kern.phys.au.dk for E_beam = 30 keV*/
        double stop_length = 0.2868* pow(10,-3); //how far the beam goes to be stopped



        /* Energy corrections in target */
        auto &from = position;
        for (auto &intersection: target.getIntersections(from, target.getCenter() /*NOT IN CENTER!*/)) {
            auto &calca = aTargetCalcs[intersection.index];
            auto &calct = tTargetCalcs[intersection.index];
            auto traveled = target.getThickness() - stop_length;
            if(i == 2 || i == 3) {
                Ea += calca->getTotalEnergyCorrection(Ea, traveled/abs(cos(from.Angle(target.getCenter()))));
                Et += calct->getTotalEnergyCorrection(Et, traveled/abs(cos(from.Angle(target.getCenter()))));
            }
            else {
                Ea += calca->getTotalEnergyCorrection(Ea, stop_length/abs(cos(from.Angle(target.getCenter()))));
                Et += calct->getTotalEnergyCorrection(Et, stop_length/abs(cos(from.Angle(target.getCenter()))));
            }
        }

        hit.Ea = Ea;
        hit.Et = Et;
        hit.Edssd = eDssd; //deposited energy in detector
        hit.BE = BE;
        hit.FE = FE;

        //Saving energies from different detectors
        if(i == 0) {
            hit.Edep0 = eDssd;
            hit.Ea0 = Ea;
            hit.Et0 = Et;
        }
        else if (i == 1){
            hit.Edep1 = eDssd;
            hit.Ea1 = Ea;
            hit.Et1 = Et;
        }
        else if (i == 2){
            hit.Edep2 = eDssd;
            hit.Ea2 = Ea;
            hit.Et2 = Et;
        }
        else if (i == 3){
            hit.Edep3 = eDssd;
            hit.Ea3 = Ea;
            hit.Et3 = Et;
        }

        //finding only


        hit.index = i;
        hits.emplace_back(move(hit));
      }
    }
  }

  void doAnalysis() {
    if (hits.empty()) return;
    //mul = hits.size();

    for(auto &hit: hits) {
        v_pos->add(hit.position);
        v_dir->add(hit.direction);
        v_theta->add(hit.theta * TMath::RadToDeg());

        v_Ea->add(hit.Ea);
        v_Et->add(hit.Et);
        v_Edssd -> add(hit.Edssd);
        v_BE->add(hit.BE);
        v_FE->add(hit.FE);

        v_Edep0 -> add(hit.Edep0);
        v_Edep1 -> add(hit.Edep1);
        v_Edep2 -> add(hit.Edep2);
        v_Edep3 -> add(hit.Edep3);

        v_Edep_alphas -> add(hit.Edep_alphas);

        v_Ea0 -> add(hit.Ea0);
        v_Ea1 -> add(hit.Ea1);
        v_Ea2 -> add(hit.Ea2);
        v_Ea3 -> add(hit.Ea3);

        v_Et0 -> add(hit.Ea0);
        v_Et1 -> add(hit.Ea1);
        v_Et2 -> add(hit.Ea2);
        v_Et3 -> add(hit.Ea3);

        v_dE->add(hit.dE);
            //v_Ecm->add(Ecm);
        v_ang->add(hit.angle * TMath::RadToDeg());

        v_i->add(static_cast<short>(hit.index));
        v_F->add(hit.fseg);
        v_B->add(hit.bseg);
        v_FT->add(hit.TF);
        v_BT->add(hit.TB);

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
        *v_Et, *v_Ea, *v_theta, *v_Edssd,
        *v_Edep0, *v_Edep1, *v_Edep2, *v_Edep3,
        *v_Ea0, *v_Ea1, *v_Ea2, *v_Ea3,
        *v_Et0, *v_Et1, *v_Et2, *v_Et3,
        *v_i, *v_FE, *v_BE,
        *v_F, *v_B, *v_Ecm,
        *v_ang, *v_pos, *v_dir,
        *v_dE, *v_FT, *v_BT, *v_Edep_alphas
    );
  }

  int NUM;
  TTree *t;
  unique_ptr<DynamicBranchVector<TVector3>> v_dir, v_pos;
  unique_ptr<DynamicBranchVector<double>> v_Edssd, v_Edep0, v_Edep1, v_Edep2, v_Edep3, v_Edep_alphas;
  unique_ptr<DynamicBranchVector<double>> v_Ea, v_Et, v_BE, v_FE, v_theta, v_dE, v_Ecm;
  unique_ptr<DynamicBranchVector<double>> v_Ea0, v_Ea1, v_Ea2, v_Ea3, v_Et0, v_Et1, v_Et2, v_Et3;
  unique_ptr<DynamicBranchVector<short>> v_i;
  unique_ptr<DynamicBranchVector<short>> v_F, v_B;
  unique_ptr<DynamicBranchVector<double>> v_ang;
  unique_ptr<DynamicBranchVector<double>> v_FT, v_BT;

  UInt_t mul{}, TPROTONS{};
  vector<Hit> hits;

  unique_ptr<EnergyLossRangeInverter> aSiCalc, tSiCalc;
  vector<unique_ptr<EnergyLossRangeInverter>> aTargetCalcs, tTargetCalcs;
  vector<double> deadlayerF, deadlayerB, deadlayerP;
  Target &target;
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