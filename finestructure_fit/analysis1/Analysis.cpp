/** Updated by Christiane Rahbek.
 *  This Analysis finds 2 coincidential instances in the DSSDs.
 *  No other detectors such as Plastics and Clovers are considered.
**/

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

    beamEnergy = 30; //keV
    beamVector = constructBeamVector(Ion(2,8), Ion(6,12), beamEnergy);
    cmBoost = -beamVector.BoostVector();

    TRITON = new ParticleType("H3");

    t = new TTree("a", "a");
    t->Branch("mul", &mul);
    t->Branch("num", &NUM);

    v_dir0 = make_unique<DynamicBranchVector<TVector3>>(*t, "dir0");
    v_pos0 = make_unique<DynamicBranchVector<TVector3>>(*t, "pos0");
    v_dir1 = make_unique<DynamicBranchVector<TVector3>>(*t, "dir1");
    v_pos1 = make_unique<DynamicBranchVector<TVector3>>(*t, "pos1");

    v_pnat = make_unique<DynamicBranchVector<TVector3>>(*t, "pnat");
    v_pnta = make_unique<DynamicBranchVector<TVector3>>(*t, "pnta");
    v_Enat = make_unique<DynamicBranchVector<double>>(*t, "Enat", "mul");
    v_Enta = make_unique<DynamicBranchVector<double>>(*t, "Enta", "mul");

    v_theta0 = make_unique<DynamicBranchVector<double>>(*t, "theta0", "mul");
    v_ang0 = make_unique<DynamicBranchVector<double>>(*t, "angle0", "mul");
    v_theta1 = make_unique<DynamicBranchVector<double>>(*t, "theta1", "mul");
    v_ang1 = make_unique<DynamicBranchVector<double>>(*t, "angle1", "mul");

    v_hitAng = make_unique<DynamicBranchVector<double>>(*t, "hitAng", "mul");

    v_Edssd0 = make_unique<DynamicBranchVector<double>>(*t, "Edssd0", "mul");
    v_Ea0 = make_unique<DynamicBranchVector<double>>(*t, "Ea0", "mul");
    v_BE0 = make_unique<DynamicBranchVector<double>>(*t, "BE0", "mul");
    v_FE0 = make_unique<DynamicBranchVector<double>>(*t, "FE0", "mul");
    v_Et0 = make_unique<DynamicBranchVector<double>>(*t, "Et0", "mul");
    v_Edssd1 = make_unique<DynamicBranchVector<double>>(*t, "Edssd1", "mul");
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

    v_F0 = make_unique<DynamicBranchVector<short>>(*t, "FI0", "mul");
    v_B0 = make_unique<DynamicBranchVector<short>>(*t, "BI0", "mul");
    v_F1 = make_unique<DynamicBranchVector<short>>(*t, "FI1", "mul");
    v_B1 = make_unique<DynamicBranchVector<short>>(*t, "BI1", "mul");


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
    for (size_t i = 0; i < output.dssdCount(); i++) { //only looking at DSSDs
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
            if(i == 1 || i == 2) { //downstream detectorer
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


        auto impulse_alpha = sqrt(2*Ea*ALPHA_MASS);
        auto impulse_triton = sqrt(2*Et*TRITON->mass);

        hit.index = i;
        hit.lVector_alpha = {impulse_alpha*hit.direction, hit.Ea+ALPHA_MASS};
        hit.lVector_triton = {impulse_triton*hit.direction, hit.Et+TRITON->mass};
        hits.emplace_back(move(hit));
      }
    }
  }


  void doAnalysis() {
      auto mult = hits.size();

      //if(hits.empty()) return;
      if (mult<2) return;

      for(size_t i = 0; i < mult; i++) {
          auto h0 = hits[i];

          auto pl0a = h0.lVector_alpha;
          auto pl0t = h0.lVector_triton;

          auto p0a = TVector3(pl0a.Px(), pl0a.Py(), pl0a.Pz());
          auto p0t = TVector3(pl0t.Px(), pl0t.Py(), pl0t.Pz());

          for(size_t j = i+1; j < mult; j++) {
              //cout << "i, j: " << i << ", " << j << endl;
              auto h1 = hits[j];

              auto pl1a = h1.lVector_alpha;
              auto pl1t = h1.lVector_triton;

              auto p1a = TVector3(pl1a.Px(), pl1a.Py(), pl1a.Pz());
              auto p1t = TVector3(pl1t.Px(), pl1t.Py(), pl1t.Pz());

              auto p_at = p0a + p1t;
              auto p_ta = p0t + p1a;

              //auto pn_at = TVector3(0,0,0) - TVector3(p_at.Px(),p_at.Py(),p_at.Pz());
              //auto pn_ta = TVector3(0,0,0) - TVector3(p_ta.Px(),p_ta.Py(),p_ta.Pz());
              //p_at.Boost(cmBoost);
              //p_ta.Boost(cmBoost);
              auto pn_at = TVector3(0,0,0) - p_at;
              auto pn_ta = TVector3(0,0,0) - p_ta;

              auto En_at = pow(pn_at.Mag(), 2)/(2*NEUTRON_MASS);
              auto En_ta = pow(pn_ta.Mag(), 2)/(2*NEUTRON_MASS);

              v_pnat->add(pn_at);
              v_pnta->add(pn_ta);

              v_Enat->add(En_at);
              v_Enta->add(En_ta);

              v_pos0->add(h0.position);
              v_dir0->add(h0.direction);
              v_theta0->add(h0.theta * TMath::RadToDeg());
              v_pos1->add(h1.position);
              v_dir1->add(h1.direction);
              v_theta1->add(h1.theta * TMath::RadToDeg());

              v_Ea0->add(h0.Ea);
              v_Et0->add(h0.Et);
              v_Edssd0->add(h0.Edssd);
              v_BE0->add(h0.BE);
              v_FE0->add(h0.FE);
              v_Ea1->add(h1.Ea);
              v_Et1->add(h1.Et);
              v_Edssd1->add(h1.Edssd);
              v_BE1->add(h1.BE);
              v_FE1->add(h1.FE);

              v_dE0->add(h0.dE);
              v_dE1->add(h1.dE);
              //v_Ecm->add(Ecm);
              v_ang0->add(h0.angle * TMath::RadToDeg());
              v_ang1->add(h1.angle * TMath::RadToDeg());

              v_i0->add(static_cast<short>(h0.index));
              v_F0->add(h0.fseg);
              v_B0->add(h0.bseg);
              v_FT0->add(h0.TF);
              v_BT0->add(h0.TB);
              v_i1->add(static_cast<short>(h1.index));
              v_F1->add(h1.fseg);
              v_B1->add(h1.bseg);
              v_FT1->add(h1.TF);
              v_BT1->add(h1.TB);

              auto hitAngle = h1.position.Angle(h0.position)*TMath::RadToDeg();

              v_hitAng->add(hitAngle);

              mul++;
          }
    }
  }

  static TLorentzVector constructBeamVector(const Ion& beam,
                                            const Ion& targetIon,
                                            double beamEnergy) {
      TLorentzVector plbeam( TVector3(0,0,sqrt(2*beamEnergy*beam.getMass())), beamEnergy+beam.getMass() );
      TLorentzVector pltarget( TVector3(0,0,0), targetIon.getMass() );
      return plbeam + pltarget;
  }

    void terminate() override {
    AbstractSortedAnalyzer::terminate();
    gDirectory->WriteTObject(t);
  }

  void clear() {
    mul = 0;
    AUSA::clear(
        *v_Et0, *v_Ea0, *v_theta0, *v_Edssd0, *v_Et1, *v_Ea1, *v_theta1, *v_Edssd1,
        *v_i0, *v_FE0, *v_BE0, *v_i1, *v_FE1, *v_BE1, *v_Enat, *v_Enta,
        *v_F0, *v_B0, *v_Ecm0, *v_F1, *v_B1, *v_Ecm1,
        *v_ang0, *v_pos0, *v_dir0, *v_ang1, *v_pos1, *v_dir1,
        *v_dE0, *v_FT0, *v_BT0, *v_dE1, *v_FT1, *v_BT1, *v_hitAng, *v_pnat, *v_pnta
    );
  }

  int NUM;
  TTree *t;
  unique_ptr<DynamicBranchVector<TVector3>> v_dir0, v_pos0, v_dir1, v_pos1;
  unique_ptr<DynamicBranchVector<TVector3>> v_pnat, v_pnta;
  unique_ptr<DynamicBranchVector<double>> v_Edssd0, v_Edssd1, v_Enat, v_Enta;
  unique_ptr<DynamicBranchVector<double>> v_Ea0, v_Et0, v_BE0, v_FE0, v_theta0, v_dE0, v_Ecm0;
  unique_ptr<DynamicBranchVector<double>> v_Ea1, v_Et1, v_BE1, v_FE1, v_theta1, v_dE1, v_Ecm1;
  unique_ptr<DynamicBranchVector<short>> v_i0, v_i1;
  unique_ptr<DynamicBranchVector<short>> v_F0, v_B0, v_F1, v_B1;
  unique_ptr<DynamicBranchVector<double>> v_ang0, v_ang1, v_hitAng;
  unique_ptr<DynamicBranchVector<double>> v_FT0, v_BT0, v_FT1, v_BT1;


  UInt_t mul{}, TPROTONS{};
  vector<Hit> hits;

  unique_ptr<EnergyLossRangeInverter> aSiCalc, tSiCalc;
  vector<unique_ptr<EnergyLossRangeInverter>> aTargetCalcs, tTargetCalcs;
  vector<double> deadlayerF, deadlayerB, deadlayerP;
  Target &target;
  ParticleType* TRITON;
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
