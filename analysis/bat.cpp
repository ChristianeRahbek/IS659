//
// Created by erik on 6/9/22.
// (Beta-)alpha-triton identification
//

#include "projectutil.h"
#include <libconfig.h++>
#include <iostream>
#include <string>
#include <ausa/json/IO.h>
#include <ausa/sort/analyzer/AbstractSortedAnalyzer.h>
#include <ausa/sort/SortedReader.h>
#include <ausa/util/StringUtil.h>
#include <ausa/util/FileUtil.h>
#include <ausa/setup/DoubleSidedSiliconDetector.h>
#include <ausa/setup/SingleSidedSiliconDetector.h>
#include <ausa/setup/SquareDSSD.h>
#include <ausa/setup/PadDetector.h>
#include <ausa/setup/SquareSSSD.h>
#include <ausa/util/DynamicBranchVector.h>
#include <ausa/eloss/Ion.h>
#include <ausa/eloss/Default.h>
#include <ausa/output/OutputConvenience.h>
#include <Math/Vector3D.h>
#include "Hit.h"
#include <TROOT.h>
#include <ctime>
#include <set>
#include <utility>
#include "IS659Detector.h"
#include <ausa/constants/Mass.h>


using namespace std;
using namespace libconfig;
using namespace AUSA;
using namespace ROOT::Math;
using namespace AUSA::Sort;
using namespace AUSA::EnergyLoss;
using namespace IS659Detector_enums;
using namespace AUSA::Constants;

#define TGATE1 1500.0 // ms
#define TGATE2 600.0  // ms

class SinglesAnalysis : public AbstractSortedAnalyzer {
public:
  SinglesAnalysis(Target &target, TFile *output, bool sim) : target(target), sim(sim) {
    output->cd(); // ensure output file is used for mid-analysis dumping and post-analysis saving

    tree = new TTree("a", "a");

    NUM = 0;
    tree->Branch("num", &NUM);
    tree->Branch("mul", &mul);

    v_id = make_unique<DynamicBranchVector<unsigned short>>(*tree, "id", "mul");

    v_dir = make_unique<DynamicBranchVector<TVector3>>(*tree, "dir");
    v_pos = make_unique<DynamicBranchVector<TVector3>>(*tree, "pos");

    v_theta = make_unique<DynamicBranchVector<double>>(*tree, "theta", "mul");
    v_phi = make_unique<DynamicBranchVector<double>>(*tree, "phi", "mul");
    v_angle = make_unique<DynamicBranchVector<double>>(*tree, "angle", "mul"); // angle of incidence w.r.t. detector surface

    v_Edep = make_unique<DynamicBranchVector<double>>(*tree, "Edep", "mul");

    // double-sided detector info
    v_FI = make_unique<DynamicBranchVector<unsigned short>>(*tree, "FI", "mul");
    v_BI = make_unique<DynamicBranchVector<unsigned short>>(*tree, "BI", "mul");
    v_FE = make_unique<DynamicBranchVector<double>>(*tree, "FE", "mul");
    v_BE = make_unique<DynamicBranchVector<double>>(*tree, "BE", "mul");
    v_FT = make_unique<DynamicBranchVector<double>>(*tree, "FT", "mul");
    v_BT = make_unique<DynamicBranchVector<double>>(*tree, "BT", "mul");

    v_Ea = make_unique<DynamicBranchVector<double>>(*tree, "Ea", "mul");
    v_Et = make_unique<DynamicBranchVector<double>>(*tree, "Et", "mul");

    tree->Branch("TPROTONS", &TPROTONS);
    tree->Branch("TSTAMP", &TSTAMP);

    tSiCalc = defaultRangeInverter("t", "Silicon"); // tAlCalc = defaultRangeInverter("t", "Aluminum");
    aSiCalc = defaultRangeInverter("a", "Silicon"); // aAlCalc = defaultRangeInverter("a", "Aluminum");
    for (auto &layer : target.getLayers()) {
      tTargetCalcs.push_back(defaultRangeInverter(Ion::predefined("t"), layer.getMaterial()));
      aTargetCalcs.push_back(defaultRangeInverter(Ion::predefined("a"), layer.getMaterial()));
    }
  }

  void setup(const SortedSetupOutput &output) override {
    AbstractSortedAnalyzer::setup(output);

    for (const auto &entry : detector_map) {
      auto id = entry.first;
      auto &det = entry.second;
      if (det->getType() == dsd) {
        deadlayerF.emplace(id, getFrontDeadLayer(output.getDssdOutput(det->getName()).detector()));
      }
    }

    // TODO: investigate possibility of enhancing numerical precision by lowering time units in xia4ids
    // unit of times are decided by the configuration file used by xia4ids - "run_unit" and "ref_unit" respectively; currently the units are 1 second and 1 millisecond respectively
    tstamp = output.getScalerOutput("TSTAMP"); // global time stamp with some reference point -- a clock which was incremented with time, constantly, during the entire experiment
    tprotons = output.getScalerOutput("TPROTONS"); // time after proton pulse hit production target
  }


  void analyze() override {
    clear();
    TSTAMP = tstamp.getValue();
    TPROTONS = tprotons.getValue();
    if (TPROTONS/1e5 < TGATE1 /*ms*/ ) { // TPROTONS:E plot indicates that no gammas from 8He are present in the chamber anymore
      findHits();
      doAnalysis();
      if (mul > 1) { tree->Fill(); }
    }
    NUM++;
  }

  void findHits() {
    for (const auto &entry : detector_map) {
      auto det = entry.second;
      unsigned short id = det->getId();
      auto type = det->getType();
      if (type == dsd) {
        findDssdHit(det);
      }
//      else if (type == IS659Type::SquareSSSD || type == IS659Type::Pad) {
//        findSsdHit(det);
        /*
         * Ain't gonna happen
         * } else { // => type == Clover
         *   findCloverHit(det);
         */
//      }
    }
  }

  void findDssdHit(IS659Detector* detector) {
    auto& o = output.getDssdOutput(detector->getName());
    auto& d = o.detector();
    auto m = AUSA::mul(o);

    for (int i = 0; i < m; i++) {
      Hit hit;
      hit.id = detector->getId();

      hit.Edep = energy(o, i);

      auto FI = fSeg(o, i);
      auto BI = bSeg(o, i);
      hit.FI = short(FI);
      hit.BI = short(BI);
      hit.FE = fEnergy(o, i);
      hit.BE = bEnergy(o, i);
      hit.FT = fTime(o, i);
      hit.BT = bTime(o, i);

      TVector3 pos = d.getUniformPixelPosition(FI, BI);
      hit.pos = pos;
      TVector3 dir = (hit.pos - target.getCenter()).Unit();
      hit.dir = dir;

      hit.theta = dir.Theta();
      hit.phi = dir.Phi();
      hit.angle = dir.Angle(-d.getNormal());

      hits.emplace_back(move(hit));
    }
  }

  void findSsdHit(IS659Detector* detector) {
    auto& o = output.getSingleOutput(detector->getName());
    auto& d = o.detector();
    auto m = AUSA::mul(o);

    for (int i = 0; i < m; i++) {
      Hit hit;
      hit.id = detector->getId();

      hit.Edep = o.energy(i);

      auto FI = o.segment(i);
      hit.FI = short(FI);
      hit.FE = hit.Edep;
      hit.FT = o.time(i);

      TVector3 pos = d.getPosition(FI);
      hit.pos = pos;
      TVector3 dir = (hit.pos - target.getCenter()).Unit();
      hit.dir = dir;

      hit.theta = dir.Theta();
      hit.phi = dir.Phi();
      hit.angle = dir.Angle(-d.getNormal());

      hits.emplace_back(move(hit));
    }
  }

  /*
   * TODO: Add key to matcher file for Clover detectors, so tolerance and low_cut can be set properly
   * Checking for hit.Edep > 0.0 for now.
   */
  void findCloverHit(IS659Detector* detector) {
    auto& o = output.getSingleOutput(detector->getName());
    auto& d = o.detector();
    auto m = AUSA::mul(o);

    Hit hit;
    hit.id = detector->getId();

    hit.Edep = 0.0;
    for (int i = 0; i < m; i++) {
      hit.Edep += o.energy(i); // addback: sum gamma energies recorded in all segments
    }
    if (hit.Edep == 0.0) return;
    hits.emplace_back(move(hit));
  }

  void doAnalysis() {
    if (hits.empty()) return;

    // only events with hits in beta-free + punchthrough-free region
    double Edep;
    for (int i = 0; i < hits.size(); /*conditional increment*/ ) {
      Edep = hits[i].Edep;
//      if (front_thresholds.at(hits[i].id).first <= Edep && Edep <= front_thresholds.at(hits[i].id).second) {
      if (true) {
        i++; // only increase i if hits[i] is not removed
        continue;
      }
      hits.erase(hits.begin() + i);
    }

    // only events with hits in different front detectors
    set<unsigned short> ids;
    for (auto &hit : hits) {
      ids.insert(hit.id);
    }

    mul = ids.size();
    if (mul < 2) return;

    for (auto &hit : hits) {
      treatDssdHit(&hit);
      auto type = detector_map.at(hit.id)->getType();
      if (type == dsd) {
        addDssdHit(&hit);
//      } else {
//        addSsdHit(&hit);
      }
    }

  }

  bool treatDssdHit(Hit* hit) {
//    double angle = detector_map.at(hit->id)->getType() == IS507Type::SquareDSSD ? hit->angle : 0.0;
    double angle = hit->angle;
    double t = deadlayerF.at(hit->id) / abs(cos(angle));
    double Ea = hit->Edep;
    double Et = hit->Edep;
    Ea += aSiCalc->getTotalEnergyCorrection(Ea, t);
    Et += tSiCalc->getTotalEnergyCorrection(Et, t);
    auto &from = hit->pos;
    for (auto &intersection : target.getIntersections(from, target.getCenter())) {
      auto &aCalc = aTargetCalcs[intersection.index];
      auto &tCalc = tTargetCalcs[intersection.index];
      Ea += aCalc->getTotalEnergyCorrection(Ea, intersection.transversed);
      Et += tCalc->getTotalEnergyCorrection(Et, intersection.transversed);
    }

    hit->Ea = Ea;
    hit->Et = Et;
    return true;
  }

//  bool treatBackHit(Hit* hit) {
//    double angle = detector_map.at(hit->id)->getType() == IS659Type::SquareDSSD ? hit->angle : 0.0;
//    double t = bdeadlayer.at(hit->id) / abs(cos(angle));
//    double E = hit->Edep;
//    E += tSiCalc->getTotalEnergyCorrection(E, t);
//    auto &from = hit->pos;
//    for (auto &intersection : target.getIntersections(from, target.getCenter())) {
//      auto &calc = tTargetCalcs[intersection.index];
//      E += calc->getTotalEnergyCorrection(E, intersection.transversed);
//    }
//
//    if (E < time_after_production_thresholds.at(hit->id)) {
//      if (TPROTONS/1e5 > TGATE2 /*ms*/ ) {
//        return false;
//      }
//    }
//    hit->E = E;
//    hit->Ecm = MFRAC*E;
//    return true;
//  }
//
//  void treatOutsideHit(Hit* hit) {
//    hit->Eg = hit->Edep;
//  }
//
//  bool treatTelescopeHit(Hit* front_hit, Hit* back_hit) {
//    if (front_hit->Edep == 0.0 || back_hit->Edep == 0.0) return false;
//
//    auto fid = front_hit->id;
//    auto bid = back_hit->id;
//    double angle = detector_map.at(fid)->getType() == IS659Type::SquareDSSD ? front_hit->angle : back_hit->angle;
//    TVector3 pos = detector_map.at(fid)->getType() == IS659Type::SquareDSSD ? front_hit->pos : back_hit->pos;
//
//    auto tF = deadlayerF.at(fid) / abs(cos(angle));
//    auto tB = deadlayerB.at(fid) / abs(cos(angle));
//    auto fcB = fcontactB.at(fid) / abs(cos(angle));
//    auto bcF = bcontactF.at(bid) / abs(cos(angle));
//    auto bt = bdeadlayer.at(bid) / abs(cos(angle));
//    double E = 0.0;
//    E += back_hit->Edep;
//    E += tSiCalc->getTotalEnergyCorrection(E, bt);
//    E += tAlCalc->getTotalEnergyCorrection(E, bcF);
//    E += tAlCalc->getTotalEnergyCorrection(E, fcB);
//    E += tSiCalc->getTotalEnergyCorrection(E, tB);
//    E += front_hit->Edep;
//    E += tSiCalc->getTotalEnergyCorrection(E, tF);
//    auto &from = pos;
//    for (auto &intersection : target.getIntersections(from, target.getCenter())) {
//      auto &calc = tTargetCalcs[intersection.index];
//      E += calc->getTotalEnergyCorrection(E, intersection.transversed);
//    }
//
//    if (E < time_after_production_thresholds.at(front_hit->id)) {
//      if (TPROTONS/1e5 > TGATE2 /*ms*/ ) {
//        return false;
//      }
//    }
//    front_hit->E = E;
//    front_hit->Ecm = MFRAC*E;
//    return true;
//  }

  void addDssdHit(Hit* hit) {
    v_id->add(hit->id);

    v_dir->add(hit->dir);
    v_pos->add(hit->pos);

    v_theta->add(hit->theta);
    v_phi->add(hit->phi);
    v_angle->add(hit->angle);

    v_Edep->add(hit->Edep);

    v_FI->add(hit->FI);
    v_BI->add(hit->BI);
    v_FE->add(hit->FE);
    v_BE->add(hit->BE);
    v_FT->add(hit->FT);
    v_BT->add(hit->BT);

    v_Ea->add(hit->Ea);
    v_Et->add(hit->Et);

    mul++;
  }

//  void addSsdHit(Hit* hit) {
//    v_id->add(hit->id);
//
//    v_dir->add(NAN_TVECTOR3);
//    v_pos->add(NAN_TVECTOR3);
//
//    v_theta->add(NAN);
//    v_phi->add(NAN);
//    v_angle->add(NAN);
//
//    v_Edep->add(hit->Edep);
//    v_fEdep->add(NAN);
//    v_bEdep->add(NAN);
//
//    v_FI->add(hit->FI);
//    v_BI->add(NAN_UINT);
//    v_FE->add(hit->FE);
//    v_BE->add(NAN);
//    v_FT->add(hit->FT);
//    v_BT->add(NAN);
//
//    v_E->add(hit->E);
//    v_Ecm->add(hit->Ecm);
//    v_Eg->add(NAN);
//
//    v_Ea->add(hit->Ea);
//    v_Eacm->add(hit->Eacm);
//
//    mul++;
//  }
//
//  void addCloverHit(Hit* hit) {
//    v_id->add(hit->id);
//
//    v_dir->add(NAN_TVECTOR3);
//    v_pos->add(NAN_TVECTOR3);
//
//    v_theta->add(NAN);
//    v_phi->add(NAN);
//    v_angle->add(NAN);
//
//    v_Edep->add(hit->Edep);
//    v_fEdep->add(NAN);
//    v_bEdep->add(NAN);
//
//    v_FI->add(NAN_UINT);
//    v_BI->add(NAN_UINT);
//    v_FE->add(NAN);
//    v_BE->add(NAN);
//    v_FT->add(NAN);
//    v_BT->add(NAN);
//
//    v_E->add(NAN);
//    v_Ecm->add(NAN);
//    v_Eg->add(hit->Eg);
//
//    mul++;
//  }
//
//  void addDssdSsdHits(Hit* front_hit, Hit* back_hit) {
//    v_id->add(front_hit->id);
//
//    v_dir->add(front_hit->dir);
//    v_pos->add(front_hit->pos);
//
//    v_theta->add(front_hit->theta);
//    v_phi->add(front_hit->phi);
//    v_angle->add(front_hit->angle);
//
//    v_Edep->add(NAN);
//    v_fEdep->add(front_hit->Edep);
//    v_bEdep->add(back_hit->Edep);
//
//    v_FI->add(front_hit->FI);
//    v_BI->add(front_hit->BI);
//    v_FE->add(front_hit->FE);
//    v_BE->add(front_hit->BE);
//    v_FT->add(front_hit->FT);
//    v_BT->add(front_hit->BT);
//
//    v_E->add(front_hit->E);
//    v_Ecm->add(front_hit->Ecm);
//    v_Eg->add(NAN);
//
//    mul++;
//  }
//
//  void addSsdDssdHits(Hit* front_hit, Hit* back_hit) {
//    v_id->add(front_hit->id);
//
//    v_dir->add(back_hit->dir);
//    v_pos->add(back_hit->pos);
//
//    v_theta->add(back_hit->theta);
//    v_phi->add(back_hit->phi);
//    v_angle->add(back_hit->angle);
//
//    v_Edep->add(NAN);
//    v_fEdep->add(front_hit->Edep);
//    v_bEdep->add(back_hit->Edep);
//
//    v_FI->add(back_hit->FI);
//    v_BI->add(back_hit->BI);
//    v_FE->add(back_hit->FE);
//    v_BE->add(back_hit->BE);
//    v_FT->add(back_hit->FT);
//    v_BT->add(back_hit->BT);
//
//    v_E->add(front_hit->E);
//    v_Ecm->add(front_hit->Ecm);
//    v_Eg->add(NAN);
//
//    mul++;
//  }

  void terminate() override {
    AbstractSortedAnalyzer::terminate();
    gDirectory->WriteTObject(tree);
  }

  void clear() {
    mul = 0;
    AUSA::clear(
        *v_id,
        *v_dir, *v_pos,
        *v_theta, *v_phi, *v_angle,
        *v_Edep,
        *v_FI, *v_BI, *v_FE, *v_BE, *v_FT, *v_BT,
        *v_Ea, *v_Et
    );
    hits.clear();
  }

  /*
   * TODO: Justify the detector map with some images/schematics and some text here.
   */
  const map<unsigned short, IS659Detector*> detector_map = {
      {0,  new IS659Detector(0, "U1", dsd, "P1")},
      {1,  new IS659Detector(1, "U2", dsd, "P2")},
      {2,  new IS659Detector(2, "U3", dsd, "P3")},
      {3,  new IS659Detector(3, "U4", dsd, "P4")},
  };

  TTree *tree;
  int NUM;
  UInt_t mul{}, TSTAMP{}, TPROTONS{};
  unique_ptr<DynamicBranchVector<unsigned short>> v_id;
  unique_ptr<DynamicBranchVector<TVector3>> v_dir, v_pos;
  unique_ptr<DynamicBranchVector<double>> v_theta, v_phi, v_angle;
  unique_ptr<DynamicBranchVector<double>> v_Edep;
  unique_ptr<DynamicBranchVector<unsigned short>> v_FI, v_BI;
  unique_ptr<DynamicBranchVector<double>> v_FE, v_BE;
  unique_ptr<DynamicBranchVector<double>> v_FT, v_BT;
  unique_ptr<DynamicBranchVector<double>> v_Ea, v_Et;

  vector<Hit> hits;

  unique_ptr<EnergyLossRangeInverter> tSiCalc, tAlCalc, aSiCalc, aAlCalc;
  vector<unique_ptr<EnergyLossRangeInverter>> tTargetCalcs, aTargetCalcs;
  map<unsigned short, double> deadlayerF;
  Target &target;

  SortedSignal tstamp, tprotons;

  bool sim;
};

string setup_path, target_path, input_path, output_dir;
bool sim;
void prepareFileIO(const string& configfile) {
  Config cfg;
  cfg.readFile(configfile.c_str());

  if (cfg.lookup("paths_relative_to_project_root")) {
    setup_path  = EUtil::getProjectRoot() + "/" + cfg.lookup("setup_file").c_str();
    target_path = EUtil::getProjectRoot() + "/" + cfg.lookup("target_file").c_str();
    input_path  = EUtil::getProjectRoot() + "/" + cfg.lookup("data_input").c_str();
    output_dir  = EUtil::getProjectRoot() + "/" + cfg.lookup("data_output_dir").c_str();
  } else {
    setup_path  = cfg.lookup("setup_file").c_str();
    target_path = cfg.lookup("target_file").c_str();
    input_path  = cfg.lookup("data_input").c_str();
    output_dir  = cfg.lookup("data_output_dir").c_str();
  }

  sim = (EUtil::getParentPath(input_path).find("_sim") != string::npos);
  if (sim) {
    setup_path = EUtil::appendToStem(setup_path, "_sim");
    output_dir += "_sim";
  }


  if (cfg.lookup("verbose")) {
    cout << "------------------------ IO configuration ------------------------" << endl;
    cout << "Setup:  " << setup_path  << endl;
    cout << "Target: " << target_path << endl;
    cout << "Input:  " << input_path  << endl;
    cout << "Output: " << output_dir  << endl;
    cout << "Sim:    " << sim         << endl;
    cout << "------------------------------------------------------------------" << endl << endl;
  }
}

// ( seq 574 576 ; seq 624 630 ; seq 641 648 ; seq 671 677 ) | parallel -u ./bpa
int main(int argc, char *argv[]) {
  prepareFileIO(EUtil::getProjectRoot() + "/analysis/bat.cfg");
  auto setup = JSON::readSetupFromJSON(setup_path);
  auto target = JSON::readTargetFromJSON(target_path);

  vector<string> input;
  int run;

  for (int i = 1; i < argc; i++) {
    run = stoi(argv[i]);
    findFilesMatchingWildcard(Form(input_path.c_str(), run), input);
  }

  system(("mkdir -p " + output_dir).c_str());
  for (auto &in : input) {
    clock_t start = clock();

    SortedReader reader{*setup};
    reader.add(in);
    reader.setVerbose(true);

    string stem = EUtil::getStem(in);
    TString outfile = (output_dir + "/" + stem + "lio.root").c_str();

    cout << "Reading from: " << in << endl;
    cout << "Printing to:  " << outfile << endl;

    TFile output(outfile, "RECREATE");
    auto analysis = make_shared<SinglesAnalysis>(target, &output, sim);
    reader.attach(analysis);
    reader.run();

    clock_t stop = clock();
    double elapsed = (double) (stop - start) / CLOCKS_PER_SEC;
    printf("\nFile: %s.\tTime elapsed: %.5f\n", (stem + "lio.root").c_str(), elapsed);
  }

  return 0;
}