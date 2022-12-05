#include <iostream>
#include <string>
#include <ausa/json/IO.h>
#include <ausa/sort/analyzer/AbstractSortedAnalyzer.h>
#include <ausa/sort/SortedReader.h>
#include <ausa/util/StringUtil.h>
#include <ausa/util/FileUtil.h>
#include <ausa/setup/DoubleSidedSiliconDetector.h>
#include <ausa/setup/SingleSidedSiliconDetector.h>
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

int main(int argc, char *argv[]) {
    string setup_path, input_path, output_dir;
    Config cfg;
    auto configfile = ANALYSIS::getProjectRoot() + "/Analysis.cfg";
    cfg.readFile(configfile.c_str());

    if (cfg.lookup("paths_relative_to_project_root")) {
        setup_path = ANALYSIS::getProjectRoot() + "/" + cfg.lookup("setup_file").c_str();
        input_path = ANALYSIS::getProjectRoot() + "/" + cfg.lookup("data_input").c_str();
        output_dir = ANALYSIS::getProjectRoot() + "/" + cfg.lookup("data_output_dir").c_str();
    } else {
        setup_path = cfg.lookup("setup_file").c_str();
        input_path = cfg.lookup("data_input").c_str();
        output_dir = cfg.lookup("data_output_dir").c_str();
    }

    TString outfile = (output_dir + "/" + "212lio.root").c_str();
    auto setup = JSON::readSetupFromJSON(setup_path);
    TFile output1(outfile, "RECREATE");
    const SortedSetupOutput output = &output1;

    for (size_t i = 0; i < output.dssdCount(); ++i) {
        auto dl = getFrontDeadLayer(output.getDssdOutput(i).detector());
        auto dlB = getBackDeadLayer(output.getDssdOutput(i).detector());
        //auto dlP = getFrontDeadLayer(output.getSingleOutput(i).detector());
        deadlayerF.push_back(dl);
        deadlayerB.push_back(dlB);
        //deadlayerP.push_back(dlP);
    }

    return 0;
}
