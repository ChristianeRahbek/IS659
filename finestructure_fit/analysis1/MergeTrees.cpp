/**
 *  Created by Christiane Rahbek on 09-03-2023.
**/

#include <TFile.h>
#include <iostream>
#include <TTree.h>
#include <TVector3.h>
#include <TClonesArray.h>
#include "projectutil.h"
#include <fstream>

using namespace std;
using namespace ROOT;

class MergeTTrees {
public:
    MergeTTrees() {
        toBeMerged = new TList();

        filePath = EUtil::getProjectRoot() + "/output/doubleDSSD/"; //change this path if I want to merge other types of analysis files
    }

    void doMerching(){
        addRuns();

        addTreesToList();

        outFile = new TFile((filePath + "mergedRuns.root").c_str(), "RECREATE");
        newTree = TTree::MergeTrees(toBeMerged);

        newTree->SetName("a");
        newTree->Write();
        outFile->Close();

        cout << "Merch was successfull!" << endl;
    }

    TList *toBeMerged;
    TTree *newTree;
    vector<int> runs;
    TFile *outFile;

    string filePath;

private:
    void addTreesToList() { // adds each run to a list
        for(auto& run : runs){
            TFile *f = new TFile((filePath + "Run" + to_string(run) + "mliolio.root").c_str());
            TTree *a = (TTree *) f->Get("a");

            toBeMerged->Add(a);
        }
    }
    void addRuns(){ //reads run number from file
        // read from run file
        string lineInFile;
        ifstream infile(EUtil::getProjectRoot() + "/runsToBeMerged.txt");

        while (getline(infile, lineInFile)) {
            istringstream iss(lineInFile);
            int run;
            if(!(iss >> run)) {
                throw runtime_error("Error parsing runsToBeMerged.txt");
            }

            runs.emplace_back(run);
            cout << "Adding run " << run << endl;
        }
    }
};



int main() {
    auto makingNewTree = new MergeTTrees();

    cout << "DoMerching..." << endl;
    makingNewTree->doMerching();

    return EXIT_SUCCESS;
}