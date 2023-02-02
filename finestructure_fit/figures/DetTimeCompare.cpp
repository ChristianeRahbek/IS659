//
// Made by Christiane 23/01-2023
//

#include "THStack.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"

void DetTimeCompare() {
    /* PRE-SETS */
    std::string runNo = "212";
    int BackIndex = 2; //only relevant if we do set back or front index to false
    int FrontIndex = 2; //only relevant if we do set back or front index to false
    bool frontIndex = true; //will sort after front index set above if true
    bool backIndex = true; //will sort after back index set above if true

    /* ********************************************************** */

    /* Preparing for for-loop */
    int fi_max, bi_max, fi_min, bi_min;
    if(frontIndex && backIndex) {
        fi_max = 17;
        bi_max = 17;
        fi_min = 1;
        bi_min = 1;
    } else if (frontIndex) {
        fi_max = 17;
        bi_max = BackIndex + 1;
        bi_min = BackIndex;
        fi_min = 1;
    } else if (backIndex) {
        fi_max = FrontIndex + 1;
        bi_max = 17;
        bi_min = 1;
        fi_min = FrontIndex;
    } else {
        fi_max = FrontIndex + 1;
        bi_max = BackIndex + 1;
        bi_min = BackIndex;
        fi_min = FrontIndex;
    }

    for(int fi = fi_min; fi < fi_max; fi++) { //looping over all possible front indexes
        for (int bi = bi_min; bi < bi_max; bi++) {

            std::string Findex = std::to_string(fi);
            std::string Bindex = std::to_string(bi);
            std::string selExtra;
            if (frontIndex && backIndex) {
                selExtra = " && FI == " + Findex + " && BI == " + Bindex;
            } else if (backIndex) {
                selExtra = " && BI == " + Bindex;
            } else if (frontIndex) {
                selExtra = " && FI == " + Findex;
            } else selExtra = ""; //to look at all data, have both front and back index be 0

            auto canv = new TCanvas("", "", 700, 700);
            canv->cd()->SetLogy();

            auto *f = new TFile(("/mnt/d/IS659/finestructure_fit/analysis1/output/Run" + runNo + "mlio.root").c_str());
            auto *tr = (TTree *) f->Get("a");

            auto h1st = new TH1F("h1st", "h1st", 500, -10, 50000);
            tr->Draw("abs(FT[0]-BT[0]) >> h1st", ("id == 0" + selExtra).c_str());
            h1st->SetLineColor(kRed);

            auto h2st = new TH1F("h2st", "h2st", 500, -10, 50000);
            tr->Draw("abs(FT[0]-BT[0]) >> h2st", ("id == 1" + selExtra).c_str());
            h2st->SetLineColor(kBlue);

            auto h3st = new TH1F("h3st", "h3st", 500, -10, 50000);
            tr->Draw("abs(FT[0]-BT[0]) >> h3st", ("id == 2" + selExtra).c_str());
            h3st->SetLineColor(kGreen);

            auto h4st = new TH1F("h4st", "h4st", 500, -10, 50000);
            tr->Draw("abs(FT[0]-BT[0]) >> h4st", ("id == 3" + selExtra).c_str());
            h4st->SetLineColor(kOrange);

            //Plotting every histogram together
            auto hs = new THStack("hs", ("Detector front/back time comparison, run " + runNo).c_str());
            hs->Add(h1st, "hist");
            hs->Add(h2st, "hist");
            hs->Add(h3st, "hist");
            hs->Add(h4st, "hist");

            auto l = new TLegend(0.7, 0.7, 0.9, 0.9);

            l->AddEntry(h1st, "Id = 0", "L");
            l->AddEntry(h2st, "Id = 1", "L");
            l->AddEntry(h3st, "Id = 2", "L");
            l->AddEntry(h4st, "Id = 3", "L");

            hs->Draw("nostack");
            l->Draw();

            if (selExtra == "") {
                canv->SaveAs(("DetTimeCompare/Run" + runNo + "/DetTimeCompareRun" + runNo + ".png").c_str());
            } else canv->SaveAs(("DetTimeCompare/Run" + runNo + "/DetTimeCompareRun" + runNo + "FI" + Findex + "BI" + Bindex + ".png").c_str());

        }
    }
}

