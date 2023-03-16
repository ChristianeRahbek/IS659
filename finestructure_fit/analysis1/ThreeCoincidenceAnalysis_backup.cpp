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

        beamEnergy = 30; //keV
        beamVector = constructBeamVector(Ion(2, 8), Ion(6, 12), beamEnergy);
        cmBoost = -beamVector.BoostVector();

        TRITON = new ParticleType("H3");

        t = new TTree("a", "a");
        t->Branch("mul", &mul);
        t->Branch("num", &NUM);

        v_dir0 = make_unique<DynamicBranchVector<TVector3>>(*t, "dir0");
        v_pos0 = make_unique<DynamicBranchVector<TVector3>>(*t, "pos0");
        v_dir1 = make_unique<DynamicBranchVector<TVector3>>(*t, "dir1");
        v_pos1 = make_unique<DynamicBranchVector<TVector3>>(*t, "pos1");
        v_dir2 = make_unique<DynamicBranchVector<TVector3>>(*t, "dir2");
        v_pos2 = make_unique<DynamicBranchVector<TVector3>>(*t, "pos2");

        v_theta0 = make_unique<DynamicBranchVector<double>>(*t, "theta0", "mul");
        v_ang0 = make_unique<DynamicBranchVector<double>>(*t, "angle0", "mul");
        v_theta1 = make_unique<DynamicBranchVector<double>>(*t, "theta1", "mul");
        v_ang1 = make_unique<DynamicBranchVector<double>>(*t, "angle1", "mul");
        v_theta2 = make_unique<DynamicBranchVector<double>>(*t, "theta2", "mul");
        v_ang2 = make_unique<DynamicBranchVector<double>>(*t, "angle2", "mul");

        v_hitAng = make_unique<DynamicBranchVector<double>>(*t, "hitAng", "mul");

        v_Edep0 = make_unique<DynamicBranchVector<double>>(*t, "Edep0", "mul");
        v_Ea0 = make_unique<DynamicBranchVector<double>>(*t, "Ea0", "mul");
        v_BE0 = make_unique<DynamicBranchVector<double>>(*t, "BE0", "mul");
        v_FE0 = make_unique<DynamicBranchVector<double>>(*t, "FE0", "mul");
        v_Et0 = make_unique<DynamicBranchVector<double>>(*t, "Et0", "mul");
        v_Edep1 = make_unique<DynamicBranchVector<double>>(*t, "Edep1", "mul");
        v_Ea1 = make_unique<DynamicBranchVector<double>>(*t, "Ea1", "mul");
        v_BE1 = make_unique<DynamicBranchVector<double>>(*t, "BE1", "mul");
        v_FE1 = make_unique<DynamicBranchVector<double>>(*t, "FE1", "mul");
        v_Et1 = make_unique<DynamicBranchVector<double>>(*t, "Et1", "mul");
        v_Edep2 = make_unique<DynamicBranchVector<double>>(*t, "Edep2", "mul");
        v_Ea2 = make_unique<DynamicBranchVector<double>>(*t, "Ea2", "mul");
        v_BE2 = make_unique<DynamicBranchVector<double>>(*t, "BE2", "mul");
        v_FE2 = make_unique<DynamicBranchVector<double>>(*t, "FE2", "mul");
        v_Et2 = make_unique<DynamicBranchVector<double>>(*t, "Et2", "mul");

        v_FT0 = make_unique<DynamicBranchVector<double>>(*t, "FT0", "mul");
        v_BT0 = make_unique<DynamicBranchVector<double>>(*t, "BT0", "mul");
        v_FT1 = make_unique<DynamicBranchVector<double>>(*t, "FT1", "mul");
        v_BT1 = make_unique<DynamicBranchVector<double>>(*t, "BT1", "mul");
        v_FT2 = make_unique<DynamicBranchVector<double>>(*t, "FT2", "mul");
        v_BT2 = make_unique<DynamicBranchVector<double>>(*t, "BT2", "mul");

        v_dE0 = make_unique<DynamicBranchVector<double>>(*t, "dE0", "mul");
        v_Ecm0 = make_unique<DynamicBranchVector<double>>(*t, "Ecm0", "mul");
        v_dE1 = make_unique<DynamicBranchVector<double>>(*t, "dE1", "mul");
        v_Ecm1 = make_unique<DynamicBranchVector<double>>(*t, "Ecm1", "mul");
        v_dE2 = make_unique<DynamicBranchVector<double>>(*t, "dE2", "mul");
        v_Ecm2 = make_unique<DynamicBranchVector<double>>(*t, "Ecm2", "mul");

        v_i0 = make_unique<DynamicBranchVector<short>>(*t, "id0", "mul");
        v_i1 = make_unique<DynamicBranchVector<short>>(*t, "id1", "mul");
        v_i2 = make_unique<DynamicBranchVector<short>>(*t, "id2", "mul");

        v_F0 = make_unique<DynamicBranchVector<short>>(*t, "FI0", "mul");
        v_B0 = make_unique<DynamicBranchVector<short>>(*t, "BI0", "mul");
        v_F1 = make_unique<DynamicBranchVector<short>>(*t, "FI1", "mul");
        v_B1 = make_unique<DynamicBranchVector<short>>(*t, "BI1", "mul");
        v_F2 = make_unique<DynamicBranchVector<short>>(*t, "FI2", "mul");
        v_B2 = make_unique<DynamicBranchVector<short>>(*t, "BI2", "mul");


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
        for (const auto &entry: detector_map) {
            auto det = entry.second;
            auto detType = det->getDetType();
            if (detType == IS659Type::SquareDSSD) findDSSDHit(det);
            else if (detType == IS659Type::Plastic) findPlasticHit(det);
            else throw std::invalid_argument("An invalid detector has been found");
        }
    }

    void findDSSDHit(IS659Detector *detector) {
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

    void findPlasticHit(IS659Detector *detector) {
        auto name = detector->getName();
        auto &o = output.getSingleOutput(name);
        auto &d = o.detector();
        auto m = AUSA::mul(o);

        for (UInt_t j = 0; j < m; j++) {
            Hit hit;

            auto id = detector->getId();

            auto Edep = o.energy(j);

            hit.fseg = short(o.segment(j)); //there are two front IDs
            hit.bseg = short(-1);

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
        cout << "multiplicity is " << mult << endl;

        //if(hits.empty()) return;
        if (mult < 2) return;
        else if(mult == 2) {
            for(size_t i = 0; i < mult; i++) {
                auto h0 = hits[i];

                for(size_t j = i + 1; j < mult; j++){
                    auto h1 = hits[j];

                    if(h0.index > 3 || h1.index > 3) break; //only DSSD coincidences

                    v_pos0->add(h0.position);
                    v_dir0->add(h0.direction);
                    v_theta0->add(h0.theta * TMath::RadToDeg());
                    v_pos1->add(h1.position);
                    v_dir1->add(h1.direction);
                    v_theta1->add(h1.theta * TMath::RadToDeg());
                    v_pos2->add(NAN_TVECTOR3);
                    v_dir2->add(NAN_TVECTOR3);
                    v_theta2->add(NAN);

                    v_Ea0->add(h0.Ea);
                    v_Et0->add(h0.Et);
                    v_Edep0->add(h0.Edep);
                    v_BE0->add(h0.BE);
                    v_FE0->add(h0.FE);
                    v_Ea1->add(h1.Ea);
                    v_Et1->add(h1.Et);
                    v_Edep1->add(h1.Edep);
                    v_BE1->add(h1.BE);
                    v_FE1->add(h1.FE);
                    v_Ea2->add(NAN);
                    v_Et2->add(NAN);
                    v_Edep2->add(NAN);
                    v_BE2->add(NAN);
                    v_FE2->add(NAN);

                    v_dE0->add(h0.dE);
                    v_dE1->add(h1.dE);
                    v_dE2->add(NAN);
                    v_ang0->add(h0.angle * TMath::RadToDeg());
                    v_ang1->add(h1.angle * TMath::RadToDeg());
                    v_ang2->add(NAN);

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
                    v_i2->add(static_cast<short>(-1));
                    v_F2->add(static_cast<short>(-1));
                    v_B2->add(static_cast<short>(-1));
                    v_FT2->add(NAN);
                    v_BT2->add(NAN);

                    double hitAngle = h1.position.Angle(h0.position) * TMath::RadToDeg();

                    v_hitAng->add(hitAngle);

                    mul++;
                }
            }
        }
        else {
            for(size_t i = 0; i < mult; i++) {
                auto h_i = hits[i];

                for(size_t j = i + 1; j < mult; j++) {
                    auto h_j = hits[j];

                    for(size_t k = j + 1; k < mult; k++) {
                        auto h_k = hits[j];
                        auto id_i = h_i.index;
                        auto id_j = h_j.index;
                        auto id_k = h_k.index;

                        if(!((id_i < 4 && (id_j < 4 || id_k < 4)) || (id_j < 4 && (id_i < 4 || id_k < 4)) ||
                          (id_k < 4 && (id_j < 4 || id_i < 4)))) continue;

                        //FROM HERE: HOW To sort them? Should I remove all coincidences with 3 DSSD matches? Then It would be easier...
                        //Lav sortering, som siger, at vi kun vil have Plastics som tilhører de fundne DSSD'er.

                        //I now choose to remove all coincidences with 3 DSSD matches
                        if(id_i < 4 && id_j < 4 && id_k < 4) continue;

                        Hit h0, h1, h2;
                        if(id_i > 3) {
                            h0 = h_j;
                            h1 = h_k;
                            h2 = h_i;
                        }
                        else if(id_j > 3) {
                            h0 = h_i;
                            h1 = h_k;
                            h2 = h_j;
                        }
                        else if(id_k > 3) {
                            h0 = h_i;
                            h1 = h_j;
                            h2 = h_k;
                        }

                        //we only want coincidences from same plastics as DSSD hits.
                        //be aware that we can throw away some perfectly fine coincidences from this...
                        auto id0 = h0.index; auto id1 = h1.index; auto id2 = h2.index;
                        if(!(((id0==0 || id1==0) && id2 == 4) || ((id0==1 || id1==1) && id2 == 5) ||
                           ((id0==2 || id1==2) && id2 == 6) || ((id0==3 || id1==3) && id2 == 7))) continue;

                        v_pos0->add(h0.position);
                        v_dir0->add(h0.direction);
                        v_theta0->add(h0.theta * TMath::RadToDeg());
                        v_pos1->add(h1.position);
                        v_dir1->add(h1.direction);
                        v_theta1->add(h1.theta * TMath::RadToDeg());
                        v_pos2->add(h2.position);
                        v_dir2->add(h2.direction);
                        v_theta2->add(h2.theta * TMath::RadToDeg());

                        v_Ea0->add(h0.Ea);
                        v_Et0->add(h0.Et);
                        v_Edep0->add(h0.Edep);
                        v_BE0->add(h0.BE);
                        v_FE0->add(h0.FE);
                        v_Ea1->add(h1.Ea);
                        v_Et1->add(h1.Et);
                        v_Edep1->add(h1.Edep);
                        v_BE1->add(h1.BE);
                        v_FE1->add(h1.FE);
                        v_Ea2->add(h2.Ea);
                        v_Et2->add(h2.Et);
                        v_Edep2->add(h2.Edep);
                        v_BE2->add(h2.BE);
                        v_FE2->add(h2.FE);

                        v_dE0->add(h0.dE);
                        v_dE1->add(h1.dE);
                        v_dE2->add(h2.dE);
                        v_ang0->add(h0.angle * TMath::RadToDeg());
                        v_ang1->add(h1.angle * TMath::RadToDeg());
                        v_ang2->add(h2.angle * TMath::RadToDeg());

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
                        v_i2->add(static_cast<short>(h2.index));
                        v_F2->add(h2.fseg);
                        v_B2->add(h2.bseg);
                        v_FT2->add(h2.TF);
                        v_BT2->add(h2.TB);

                        double hitAngle;
                        if (h0.index < 4 && h1.index < 4) {
                            hitAngle = h1.position.Angle(h0.position) * TMath::RadToDeg();
                        }
                        else if (h0.index < 4 && h2.index < 4) {
                            hitAngle = h2.position.Angle(h0.position) * TMath::RadToDeg();
                        }
                        else if (h1.index < 4 && h2.index < 4) {
                            hitAngle = h2.position.Angle(h1.position) * TMath::RadToDeg();
                        }
                        else hitAngle = NAN;

                        v_hitAng->add(hitAngle);

                        mul++;

                    }
                }
            }

        }

        /*
        for (size_t i = 0; i < mult; i++) {
            auto h0 = hits[i];

            for (size_t j = i + 1; i < mult; i++) {
                auto h1 = hits[j];

                //
                //removing ability to look at 2 coincidences
                //
                if (mult < 3) {
                    v_pos0->add(h0.position);
                    v_dir0->add(h0.direction);
                    v_theta0->add(h0.theta * TMath::RadToDeg());
                    v_pos1->add(h1.position);
                    v_dir1->add(h1.direction);
                    v_theta1->add(h1.theta * TMath::RadToDeg());
                    v_pos2->add(NAN_TVECTOR3);
                    v_dir2->add(NAN_TVECTOR3);
                    v_theta2->add(NAN);

                    v_Ea0->add(h0.Ea);
                    v_Et0->add(h0.Et);
                    v_Edep0->add(h0.Edep);
                    v_BE0->add(h0.BE);
                    v_FE0->add(h0.FE);
                    v_Ea1->add(h1.Ea);
                    v_Et1->add(h1.Et);
                    v_Edep1->add(h1.Edep);
                    v_BE1->add(h1.BE);
                    v_FE1->add(h1.FE);
                    v_Ea2->add(NAN);
                    v_Et2->add(NAN);
                    v_Edep2->add(NAN);
                    v_BE2->add(NAN);
                    v_FE2->add(NAN);

                    v_dE0->add(h0.dE);
                    v_dE1->add(h1.dE);
                    v_dE2->add(NAN);
                    v_ang0->add(h0.angle * TMath::RadToDeg());
                    v_ang1->add(h1.angle * TMath::RadToDeg());
                    v_ang2->add(NAN);

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
                    v_i2->add(static_cast<short>(-1));
                    v_F2->add(static_cast<short>(-1));
                    v_B2->add(static_cast<short>(-1));
                    v_FT2->add(NAN);
                    v_BT2->add(NAN);

                    double hitAngle;
                    if (h0.index < 4 && h1.index < 4) {
                        hitAngle = h1.position.Angle(h0.position) * TMath::RadToDeg();
                    }
                    else hitAngle = NAN;

                    v_hitAng->add(hitAngle);

                    mul++;
                }

                for (size_t k = j + 1; k < mult; k++) {
                    auto h2 = hits[k];

                    if(abs(h0.TF-h1.TF)>1500 || abs(h2.TF-h1.TF)>1500 || abs(h0.TF-h2.TF)>1500) continue;

                    v_pos0->add(h0.position);
                    v_dir0->add(h0.direction);
                    v_theta0->add(h0.theta * TMath::RadToDeg());
                    v_pos1->add(h1.position);
                    v_dir1->add(h1.direction);
                    v_theta1->add(h1.theta * TMath::RadToDeg());
                    v_pos2->add(h2.position);
                    v_dir2->add(h2.direction);
                    v_theta2->add(h2.theta * TMath::RadToDeg());


                    v_Ea0->add(h0.Ea);
                    v_Et0->add(h0.Et);
                    v_Edep0->add(h0.Edep);
                    v_BE0->add(h0.BE);
                    v_FE0->add(h0.FE);
                    v_Ea1->add(h1.Ea);
                    v_Et1->add(h1.Et);
                    v_Edep1->add(h1.Edep);
                    v_BE1->add(h1.BE);
                    v_FE1->add(h1.FE);
                    v_Ea2->add(h2.Ea);
                    v_Et2->add(h2.Et);
                    v_Edep2->add(h2.Edep);
                    v_BE2->add(h2.BE);
                    v_FE2->add(h2.FE);

                    v_dE0->add(h0.dE);
                    v_dE1->add(h1.dE);
                    v_dE2->add(h2.dE);
                    v_ang0->add(h0.angle * TMath::RadToDeg());
                    v_ang1->add(h1.angle * TMath::RadToDeg());
                    v_ang2->add(h2.angle * TMath::RadToDeg());

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
                    v_i2->add(static_cast<short>(h2.index));
                    v_F2->add(h2.fseg);
                    v_B2->add(h2.bseg);
                    v_FT2->add(h2.TF);
                    v_BT2->add(h2.TB);

                    double hitAngle;
                    if (h0.index < 4 && h1.index < 4) {
                        hitAngle = h1.position.Angle(h0.position) * TMath::RadToDeg();
                    }
                    else if (h0.index < 4 && h2.index < 4) {
                        hitAngle = h2.position.Angle(h0.position) * TMath::RadToDeg();
                    }
                    else if (h1.index < 4 && h2.index < 4) {
                        hitAngle = h2.position.Angle(h1.position) * TMath::RadToDeg();
                    }
                    else hitAngle = NAN;

                    v_hitAng->add(hitAngle);

                    mul++;
                }
            }
        }
        */
    }

    void SortByType(Hit h0, Hit h1, Hit h2) {
        auto id0 = h0.index; auto id1 = h1.index; auto id2 = h2.index;
        if((id0 > 4 && id1 > 4) || (id0 > 4 && id2 > 4)) {}
    }

    static TLorentzVector constructBeamVector(const Ion &beam,
                                              const Ion &targetIon,
                                              double beamEnergy) {
        TLorentzVector plbeam(TVector3(0, 0, sqrt(2 * beamEnergy * beam.getMass())), beamEnergy + beam.getMass());
        TLorentzVector pltarget(TVector3(0, 0, 0), targetIon.getMass());
        return plbeam + pltarget;
    }

    void terminate() override {
        AbstractSortedAnalyzer::terminate();
        gDirectory->WriteTObject(t);
    }

    void clear() {
        mul = 0;
        AUSA::clear(
                *v_Et0, *v_Ea0, *v_theta0, *v_Edep0, *v_Et1, *v_Ea1, *v_theta1, *v_Edep1,
                *v_i0, *v_FE0, *v_BE0, *v_i1, *v_FE1, *v_BE1,
                *v_F0, *v_B0, *v_Ecm0, *v_F1, *v_B1, *v_Ecm1,
                *v_ang0, *v_pos0, *v_dir0, *v_ang1, *v_pos1, *v_dir1,
                *v_dE0, *v_FT0, *v_BT0, *v_dE1, *v_FT1, *v_BT1, *v_hitAng,
                *v_dir2, *v_pos2, *v_Edep2, *v_Ea2, *v_Et2, *v_BE2, *v_FE2, *v_theta2, *v_dE2,
                *v_Ecm2, *v_i2, *v_F2, *v_B2, *v_ang2, *v_FT2, *v_BT2
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
    unique_ptr<DynamicBranchVector<TVector3>> v_dir0, v_pos0, v_dir1, v_pos1, v_dir2, v_pos2;
    unique_ptr<DynamicBranchVector<double>> v_Edep0, v_Edep1, v_Edep2;
    unique_ptr<DynamicBranchVector<double>> v_Ea0, v_Et0, v_BE0, v_FE0, v_theta0, v_dE0, v_Ecm0;
    unique_ptr<DynamicBranchVector<double>> v_Ea1, v_Et1, v_BE1, v_FE1, v_theta1, v_dE1, v_Ecm1;
    unique_ptr<DynamicBranchVector<double>> v_Ea2, v_Et2, v_BE2, v_FE2, v_theta2, v_dE2, v_Ecm2;
    unique_ptr<DynamicBranchVector<short>> v_i0, v_i1, v_i2;
    unique_ptr<DynamicBranchVector<short>> v_F0, v_B0, v_F1, v_B1, v_F2, v_B2;
    unique_ptr<DynamicBranchVector<double>> v_ang0, v_ang1, v_ang2, v_hitAng;
    unique_ptr<DynamicBranchVector<double>> v_FT0, v_BT0, v_FT1, v_BT1, v_FT2, v_BT2;


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
        output_dir = EUtil::getProjectRoot() + "/" + cfg.lookup("data_output_dir_3coincidence").c_str();
    } else {
        setup_path = cfg.lookup("setup_file").c_str();
        target_path = cfg.lookup("target_file").c_str();
        input_path = cfg.lookup("data_input").c_str();
        output_dir = cfg.lookup("data_output_dir_3coincidence ").c_str();
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
