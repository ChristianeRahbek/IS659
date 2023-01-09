#include <TVector3.h>
#include <TLorentzVector.h>
#include "ParticleType.h"

struct Hit {
    double deposited;
    double paddeposited;

    double Ea, Et, FE, BE, Edssd, dE,
           EBeta, fbdiff, Ecm, Edep_alphas,
           Edep0, Edep1, Edep2, Edep3, /*Deposited energy in each detector*/
           Ea0, Ea1, Ea2, Ea3, Et0, Et1, Et2, Et3; //alpha and triton energies in each detector


    double TF, TB, TPad, T;

    double tarTrav, detTrav;
    double Ectarget;
    TVector3 direction, position, origin;
    double theta;
    double angle;
    double phi;
    int numInDet, detectorMul;
    short fseg, bseg;
    std::vector<ParticleType> type;
    bool canBeAlpha;
    bool canBeBeta;
    size_t index;

    TLorentzVector lVector, lVectorBeta, lVectorAlpha;
};
