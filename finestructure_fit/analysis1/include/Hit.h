#include <TVector3.h>
#include <TLorentzVector.h>
#include "ParticleType.h"


struct Hit {
    double deposited;
    double paddeposited;

    double Ea, Et, FE, BE, Edssd, dE,
           EBeta, fbdiff, Ecm, Edep_alphas, EPlastic;


    double TF, TB, TPlastic, T;

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

    TLorentzVector lVector_alpha, lVector_triton, lVectorBeta, lVectorAlpha;
};
