//
// Created by erik on 6/9/22.
//

#ifndef IS659_HIT_H
#define IS659_HIT_H

#include <TVector3.h>
#include <TLorentzVector.h>

struct Hit {
  unsigned short id;
  TVector3 dir, pos;
  double theta, phi, angle;
  double Edep;
  unsigned short FI, BI;
  double FE, BE;
  double FT, BT;
  double Ea, Et;
};

#endif //IS659_HIT_H