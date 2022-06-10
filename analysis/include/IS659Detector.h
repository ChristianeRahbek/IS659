//
// Created by erik on 11/05/2021.
//

#ifndef IS659_DETECTOR_H
#define IS659_DETECTOR_H

#include <utility>

namespace IS659Detector_enums {
  enum IS659Type {
    dsd,
    plast,
    clov,
    indie
  };
}

using namespace std;
using namespace IS659Detector_enums;

class IS659Detector {

public:
  IS659Detector(unsigned short id, string name, IS659Type type, string partner = "") :
  id(id), name(move(name)), type(type), partner(move(partner)) {}

  unsigned short getId() const { return id; }
  string getName() const { return name; }
  IS659Type getType() const { return type; }
  string getPartner() const { return partner; }
  bool hasPartner() const { return !partner.empty(); }

private:
  unsigned short id;
  string name;
  IS659Type type;
  string partner;
};

#endif //IS659_DETECTOR_H