//
// Created by Christiane on 27-02-2023.
// Inspired by https://gitlab.au.dk/ausa/erik/is507_v3/-/blob/master/analysis/include/IS507Detector.h
//

#ifndef ANALYSIS1_IS659DETECTOR_H
#define ANALYSIS1_IS659DETECTOR_H

#include <utility>
#include <string>

namespace IS659Detector_enums {
    enum IS659Type{
        SquareDSSD,
        Plastic
    };
}

using namespace std;
using namespace IS659Detector_enums;

class IS659Detector {
public:
    IS659Detector(unsigned short id, string name, IS659Type detType):
    id(id), name(move(name)), detType(detType) {}

    unsigned short getId() const {return id;}
    string getName() const {return name;}
    IS659Type getDetType() const {return detType;}

private:
    unsigned short id;
    string name;
    IS659Type detType;
};

#endif //ANALYSIS1_IS659DETECTOR_H
