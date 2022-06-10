cmake_minimum_required(VERSION 3.15)
project(calibration)

find_program(Calibrator_FOUND Calibrator)
if(NOT Calibrator_FOUND)
  message(WARNING "System-wide installation of Calibrator from the Subatomic Group at Aarhus University is required "
                  " for the 'calibration' subproject to work. Get it here: https://gitlab.au.dk/ausa/calibrator. "
                  "Skipping this sub project for now...")
  return()
endif()
link_libraries(AUSA::AUSA config++)

# TODO: script calibration of files
# Calibrator -s setup/setup.json -k setup/calibration/4a_simple.json -r 2 -b 13 -d U1 --gui data/unpacked/Run032.root
# see also top of *.cal files in setup dir