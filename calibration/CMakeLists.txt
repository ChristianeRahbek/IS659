cmake_minimum_required(VERSION 3.15)
project(calibration)

find_program(Calibrator_FOUND Calibrator)
if(NOT Calibrator_FOUND)
  message(WARNING "System-wide installation of Calibrator from the Subatomic Group at Aarhus University is required "
                  " for the 'calibration' sub project to work. Get it here: https://gitlab.au.dk/ausa/calibrator. "
                  "Skipping this sub project for now...")
endif()
link_libraries(AUSA::AUSA config++)