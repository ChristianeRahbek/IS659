cmake_minimum_required(VERSION 3.15)
project(sorting)

find_program(Sorter_FOUND Sorter)
if(NOT Calibrator_FOUND)
  message(WARNING "System-wide installation of Sorter from the Subatomic Group at Aarhus University is required "
                  " for the 'sorting' sub project to work. Get it here: https://gitlab.au.dk/ausa/sorter. "
                  "Skipping this subproject for now...")
  return()
endif()
link_libraries(AUSA::AUSA config++)

# TODO: script sorting of files
# Sorter -s setup/setup.json -m setup/matcher.json -o data/sorted/ data/unpacked/Run159.root