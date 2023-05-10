// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Plotting_fct.h"
int main(int argc, char *argv[]) {
  srand_random_set();
  HumanTime time1;
  try {
    FullNamelist eFull = NAMELIST_CREATE_DEFAULT_SETTING();
    if (argc != 2) {
      std::cerr << "ROMS_CreateDefaultSetting [file.nml]\n";
      std::cerr
          << "with file.nml the file describing the default variable chosen\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    CreateDefaultInputFiles(eFull);
    std::cerr << "Normal termination of ROMS_CreateDefaultSetting\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in ROMS_CreateDefaultSetting\n";
    exit(e.eVal);
  }
  runtime(time1);
}
