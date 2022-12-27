// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "River.h"
int main(int argc, char *argv[]) {
  srand_random_set();
  SingletonTime time1;
  try {
    FullNamelist eFull = NAMELIST_GetStandard_ComputeRiverForcing_ROMS();
    if (argc != 2) {
      std::cerr << "RIVER_MakeRiverFile [file.nml]\n";
      std::cerr << "with file.nml the file describing the chosen options\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    CreateRiverFile(eFull);
    std::cerr << "Normal termination of RIVER_MakeRiverFile\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in RIVER_MakeRiverFile\n";
    exit(e.eVal);
  }
  runtime(time1);
}
