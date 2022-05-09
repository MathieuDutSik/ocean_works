// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Satellite.h"

int main(int argc, char *argv[]) {
  srand_random_set();
  SingletonTime time1;
  try {
    FullNamelist eFull = NAMELIST_GetStandardSST_COMPARISON();
    if (argc != 2) {
      std::cerr << "MEAS_SST_Comparison [file.nml]\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    Process_sst_Comparison_Request(eFull);
    std::cerr << "Normal termination of SAT_AltimeterComparison\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in SAT_AltimeterComparison\n";
    exit(e.eVal);
  }
  runtime(time1);
}
