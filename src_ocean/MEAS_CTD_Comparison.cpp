// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Satellite.h"

int main(int argc, char *argv[]) {
  srand_random_set();
  HumanTime time1;
  try {
    FullNamelist eFull = NAMELIST_GetStandardCTD_COMPARISON();
    if (argc != 2) {
      std::cerr << "MEAS_CTD_Comparison [file.nml]\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    Process_ctd_Comparison_Request(eFull);
    std::cerr << "Normal termination of MEAS_CTD_Comparison\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in MEAS_CTD_Comparison\n";
    exit(e.eVal);
  }
  runtime(time1);
}
