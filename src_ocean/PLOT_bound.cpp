// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "ROMSboundary.h"
int main(int argc, char *argv[]) {
  srand_random_set();
  std::cerr << std::fixed;
  std::cerr << std::setprecision(9);
  SingletonTime time1;
  try {
    FullNamelist eFull = NAMELIST_GetStandardPLOT_BOUNDARY();
    if (argc != 2) {
      std::cerr << "PLOT_bound [file.nml]\n";
      std::cerr << "with file.nml the file describing the plotting routines\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    BOUND_Plotting_Function(eFull);
    std::cerr << "Normal termination of PLOT_bound\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in PLOT_bound\n";
    exit(e.eVal);
  }
  runtime(time1);
}
