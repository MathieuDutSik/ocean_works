// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Transect.h"
int main(int argc, char *argv[]) {
  srand_random_set();
  HumanTime time1;
  try {
    FullNamelist eFull = NAMELIST_GetStandard_PlotTransect();
    if (argc != 2) {
      std::cerr << "PLOT_transect [file.nml]\n";
      std::cerr << "with file.nml the file describing the plotting routines\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    TRANSECT_Plot(eFull);
    std::cerr << "Normal termination of PLOT_transect\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in PLOT_transect\n";
    exit(e.eVal);
  }
  runtime(time1);
}
