// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Floats.h"
int main(int argc, char *argv[]) {
  srand_random_set();
  HumanTime time1;
  try {
    FullNamelist eFull = NAMELIST_GetStandardPLOT_DRIFTER_TRACK();
    if (argc != 2) {
      std::cerr << "PLOT_drifters [file.nml]\n";
      std::cerr << "with file.nml the file describing the plotting routines\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    ICHTHYOP_PlotTrajectories(eFull);
    std::cerr << "Normal termination of PLOT_drifters\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in PLOT_drifters\n";
    exit(e.eVal);
  }
  runtime(time1);
}
