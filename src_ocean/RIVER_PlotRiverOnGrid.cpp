// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "River.h"
int main(int argc, char *argv[]) {
  srand_random_set();
  std::cerr << std::fixed;
  std::cerr << std::setprecision(9);
  HumanTime time1;
  try {
    FullNamelist eFull = NAMELIST_PLOT_River();
    if (argc != 2) {
      std::cerr << "RIVER_PlotRiverOnGrid [file.nml]\n";
      std::cerr << "with file.nml the file describing the plotting routines\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    //
    PlotRiverInformation(eFull);
    std::cerr << "Normal termination of RIVER_PlotRiverOnGrid\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in RIVER_PlotRiverOnGrid\n";
    exit(e.eVal);
  }
  runtime(time1);
}
