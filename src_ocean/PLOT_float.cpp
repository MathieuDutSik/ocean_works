// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Floats.h"

int main(int argc, char *argv[]) {
  srand_random_set();
  HumanTime time1;
  try {
    FullNamelist eFull = NAMELIST_GetStandard_PlotRomsFloats();
    if (argc != 2) {
      std::cerr << "PLOT_float [file.nml]\n";
      std::cerr << "with file.nml the file describing the plotting routines\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    PLOT_ROMS_float(eFull);
    std::cerr << "Normal termination of PLOT_float\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in PLOT_float\n";
    exit(e.eVal);
  }
  runtime(time1);
}
