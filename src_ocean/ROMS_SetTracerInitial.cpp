// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Plotting_fct.h"
int main(int argc, char *argv[]) {
  srand_random_set();
  HumanTime time1;
  try {
    FullNamelist eFull = NAMELIST_SET_VARIABLE_INITIAL_ROMS();
    if (argc != 2) {
      std::cerr << "ROMS_SetTracerInitial [file.nml]\n";
      std::cerr << "with file.nml the file describing the plotting routines\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    SetNetcdfInitial(eFull);
    std::cerr << "Normal termination of ROMS_SetTracerInitial\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in ROMS_SetTracerInitial\n";
    exit(e.eVal);
  }
  runtime(time1);
}
