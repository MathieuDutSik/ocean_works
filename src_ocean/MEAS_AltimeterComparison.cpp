// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Satellite.h"

#ifdef DEBUG_TRAP_NAN
#include <fenv.h>
#endif

int main(int argc, char *argv[]) {
  srand_random_set();
#ifdef DEBUG_TRAP_NAN
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif
  HumanTime time1;
  try {
    FullNamelist eFull = NAMELIST_GetStandardALTIMETRY_COMPARISON();
    if (argc != 2) {
      std::cerr << "SAT_AltimeterComparison [alti.nml]\n";
      std::cerr << "with alti.nml the file describing the chosen options\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    Process_Altimetry_Comparison_Request(eFull);
    std::cerr << "Normal termination of SAT_AltimeterComparison\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in SAT_AltimeterComparison\n";
    exit(e.eVal);
  }
  runtime(time1);
}
