// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Satellite.h"
int main(int argc, char *argv[]) {
  srand_random_set();
  HumanTime time1;
  try {
    FullNamelist eFull = NAMELIST_Comparison_Altimetry_Source();
    if (argc != 2) {
      std::cerr << "SAT_ComparisonAltimetrySource [alti.nml]\n";
      std::cerr << "with alti.nml the file describing the chosen options\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    Process_Comparison_Altimetry_Sources(eFull);
    std::cerr << "Normal termination of SAT_ComparisonAltimetrySource\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in SAT_ComparisonAltimetrySource\n";
    exit(e.eVal);
  }
  runtime(time1);
}
