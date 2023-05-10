// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "AquaSatellite.h"
int main(int argc, char *argv[]) {
  srand_random_set();
  std::cerr << std::fixed;
  std::cerr << std::setprecision(9);
  HumanTime time1;
  try {
    FullNamelist eFull = NAMELIST_GetStandardAQUA();
    if (argc != 2) {
      std::cerr << "SAT_AquaOperation [file.nml]\n";
      std::cerr
          << "with file.nml the file describing the input and plotting used\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    AquaDownloading(eFull);
    std::cerr << "Normal termination of SAT_AquaOperation\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in SAT_AquaOperation\n";
    exit(e.eVal);
  }
  runtime(time1);
}
