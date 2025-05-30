// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Model_interpolation.h"
int main(int argc, char *argv[]) {
  srand_random_set();
  std::cerr << std::fixed;
  std::cerr << std::setprecision(9);
  HumanTime time1;
  try {
    FullNamelist eFull = NAMELIST_GetStandard_CREATE_sflux();
    if (argc != 2) {
      std::cerr << "DATA_CreateSflux [file.nml]\n";
      std::cerr << "with file.nml the file describing the choices made\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    CREATE_sflux_files(eFull);
    std::cerr << "Normal termination of DATA_CreateSflux\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in DATA_CreateSflux\n";
    exit(e.eVal);
  }
  runtime(time1);
}
