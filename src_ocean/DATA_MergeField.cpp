// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Model_interpolation.h"
int main(int argc, char *argv[]) {
  srand_random_set();
  HumanTime time1;
  try {
    FullNamelist eFull = NAMELIST_GetStandardMODEL_MERGING();
    if (argc != 2) {
      std::cerr << "DATA_MergeField [file.nml]\n";
      std::cerr
          << "with file.nml the file describing the interpolation process\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    INTERPOL_field_Function(eFull);
    std::cerr << "Normal termination of DATA_MergeField\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in DATA_MergeField\n";
    exit(e.eVal);
  }
  runtime(time1);
}
