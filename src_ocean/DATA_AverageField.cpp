// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Model_interpolation.h"
int main(int argc, char *argv[]) {
  srand_random_set();
  SingletonTime time1;
  try {
    FullNamelist eFull = NAMELIST_InfileAveraging();
    if (argc != 2) {
      std::cerr << "DATA_AverageField [file.nml]\n";
      std::cerr
          << "with file.nml the file describing the interpolation process\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    Average_field_Function(eFull);
    std::cerr << "Normal termination of DATA_AverageField\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in DATA_AverageField\n";
    exit(e.eVal);
  }
  runtime(time1);
}
