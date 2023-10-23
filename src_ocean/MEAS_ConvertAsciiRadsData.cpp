// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Satellite.h"

int main(int argc, char *argv[]) {
  srand_random_set();
  HumanTime time1;
  try {
    if (argc != 3) {
      std::cerr << "MEAS_ConvertAsciiRadsData [PrefixIn] [PrefixOut]\n";
      return -1;
    }
    std::string PrefixI = argv[1];
    std::string PrefixO = argv[2];
    int method = 1;
    RadsAscToNetcdf(PrefixI, PrefixO, method);
    std::cerr << "Normal termination of MEAS_ConvertAsciiRadsData\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in MEAS_ConvertAsciiRadsData\n";
    exit(e.eVal);
  }
  runtime(time1);
}
