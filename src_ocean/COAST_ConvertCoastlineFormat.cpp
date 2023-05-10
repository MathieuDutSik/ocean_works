// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Basic_shapelib.h"

int main(int argc, char *argv[]) {
  srand_random_set();
  HumanTime time1;
  try {
    if (argc != 2) {
      std::cerr << "COAST_ConvertCoastlineFormat [fileIn] [fileOut]\n";
      return -1;
    }
    std::string eFileI = argv[1];
    std::string eFileO = argv[2];
    ShapefileData eShape = ReadCoastlineInformation(eFileI);
    WriteCoastlineInformation(eFileO, eShape);
    std::cerr << "Normal termination of COAST_ConvertCoastlineFormat\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in COAST_ConvertCoastlineFormat\n";
    exit(e.eVal);
  }
  runtime(time1);
}
