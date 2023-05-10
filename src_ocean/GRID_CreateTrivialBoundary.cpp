// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Model_grids.h"
int main(int argc, char *argv[]) {
  srand_random_set();
  HumanTime time1;
  try {
    if (argc != 3) {
      std::cerr << "GRID_ConvertGrid [BoundFileIN] [BoundOUT]\n"
      std::cerr << "with GridFileIN    the input grid\n";
      std::cerr << "with BoundFileIN   the input boundary (put unset if not "
                   "available)\n";
      std::cerr << " and GridFileOUT   the output grid\n";
      std::cerr << " BathyChange values are:\n";
      std::cerr << "  0: no operation\n";
      std::cerr << "  1: changing the sign of the bathymetry\n";
      std::cerr << "  2: setting up the bathymetry to constant equal to zero\n";
      return -1;
    }
    std::string GridFileIN = argv[1];
    std::string BoundFileIN = "unset";
    std::string GridFileOUT = argv[2];
    std::cerr << " GridFileIN = " << GridFileIN << "\n";
    std::cerr << "GridFileOUT = " << GridFileOUT << "\n";
    std::cerr << "BathyChange = " << BathyChange << "\n";
    if (BathyChange != 0 && BathyChange != 1 && BathyChange != 2) {
      std::cerr << "We should have SignChange = 0 or 1 or 2\n";
      throw TerminalException{1};
    }
    GridArray GrdArr = ReadUnstructuredGrid(GridFileIN, BoundFileIN);
    WriteUnstructuredGrid(GridFileOUT, GrdArr);
    std::cerr << "Normal termination of GRID_ConvertGrid\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRID_ConvertGrid\n";
    exit(e.eVal);
  }
  runtime(time1);
}
