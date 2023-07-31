// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Model_grids.h"
int main(int argc, char *argv[]) {
  srand_random_set();
  HumanTime time1;
  try {
    if (argc != 3) {
      std::cerr << "GRID_CreateTrivialBoundary [BoundFileIN] [BoundOUT]\n";
      std::cerr << "with GridFileIN    the input grid\n";
      std::cerr << "with BoundOUT   the input boundary with a gr3 extension\n";
      return -1;
    }
    std::string GridFileIN = argv[1];
    std::string BoundFileIN = "unset";
    std::string BoundFileOUT = argv[2];
    std::cerr << "  GridFileIN = " << GridFileIN << "\n";
    std::cerr << "BoundFileOUT = " << BoundFileOUT << "\n";
    GridArray GrdArr = ReadUnstructuredGrid(GridFileIN, BoundFileIN);
    GrdArr.IOBP = GetTrivialIOBP(GrdArr);
    WriteWWMboundaryGR3(BoundFileOUT, GrdArr);
    std::cerr << "Normal termination of GRID_CreateTrivialBoundary\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRID_CreateTrivialBoundary\n";
    exit(e.eVal);
  }
  runtime(time1);
}
