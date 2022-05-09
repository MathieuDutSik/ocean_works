// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Model_grids.h"
int main(int argc, char *argv[]) {
  srand_random_set();
  SingletonTime time1;
  try {
    if (argc != 3) {
      std::cerr << "GRID_ComputeWindRot2grid [fileIN] [fileOUT]\n";
      std::cerr << "with fileIN  the input  grid\n";
      std::cerr << " and fileOUT the output grid\n";
      return -1;
    }
    std::string GridFileIN = argv[1];
    std::string GridFileOUT = argv[2];
    std::cerr << " GridFileIN = " << GridFileIN << "\n";
    std::cerr << "GridFileOUT = " << GridFileOUT << "\n";
    std::string BoundFile = "unset";
    GridArray GrdArr = ReadUnstructuredGrid(GridFileIN, BoundFile);
    if (!GrdArr.GrdArrRho.DEP) {
      std::cerr << "We should have DEP initialized\n";
      throw TerminalException{1};
    }
    MyMatrix<double> &DEP = *GrdArr.GrdArrRho.DEP;
    int nbVert = DEP.size();
    for (int iVert = 0; iVert < nbVert; iVert++)
      DEP(iVert, 0) = 0;
    WriteUnstructuredGrid(GridFileOUT, GrdArr);
    std::cerr << "Normal termination of GRID_ComputeWindRot2grid\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRID_ComputeWindRot2grid\n";
    exit(e.eVal);
  }
  runtime(time1);
}
