// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Basic_shapelib.h"
#include "Model_grids.h"
#include "Timings.h"

int main(int argc, char *argv[]) {
  srand_random_set();
  std::cerr << std::setprecision(12);
  SingletonTime time1;
  try {
    if (argc != 3) {
      std::cerr << "COAST_ConvertMeshToCoast [fileIn] [fileOut]\n";
      return -1;
    }
    std::string eFileI = argv[1];
    std::string eFileO = argv[2];
    std::string BoundFile = "unset";
    GridArray GrdArr = ReadUnstructuredGrid(eFileI, BoundFile);
    int nbVert = GrdArr.GrdArrRho.LON.rows();
    std::vector<std::vector<int>> eListList =
        GetListBoundaryCycles(GrdArr.INE, nbVert);
    ShapefileData eSHP =
        ExtractShapefileFromUnstructuredMesh(GrdArr, eListList);
    WriteCoastlineInformation(eFileO, eSHP);
    std::cerr << "Normal termination of COAST_ConvertMeshToCoast\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in COAST_ConvertMeshToCoast\n";
    exit(e.eVal);
  }
  runtime(time1);
}
