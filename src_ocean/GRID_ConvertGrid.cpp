// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Model_grids.h"
int main(int argc, char *argv[]) {
  srand_random_set();
  HumanTime time1;
  try {
    if (argc != 5) {
      std::cerr << "GRID_ConvertGrid [GridFileIN] [BoundFileIN] [GridFileOUT] "
                   "[BathyChange]\n";
      std::cerr << "with GridFileIN    the input grid\n";
      std::cerr << "with BoundFileIN   the input boundary (put unset if not "
                   "available)\n";
      std::cerr << " and GridFileOUT   the output grid\n";
      std::cerr << " Operation being done:\n";
      std::cerr << "  0    : no operation (just file format conversion)\n";
      std::cerr << "  1    : changing the sign of the bathymetry\n";
      std::cerr << "  2    : setting up the bathymetry to constant equal to zero\n";
      std::cerr << "  3:KM : Reduction by minimum distance. Distance is in km\n";
      return -1;
    }
    std::string GridFileIN = argv[1];
    std::string BoundFileIN = argv[2];
    std::string GridFileOUT = argv[3];
    std::string operation = argv[4];
    std::cerr << " GridFileIN = " << GridFileIN << "\n";
    std::cerr << "GridFileOUT = " << GridFileOUT << "\n";
    std::cerr << "operation = " << operation << "\n";
    GridArray GrdArr = ReadUnstructuredGrid(GridFileIN, BoundFileIN);
    double minLon = GrdArr.GrdArrRho.LON.minCoeff();
    double maxLon = GrdArr.GrdArrRho.LON.maxCoeff();
    double minLat = GrdArr.GrdArrRho.LAT.minCoeff();
    double maxLat = GrdArr.GrdArrRho.LAT.maxCoeff();
    std::cerr << "LON(min/max)=" << minLon << " / " << maxLon << "\n";
    std::cerr << "LAT(min/max)=" << minLat << " / " << maxLat << "\n";
    if (!GrdArr.GrdArrRho.DEP) {
      std::cerr << "DEP is not assigned\n";
      throw TerminalException{1};
    }
    auto writing=[&]() -> void {
      if (operation == "0") {
        return WriteUnstructuredGrid(GridFileOUT, GrdArr);
      }
      if (operation == "1") {
        MyMatrix<double> &DEP = *GrdArr.GrdArrRho.DEP;
        int nbNode = DEP.size();
        for (int iNode = 0; iNode < nbNode; iNode++)
          DEP(iNode) = -DEP(iNode);
        return WriteUnstructuredGrid(GridFileOUT, GrdArr);
      }
      if (operation == "2") {
        MyMatrix<double> &DEP = *GrdArr.GrdArrRho.DEP;
        int nbNode = DEP.size();
        for (int iNode = 0; iNode < nbNode; iNode++)
          DEP(iNode) = 0;
        return WriteUnstructuredGrid(GridFileOUT, GrdArr);
      }
      std::optional<std::string> opt_dist_reduction =
        get_postfix(operation, "3:");
      if (opt_dist_reduction) {
        std::string str_distkm = *opt_dist_reduction;
        std::cerr << "str_distkm=" << str_distkm << "\n";
        double CritDistKM = ParseScalar<double>(str_distkm);
        std::cerr << "CritDistKM=" << CritDistKM << "\n";
        GridArray GrdArrRet = MergeNeighboringVertices(GrdArr, CritDistKM);
        return WriteUnstructuredGrid(GridFileOUT, GrdArrRet);
      }
      std::cerr << "Failed to find a matching entry\n";
      throw TerminalException{1};
    };
    writing();
    std::cerr << "Normal termination of GRID_ConvertGrid\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRID_ConvertGrid\n";
    exit(e.eVal);
  }
  runtime(time1);
}
