// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Model_grids.h"

int main(int argc, char *argv[]) {
  srand_random_set();
  SingletonTime time1;
  try {
    FullNamelist eFull = NAMELIST_NodeElimination();
    if (argc != 2) {
      std::cerr << "GRID_EliminateNodes File.nml\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);

    SingleBlock const& BlockPROC = eFull.ListBlock.at("PROC");
    std::string GridFileIn = BlockPROC.ListStringValues.at("GridFileIn");
    std::string GridFileOut = BlockPROC.ListStringValues.at("GridFileOut");
    std::vector<double> ListLon = BlockPROC.ListListDoubleValues.at("lon");
    std::vector<double> ListLat = BlockPROC.ListListDoubleValues.at("lat");
    bool keep_biggest = BlockPROC.ListBoolValues.at("KeepBiggestConnected");
    GridArray GrdArr = ReadUnstructuredGrid(GridFileIn, "unset");
    double lon0 = ListLon[0];
    double lon1 = ListLon[1];
    double lat0 = ListLat[0];
    double lat1 = ListLat[1];
    double d_lon0 = lon1 - lon0;
    double d_lat0 = lat1 - lat0;
    int mnp = GrdArr.GrdArrRho.LON.rows();
    std::vector<int> ListStatus(mnp);
    int sumPlus = 0;
    for (int i=0; i<mnp; i++) {
      double lon = GrdArr.GrdArrRho.LON(i,0);
      double lat = GrdArr.GrdArrRho.LAT(i,0);
      double d_lon1 = lon - lon0;
      double d_lat1 = lat - lat0;
      double det = d_lon1 * d_lat0 - d_lon0 * d_lat1;
      int val;
      if (det > 0)
        val = 1;
      else
        val = 0;
      sumPlus += val;
      ListStatus[i] = val;
    }
    std::cerr << "Conclusion of refining mnp=" << mnp << " sumPlus=" << sumPlus << "\n";
    GridArray GrdArrRed = SelectSubsetVertices(GrdArr, ListStatus);
    if (keep_biggest) {
      GrdArrRed = KeepLargestConnectedComponent(GrdArrRed);
    }
    WriteUnstructuredGrid(GridFileOut, GrdArrRed);
    std::cerr << "Normal termination of GRID_EliminateNodes\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRID_EliminateNodes\n";
    exit(e.eVal);
  }
  runtime(time1);
}
