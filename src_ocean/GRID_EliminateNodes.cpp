// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Model_grids.h"

int main(int argc, char *argv[]) {
  srand_random_set();
  HumanTime time1;
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
    std::string SegmentFile = BlockPROC.ListStringValues.at("SegmentFile");
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
    std::cerr << "GridFileOut=" << GridFileOut << "\n";
    if (GridFileOut != "unset") {
      std::cerr << "Writing to file\n";
      WriteUnstructuredGrid(GridFileOut, GrdArrRed);
    }
    //
    if (SegmentFile != "unset") {
      int mnp_red = GrdArr.GrdArrRho.LON.rows();
      auto get_nearest_node=[&](double const& e_lon, double const& e_lat) -> int {
        double min_dist = std::numeric_limits<double>::max();
        int i_near = -1;
        for (int i=0; i<mnp_red; i++) {
          double delta_lon = std::abs(GrdArr.GrdArrRho.LON(i,0) - e_lon);
          double delta_lat = std::abs(GrdArr.GrdArrRho.LAT(i,0) - e_lat);
          double dist = delta_lon + delta_lat;
          if (dist < min_dist) {
            min_dist = dist;
            i_near = i;
          }
        }
        return i_near;
      };
      std::vector<std::vector<int>> ListCyc = GetListBoundaryCycles(GrdArrRed.INE, mnp_red);
      std::cerr << "|ListCyc|=" << ListCyc.size() << "\n";
      std::vector<int> eCyc = GetLongestCycle(ListCyc);
      std::cerr << "|eCyc|=" << eCyc.size() << "\n";
      int idx0 = get_nearest_node(lon0, lat0);
      int idx1 = get_nearest_node(lon1, lat1);
      std::vector<int> eSegment = GetShortestSegment(eCyc, idx0, idx1);
      //
      std::ofstream os(SegmentFile);
      size_t len = eSegment.size();
      os << len << "\n";
      for (size_t u=0; u<len; u++) {
        int pos = eSegment[u];
        double eLon = GrdArr.GrdArrRho.LON(pos,0);
        double eLat = GrdArr.GrdArrRho.LAT(pos,0);
        os << " " << eLon << " " << eLat << "\n";
      }
    }
    std::cerr << "Normal termination of GRID_EliminateNodes\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRID_EliminateNodes\n";
    exit(e.eVal);
  }
  runtime(time1);
}
