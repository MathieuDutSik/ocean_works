#include "Model_grids.h"
int main(int argc, char *argv[]) {
  srand_random_set();
  std::chrono::time_point<std::chrono::system_clock> time1 =
      std::chrono::system_clock::now();
  try {
    if (argc != 5) {
      std::cerr << "GRID_ScaleXYgrid [GridFileIN] [Spherical] [GridFileOUT] "
                   "[FactMult]\n";
      std::cerr << "with GridFileIN    the input grid\n";
      std::cerr << "with Spherical being 0 or 1\n";
      std::cerr << " and GridFileOUT   the output grid\n";
      std::cerr << " BathyChange values are:\n";
      std::cerr << "  0: no operation\n";
      std::cerr << "  1: changing the sign of the bathymetry\n";
      std::cerr << "  2: setting up the bathymetry to constant equal to zero\n";
      return -1;
    }
    std::string GridFileIN = argv[1];
    int Spherical;
    sscanf(argv[2], "%d", &Spherical);
    std::string GridFileOUT = argv[3];
    double FactMult;
    sscanf(argv[4], "%lf", &FactMult);
    std::cerr << " GridFileIN = " << GridFileIN << "\n";
    std::cerr << "GridFileOUT = " << GridFileOUT << "\n";
    std::cerr << "   FactMult = " << FactMult << "\n";
    std::string BoundFileIN = "unset";
    GridArray GrdArr = ReadUnstructuredGrid(GridFileIN, BoundFileIN);
    double minLon = GrdArr.GrdArrRho.LON.minCoeff();
    double maxLon = GrdArr.GrdArrRho.LON.maxCoeff();
    double minLat = GrdArr.GrdArrRho.LAT.minCoeff();
    double maxLat = GrdArr.GrdArrRho.LAT.maxCoeff();
    double distLon1, distLon2, distLat1, distLat2;
    if (Spherical == 1) {
      distLon1 = GeodesicDistanceKM(minLon, minLat, maxLon, minLat);
      distLon2 = GeodesicDistanceKM(minLon, maxLat, maxLon, maxLat);
      distLat1 = GeodesicDistanceKM(minLon, minLat, minLon, maxLat);
      distLat2 = GeodesicDistanceKM(maxLon, minLat, maxLon, maxLat);
    } else {
      distLon1 = (maxLon - minLon) / double(1000);
      distLon2 = distLon1;
      distLat1 = (maxLat - minLat) / double(1000);
      distLat2 = distLat1;
    }
    std::cerr << "Distances from lowest longitude to highest 1: " << distLon1
              << " 2: " << distLon2 << "\n";
    std::cerr << "Distances from lowest  latitude to highest 1: " << distLat1
              << " 2: " << distLat2 << "\n";
    //
    int nbNode = GrdArr.GrdArrRho.LON.size();
    for (int iNode = 0; iNode < nbNode; iNode++) {
      GrdArr.GrdArrRho.LON(iNode) =
          FactMult * (GrdArr.GrdArrRho.LON(iNode) - minLon);
      GrdArr.GrdArrRho.LAT(iNode) =
          FactMult * (GrdArr.GrdArrRho.LAT(iNode) - minLat);
    }
    //
    WriteUnstructuredGrid(GridFileOUT, GrdArr);
    std::cerr << "Normal termination of GRID_ScaleXYgrid\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRID_ScaleXYgrid\n";
    exit(e.eVal);
  }
  std::chrono::time_point<std::chrono::system_clock> time2 =
      std::chrono::system_clock::now();
  std::cerr
      << "runtime = "
      << std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count()
      << "\n";
}
