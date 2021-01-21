#include "Model_grids.h"
int main(int argc, char *argv[])
{
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
  try {
    if (argc != 5) {
      std::cerr << "GRID_ConvertGrid [GridFileIN] [BoundFileIN] [GridFileOUT] [BathyChange]\n";
      std::cerr << "with GridFileIN    the input grid\n";
      std::cerr << "with BoundFileIN   the input boundary (put unset if not available)\n";
      std::cerr << " and GridFileOUT   the output grid\n";
      std::cerr << " BathyChange values are:\n";
      std::cerr << "  0: no operation\n";
      std::cerr << "  1: changing the sign of the bathymetry\n";
      std::cerr << "  2: setting up the bathymetry to constant equal to zero\n";
      return -1;
    }
    std::string GridFileIN  = argv[1];
    std::string BoundFileIN = argv[2];
    std::string GridFileOUT = argv[3];
    int BathyChange;
    sscanf(argv[4], "%d", &BathyChange);
    std::cerr << " GridFileIN = " << GridFileIN << "\n";
    std::cerr << "GridFileOUT = " << GridFileOUT << "\n";
    std::cerr << "BathyChange = " << BathyChange << "\n";
    if (BathyChange != 0 && BathyChange != 1 && BathyChange != 2) {
      std::cerr << "We should have SignChange = 0 or 1 or 2\n";
      throw TerminalException{1};
    }
    GridArray GrdArr=ReadUnstructuredGrid(GridFileIN, BoundFileIN);
    double minLon=GrdArr.GrdArrRho.LON.minCoeff();
    double maxLon=GrdArr.GrdArrRho.LON.maxCoeff();
    double minLat=GrdArr.GrdArrRho.LAT.minCoeff();
    double maxLat=GrdArr.GrdArrRho.LAT.maxCoeff();
    std::cerr << "LON(min/max)=" << minLon << " / " << maxLon << "\n";
    std::cerr << "LAT(min/max)=" << minLat << " / " << maxLat << "\n";
    if (BathyChange == 1) {
      int nbNode=GrdArr.GrdArrRho.DEP.size();
      for (int iNode=0; iNode<nbNode; iNode++)
	GrdArr.GrdArrRho.DEP(iNode) = - GrdArr.GrdArrRho.DEP(iNode);
    }
    if (BathyChange == 2) {
      int nbNode=GrdArr.GrdArrRho.DEP.size();
      for (int iNode=0; iNode<nbNode; iNode++)
	GrdArr.GrdArrRho.DEP(iNode) = 0;
    }
    WriteUnstructuredGrid(GridFileOUT, GrdArr);
    std::cerr << "Normal termination of GRID_ConvertGrid\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Error in GRID_ConvertGrid\n";
    exit(e.eVal);
  }
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "runtime = " << std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count() << "\n";
}
