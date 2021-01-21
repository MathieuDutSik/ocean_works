#include "Model_grids.h"
int main(int argc, char *argv[])
{
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
  try {
    if (argc != 3) {
      std::cerr << "GRID_ComputeWindRot2grid [fileIN] [fileOUT]\n";
      std::cerr << "with fileIN  the input  grid\n";
      std::cerr << " and fileOUT the output grid\n";
      return -1;
    }
    std::string GridFileIN = argv[1];
    std::string GridFileOUT= argv[2];
    std::cerr << " GridFileIN = " << GridFileIN << "\n";
    std::cerr << "GridFileOUT = " << GridFileOUT << "\n";
    std::string BoundFile="unset";
    GridArray GrdArr=ReadUnstructuredGrid(GridFileIN, BoundFile);
    int nbVert=GrdArr.GrdArrRho.LON.size();
    for (int iVert=0; iVert<nbVert; iVert++)
      GrdArr.GrdArrRho.DEP(iVert,0)=0;
    WriteUnstructuredGrid(GridFileOUT, GrdArr);
    std::cerr << "Normal termination of GRID_ComputeWindRot2grid\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Error in GRID_ComputeWindRot2grid\n";
    exit(e.eVal);
  }
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "runtime = " << std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count() << "\n";
}
