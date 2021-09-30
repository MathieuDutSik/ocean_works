#include "Model_grids.h"
int main(int argc, char *argv[])
{
  srand_random_set();
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
  try {
    if (argc != 3) {
      std::cerr << "GRID_SetBathymetry [fileGrid] [TheDep]\n";
      std::cerr << "with fileGrid  the grid file\n";
      std::cerr << " and TheDep the chosen bathymetry\n";
      return -1;
    }
    std::string GridFile = argv[1];
    std::string TheDepStr= argv[2];
    double TheDep;
    std::istringstream(TheDepStr) >> TheDep;
    std::cerr << " GridFile = " << GridFile << "\n";
    std::cerr << " TheDep=" << TheDep << "\n";
    std::string BoundFile="unset";
    GridArray GrdArr=ReadUnstructuredGrid(GridFile, BoundFile);
    int mnp=GrdArr.GrdArrRho.LON.size();
    for (int i=0; i<mnp; i++)
      GrdArr.GrdArrRho.DEP(i,0)=TheDep;
    WriteUnstructuredGrid(GridFile, GrdArr);
    std::cerr << "Normal termination of GRID_SetBathymetry\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Error in GRID_SetBathymetry\n";
    exit(e.eVal);
  }
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "runtime = " << std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count() << "\n";
}
