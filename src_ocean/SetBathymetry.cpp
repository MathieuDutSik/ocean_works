#include "Model_grids.h"
int main(int argc, char *argv[])
{
  try {
    if (argc != 3) {
      std::cerr << "SetBathymetry is used as\n";
      std::cerr << "SetBathymetry [fileGrid] [TheDep]\n";
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
    std::cerr << "Normal termination of SetBathymetry\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Error in SetBathymetry\n";
    exit(e.eVal);
  }
}
