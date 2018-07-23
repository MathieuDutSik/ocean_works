#include "Model_grids.h"
int main(int argc, char *argv[])
{
  try {
    if (argc != 3) {
      std::cerr << "ComputeWindRot2grid is used as\n";
      std::cerr << "ComputeWindRot2grid [fileIN] [fileOUT]\n";
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
    std::cerr << "Normal termination of ComputeWindRot2grid\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
