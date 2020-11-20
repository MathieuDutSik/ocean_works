#include "Basic_shapelib.h"
#include "Model_grids.h"

int main(int argc, char *argv[])
{
  std::cerr << std::setprecision(12);
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
  try {
    if (argc != 3) {
      std::cerr << "COAST_ConvertMeshToCoast is used as\n";
      std::cerr << "COAST_ConvertMeshToCoast [fileIn] [fileOut]\n";
      return -1;
    }
    std::string eFileI=argv[1];
    std::string eFileO=argv[2];
    std::string BoundFile="unset";
    GridArray GrdArr=ReadUnstructuredGrid(eFileI, BoundFile);
    int nbVert=GrdArr.GrdArrRho.LON.rows();
    std::vector<std::vector<int>> eListList=GetListBoundaryCycles(GrdArr.INE, nbVert);
    ShapefileData eSHP=ExtractShapefileFromUnstructuredMesh(GrdArr, eListList);
    WriteCoastlineInformation(eFileO, eSHP);
    std::cerr << "Normal termination of COAST_ConvertMeshToCoast\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Error in COAST_ConvertMeshToCoast\n";
    exit(e.eVal);
  }
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "runtime = " << std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count() << "\n";
}
