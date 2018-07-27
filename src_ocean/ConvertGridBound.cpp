#include "Model_grids.h"
int main(int argc, char *argv[])
{
  try {
    if (argc != 5) {
      std::cerr << "ConvertGridBound is used as\n";
      std::cerr << "ConvertGridBound [fileIN] [boundIN] [fileOUT] [boundOUT]\n";
      std::cerr << " with fileIN  the input  grid\n";
      std::cerr << "     boundIN  the input  boundary\n";
      std::cerr << "     fileOUT  the output grid\n";
      std::cerr << "and boundOUT  the output boundary\n";
      return -1;
    }
    std::string GridFileIN  = argv[1];
    std::string BoundFileIN = argv[2];
    std::string GridFileOUT = argv[3];
    std::string BoundFileOUT= argv[4];
    std::cerr << "  GridFileIN = " << GridFileIN << "\n";
    std::cerr << " BoundFileIN = " << BoundFileIN << "\n";
    std::cerr << " GridFileOUT = " << GridFileOUT << "\n";
    std::cerr << "BoundFileOUT = " << BoundFileOUT << "\n";
    GridArray GrdArr=ReadUnstructuredGrid(GridFileIN, BoundFileIN);
    WriteUnstructuredGridTot(GridFileOUT, BoundFileOUT, GrdArr);
    std::cerr << "Normal termination of ConvertGridBound\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Error in ConvertGridBound\n";
    exit(e.eVal);
  }
}
