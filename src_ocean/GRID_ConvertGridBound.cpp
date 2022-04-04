#include "Model_grids.h"
int main(int argc, char *argv[]) {
  srand_random_set();
  SingletonTime time1;
  try {
    if (argc != 5) {
      std::cerr
          << "GRID_ConvertGridBound [fileIN] [boundIN] [fileOUT] [boundOUT]\n";
      std::cerr << " with fileIN  the input  grid\n";
      std::cerr << "     boundIN  the input  boundary\n";
      std::cerr << "     fileOUT  the output grid\n";
      std::cerr << "and boundOUT  the output boundary\n";
      return -1;
    }
    std::string GridFileIN = argv[1];
    std::string BoundFileIN = argv[2];
    std::string GridFileOUT = argv[3];
    std::string BoundFileOUT = argv[4];
    std::cerr << "  GridFileIN = " << GridFileIN << "\n";
    std::cerr << " BoundFileIN = " << BoundFileIN << "\n";
    std::cerr << " GridFileOUT = " << GridFileOUT << "\n";
    std::cerr << "BoundFileOUT = " << BoundFileOUT << "\n";
    GridArray GrdArr = ReadUnstructuredGrid(GridFileIN, BoundFileIN);
    WriteUnstructuredGridTot(GridFileOUT, BoundFileOUT, GrdArr);
    std::cerr << "Normal termination of GRID_ConvertGridBound\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRID_ConvertGridBound\n";
    exit(e.eVal);
  }
  runtime(time1);
}
