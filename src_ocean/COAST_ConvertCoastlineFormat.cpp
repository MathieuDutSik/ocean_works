#include "Basic_shapelib.h"

int main(int argc, char *argv[]) {
  srand_random_set();
  std::chrono::time_point<std::chrono::system_clock> time1 =
      std::chrono::system_clock::now();
  try {
    if (argc != 2) {
      std::cerr << "COAST_ConvertCoastlineFormat [fileIn] [fileOut]\n";
      return -1;
    }
    std::string eFileI = argv[1];
    std::string eFileO = argv[2];
    ShapefileData eShape = ReadCoastlineInformation(eFileI);
    WriteCoastlineInformation(eFileO, eShape);
    std::cerr << "Normal termination of COAST_ConvertCoastlineFormat\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in COAST_ConvertCoastlineFormat\n";
    exit(e.eVal);
  }
  std::chrono::time_point<std::chrono::system_clock> time2 =
      std::chrono::system_clock::now();
  std::cerr
      << "runtime = "
      << std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count()
      << "\n";
}
