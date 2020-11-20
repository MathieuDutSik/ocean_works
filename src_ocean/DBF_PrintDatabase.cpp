#include "Basic_shapelib.h"

int main(int argc, char *argv[])
{
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
  try {
    if (argc != 2) {
      std::cerr << "SHP_GetShapeType is used as\n";
      std::cerr << "SHP_GetShapeType [fileIn.shp]\n";
      return -1;
    }
    std::string eFileName=argv[1];
    DBF_PrintDatabaseInfo(std::cout, eFileName);
    std::cerr << "Normal termination of DBF_PrintSatabase\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Error in DBF_PrintSatabase\n";
    exit(e.eVal);
  }
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "runtime = " << std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count() << "\n";
}
