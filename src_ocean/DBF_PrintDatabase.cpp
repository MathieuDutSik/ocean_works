#include "Basic_shapelib.h"

int main(int argc, char *argv[])
{
  srand_random_set();
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
  try {
    if (argc != 2) {
      std::cerr << "DBF_PrintDatabase [fileIn.shp]\n";
      return -1;
    }
    std::string eFileName=argv[1];
    DBF_PrintDatabaseInfo(std::cout, eFileName);
    std::cerr << "Normal termination of DBF_PrintDatabase\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Error in DBF_PrintDatabase\n";
    exit(e.eVal);
  }
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "runtime = " << std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count() << "\n";
}
