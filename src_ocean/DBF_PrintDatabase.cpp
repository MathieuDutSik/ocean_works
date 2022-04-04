#include "Basic_shapelib.h"

int main(int argc, char *argv[]) {
  srand_random_set();
  SingletonTime time1;
  try {
    if (argc != 2) {
      std::cerr << "DBF_PrintDatabase [fileIn.shp]\n";
      return -1;
    }
    std::string eFileName = argv[1];
    DBF_PrintDatabaseInfo(std::cout, eFileName);
    std::cerr << "Normal termination of DBF_PrintDatabase\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in DBF_PrintDatabase\n";
    exit(e.eVal);
  }
  runtime(time1);
}
