#include "Basic_shapelib.h"

int main(int argc, char *argv[])
{
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
}
