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
    int eShapeType=SHP_GetShapeType(eFileName);
    std::cerr << "eFileName=" << eFileName << " eShapeType=" << eShapeType << "\n";
    std::cerr << "Normal Termination of SHP_GetShapeType\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Error in SHP_GetShapeType\n";
    exit(e.eVal);
  }
}
