#include "Basic_shapelib.h"

int main(int argc, char *argv[])
{
  try {
    if (argc != 2) {
      std::cerr << "SHP_PrintShapeType is used as\n";
      std::cerr << "SHP_PrintShapeType [fileIn.shp]\n";
      return -1;
    }
    std::string eFileName=argv[1];
    ShapefileData eShape=SHP_ReadShapefile(eFileName);
    std::cerr << "After SHP_ReadShapeFile\n";
    PrintShapefileData(std::cout, eShape);
    std::cerr << "Normal termination of SHP_PrintShapeType\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Error in SHP_PrintShapeType\n";
    exit(e.eVal);
  }
}
