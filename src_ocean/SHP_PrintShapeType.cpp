// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Basic_shapelib.h"

int main(int argc, char *argv[]) {
  srand_random_set();
  SingletonTime time1;
  try {
    if (argc != 2) {
      std::cerr << "SHP_PrintShapeType [fileIn.shp]\n";
      return -1;
    }
    std::string eFileName = argv[1];
    ShapefileData eShape = SHP_ReadShapefile(eFileName);
    std::cerr << "After SHP_ReadShapeFile\n";
    PrintShapefileData(std::cout, eShape);
    std::cerr << "Normal termination of SHP_PrintShapeType\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in SHP_PrintShapeType\n";
    exit(e.eVal);
  }
  runtime(time1);
}
