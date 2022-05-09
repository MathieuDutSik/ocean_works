// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Basic_shapelib.h"

int main(int argc, char *argv[]) {
  srand_random_set();
  SingletonTime time1;
  try {
    if (argc != 2) {
      std::cerr << "SHP_GetShapeType [fileIn.shp]\n";
      return -1;
    }
    std::string eFileName = argv[1];
    int eShapeType = SHP_GetShapeType(eFileName);
    std::cerr << "eFileName=" << eFileName << " eShapeType=" << eShapeType
              << "\n";
    std::cerr << "Normal Termination of SHP_GetShapeType\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in SHP_GetShapeType\n";
    exit(e.eVal);
  }
  runtime(time1);
}
