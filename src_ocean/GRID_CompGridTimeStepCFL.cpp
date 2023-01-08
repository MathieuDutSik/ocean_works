// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Model_grids.h"
int main(int argc, char *argv[]) {
  srand_random_set();
  SingletonTime time1;
  try {
    if (argc != 4) {
      std::cerr
          << "GRID_CompGridTimeStepCFL [ModelName] [GridFile] [HisPrefix]\n";
      std::cerr << "with ModelName     the chosen model\n";
      std::cerr << "with GridFile      the chosen grid\n";
      std::cerr << " and HisPrefix     the history prefix\n";
      return -1;
    }
    std::string ModelName = argv[1];
    std::string GridFile = argv[2];
    std::string HisPrefix = argv[3];
    TripleModelDesc eTriple{ModelName, GridFile, "unset", HisPrefix, {}};
    GridArray GrdArr = PRE_RETRIEVE_GRID_ARRAY(eTriple);
    std::cerr << "We have GrdArr\n";
    DataCFL rec = ComputeTimeStepCFL(GrdArr);
    std::cerr << "TimeStepCFL = " << rec.MinTimeStep
              << " MinDist=" << rec.MinDist << " AvgDist=" << rec.AvgDist
              << "\n";
    //    double MaxDist = ComputeMaxDistance(GrdArr);
    //    std::cerr << "MaxDist=" << MaxDist << "\n";
    QuadArray Qarr = GetQuadArray(GrdArr);
    std::cerr << "MinLon=" << Qarr.MinLon << "\n";
    std::cerr << "MaxLon=" << Qarr.MaxLon << "\n";
    std::cerr << "MinLat=" << Qarr.MinLat << "\n";
    std::cerr << "MaxLat=" << Qarr.MaxLat << "\n";
    std::cerr << "Normal termination of GRID_CompGridTimeStepCFL\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRID_CompGridTimeStepCFL\n";
    exit(e.eVal);
  }
  runtime(time1);
}
