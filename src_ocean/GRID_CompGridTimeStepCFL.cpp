#include "Model_grids.h"
int main(int argc, char *argv[])
{
  srand_random_set();
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
  try {
    if (argc != 4) {
      std::cerr << "GRID_CompGridTimeStepCFL [ModelName] [GridFile] [HisPrefix]\n";
      std::cerr << "with ModelName     the chosen model\n";
      std::cerr << "with GridFile      the chosen grid\n";
      std::cerr << " and HisPrefix     the history prefix\n";
      return -1;
    }
    std::string ModelName   = argv[1];
    std::string GridFile    = argv[2];
    std::string HisPrefix   = argv[3];
    TripleModelDesc eTriple{ModelName, GridFile, "unset", HisPrefix, {}};
    GridArray GrdArr=PRE_RETRIEVE_GRID_ARRAY(eTriple);
    std::cerr << "We have GrdArr\n";
    DataCFL rec = ComputeTimeStepCFL(GrdArr);
    std::cerr << "TimeStepCFL = " << rec.MinTimeStep << " MinDist=" << rec.MinDist << " AvgDist=" << rec.AvgDist << "\n";
    std::cerr << "Normal termination of GRID_CompGridTimeStepCFL\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Error in GRID_CompGridTimeStepCFL\n";
    exit(e.eVal);
  }
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "runtime = " << std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count() << "\n";
}
