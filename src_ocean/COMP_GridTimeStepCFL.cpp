#include "Model_grids.h"
int main(int argc, char *argv[])
{
  try {
    if (argc != 4) {
      std::cerr << "COMP_GridTimeStepCFL is used as\n";
      std::cerr << "COMP_GridTimeStepCFL [ModelName] [GridFile] [HisPrefix]\n";
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
    double TimeStepCFL = ComputeTimeStepCFL(GrdArr);
    std::cerr << "TimeStepCFL = " << TimeStepCFL << "\n";
    std::cerr << "Normal termination of ConvertGrid\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Error in ConvertGrid\n";
    exit(e.eVal);
  }
}
