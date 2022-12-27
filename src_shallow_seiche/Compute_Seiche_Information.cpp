#include "seiche_numeric_eigenvalue.h"
int main(int argc, char *argv[]) {
  std::cerr << std::fixed;
  std::cerr << std::setprecision(9);
  try {
    FullNamelist eFull = NAMELIST_SEICHE_Eigen();
    if (argc != 2) {
      std::cerr << "Compute_Seiche_Information is used as\n";
      std::cerr << "Compute_Seiche_Information [file.nml]\n";
      std::cerr << "with file.nml the file describing the plotting routines\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    //
    SingleBlock eBlCOMP = eFull.ListBlock.at("COMP");
    // Loading grid
    std::string GridFile = eBlCOMP.ListStringValues.at("GridFile");
    std::string BoundFile = "unset";
    GridArray GrdArr = ReadUnstructuredGrid(GridFile, BoundFile);
    // Loading numerical parameters
    double h0 = eBlCOMP.ListDoubleValues.at("h0");
    int maxNbEig = eBlCOMP.ListIntValues.at("maxNbEig");
    // Loading the OutFile
    std::string OutFile = eBlCOMP.ListStringValues.at("OutFile");
    bool RescaleEigenvector = eBlCOMP.ListBoolValues.at("RescaleEigenvector");
    bool UseSI_Bmatrix = eBlCOMP.ListBoolValues.at("UseSI_Bmatrix");
    int ncv = eBlCOMP.ListIntValues.at("ncv");
    //
    // Now computing
    //
    std::vector<PeriodicSolution> ListSol = ComputeEigenvaluesSWE1(
        h0, maxNbEig, GrdArr, RescaleEigenvector, UseSI_Bmatrix, ncv);
    //
    // Outputing it
    //
    WriteSeicheInfoAsNetcdfFile(ListSol, OutFile);
    std::cerr << "Normal termination of Compute_Seiche_Information\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in Compute_Seiche_Information\n";
    exit(e.eVal);
  }
}
