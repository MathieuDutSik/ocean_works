//#include "NamelistExampleOcean.h"
#include "Plotting_fct.h"

int main(int argc, char *argv[])
{
  srand_random_set();
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
  try {
    FullNamelist eFull = NAMELIST_GetStandard_PlotGrid();
    if (argc != 2) {
      std::cerr << "PLOT_grid [file.nml]\n";
      std::cerr << "with file.nml the file describing the plotting routines\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull);
      return -1;
    }
    std::string eFileName=argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    //
    TripleModelDesc eTriple = Retrieve_triple_from_array(eFull);
    GridArray GrdArr=RETRIEVE_GRID_ARRAY(eTriple);
    PermanentInfoDrawing ePerm=GET_PERMANENT_INFO(eFull);
    ePerm.eDrawArr = CommonAssignation_DrawArr(eFull);
    std::cerr << "ePerm.eDrawArr=\n";
    PrintDrawArray(std::cerr, ePerm.eDrawArr);
    NCLcaller<GeneralType> eCall(ePerm.NPROC);
    GRID_PLOTTING(GrdArr, eTriple.GridFile, eCall, ePerm);
    //
    std::cerr << "Normal termination of PLOT_grid\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Error in PLOT_grid\n";
    exit(e.eVal);
  }
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "runtime = " << std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count() << "\n";
}
