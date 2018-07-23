//#include "Model_interpolation.h"

#include "Floats.h"

int main(int argc, char *argv[])
{
  try {
    FullNamelist eFull=NAMELIST_GetStandard_CREATE_LTransInput();
    if (argc != 2) {
      std::cerr << "CREATE_RomsOfflineInput is used as\n";
      std::cerr << "CREATE_RomsOfflineInput [file.nml]\n";
      std::cerr << "with file.nml the file describing the choices made\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull);
      return -1;
    }
    std::string eFileName=argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    SingleBlock eBlPROC=eFull.ListBlock.at("PROC");
    std::string eFile=eBlPROC.ListStringValues.at("FloatFile");
    std::string eModelName=eBlPROC.ListStringValues.at("MODELNAME");
    CHECK_Model_Allowedness(eModelName);
    std::string GridFile=eBlPROC.ListStringValues.at("GridFile");
    std::string HisPrefix=eBlPROC.ListStringValues.at("HisPrefix");
    TripleModelDesc eTriple{eModelName, GridFile, "unset", HisPrefix, {}};
    GridArray GrdArr=RETRIEVE_GRID_ARRAY(eTriple);
    double DistanceKM=eBlPROC.ListDoubleValues.at("MaxDistKM");
    std::vector<double> ListLONpt=eBlPROC.ListListDoubleValues.at("ListLonFloat");
    std::vector<double> ListLATpt=eBlPROC.ListListDoubleValues.at("ListLatFloat");
    std::vector<double> ListDepth=eBlPROC.ListListDoubleValues.at("ListDepth");
    std::vector<int> ListTime=eBlPROC.ListListIntValues.at("ListTime");
    ROMS_CreateLTransFile(eFile, GrdArr, DistanceKM, ListLONpt, ListLATpt, ListDepth, ListTime);
    std::cerr << "Normal termination of CREATE_RomsOfflineInput\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
