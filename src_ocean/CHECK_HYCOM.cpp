#include "Plotting_fct.h"
#include "Basic_netcdf.h"

int main(int argc, char *argv[])
{
  std::cerr << std::fixed;
  std::cerr << std::setprecision(9);
  try {
    if (argc != 2) {
      std::cerr << "CHECK_HYCOM is used as\n";
      std::cerr << "CHECK_HYCOM [ePrefix]\n";
      throw TerminalException{1};
    }
    std::string ePrefix = argv[1];
    std::vector<std::string> ListFile=FILE_DirectoryMatchingPrefixExtension(ePrefix, ".nc");
    if (ListFile.size() == 0) {
      std::cerr << "We found 0 files in ePrefix=" << ePrefix << "\n";
      throw TerminalException{1};
    }
    std::string eFilePrev = ListFile[0];
    Eigen::Tensor<int,3> eTens = NC_Read3Dvariable_Mask_file(eFilePrev, "surf_el");
    MyMatrix<int> MSK = StrictProjectionMask(eTens, 0);
    for (auto &eFile : ListFile) {
      //
      // Checking if masks changes as things move on
      //
      Eigen::Tensor<int,3> eTens_B = NC_Read3Dvariable_Mask_file(eFilePrev, "surf_el");
      MyMatrix<int> MSK_B = StrictProjectionMask(eTens_B, 0);
      MyMatrix<int> MSK_diff = MSK - MSK_B;
      int nbRow=MSK_diff.rows();
      int nbCol=MSK_diff.cols();
      int nbErrorMSK=0;
      for (int iRow=0; iRow<nbRow; iRow++)
	for (int iCol=0; iCol<nbCol; iCol++)
	  if (MSK_diff(iRow, iCol) != 0)
	    nbErrorMSK++;
      //
      // Checking that salinity mask corresponds to the ones of surface eleva
      //
      netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
      netCDF::NcVar data=dataFile.getVar("salinity");
      MyVector<int> StatusFill=NC_ReadVariable_StatusFill_data(data);
      std::vector<size_t> ListDim = NC_ReadVariable_listdim(data);
      int nbTime=ListDim[0];
      int s_vert=ListDim[1];
      int eta=ListDim[2];
      int xi=ListDim[3];
      int idx=0;
      int nbError=0;
      for (int iTime=0; iTime<nbTime; iTime++)
	for (int iS=0; iS<s_vert; iS++)
	  for (int i=0; i<eta; i++)
	    for (int j=0; j<xi; j++) {
	      if (MSK(i,j) == 1) {
		if (StatusFill(idx) != 1) {
		  //		    std::cerr << "Found inconsistency at i=" << i << " j=" << j << " iTime=" << iTime << " iS=" << iS << "\n";
		  nbError++;
		}
	      }
	      idx++;
	    }
      //
      std::cerr << "eFile=" << eFile << " nbError=" << nbError << " nbErrorMSK=" << nbErrorMSK << "\n";
    }
  }
  catch (TerminalException const& e) {
    std::cerr << "Error in PLOT_results\n";
    exit(e.eVal);
  }
}
