#include "Plotting_fct.h"
#include "Basic_netcdf.h"

int main(int argc, char *argv[])
{
  std::cerr << std::fixed;
  std::cerr << std::setprecision(9);
  try {
    std::string ePrefix = "./P1/out_";
    std::string eFilePrev = ePrefix + "2016_6_1.nc";
    Eigen::Tensor<int,3> eTens = NC_Read3Dvariable_Mask_file(eFilePrev, "surf_el");
    MyMatrix<int> MSK = StrictProjectionMask(eTens, 0);
    for (int iMonth=6; iMonth<=12; iMonth++) {
      int monlen = MONTH_LEN(2016,iMonth);
      for (int iDay=1; iDay<=monlen; iDay++) {
	std::stringstream s;
	s << ePrefix << "2016_" << iMonth << "_" << iDay << ".nc";
	std::string eFile = s.str();
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
	std::cerr << "eFile=" << eFile << " nbError=" << nbError << "\n";
      }
    }

    
  }
  catch (TerminalException const& e) {
    std::cerr << "Error in PLOT_results\n";
    exit(e.eVal);
  }
}
