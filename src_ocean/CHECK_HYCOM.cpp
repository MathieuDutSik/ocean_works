// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Basic_netcdf.h"
#include "Plotting_fct.h"

int main(int argc, char *argv[]) {
  srand_random_set();
  std::cerr << std::fixed;
  std::cerr << std::setprecision(9);
  HumanTime time1;
  try {
    if (argc != 2) {
      std::cerr << "CHECK_HYCOM [ePrefix]\n";
      return -1;
    }
    std::string ePrefix = argv[1];
    std::vector<std::string> ListFile =
        FILE_DirectoryMatchingPrefixExtension(ePrefix, ".nc");
    if (ListFile.size() == 0) {
      std::cerr << "We found 0 files in ePrefix=" << ePrefix << "\n";
      throw TerminalException{1};
    }
    std::string eFilePrev = ListFile[0];
    Eigen::Tensor<int, 3> eTens =
        NC_Read3Dvariable_Mask_file(eFilePrev, "surf_el");
    MyMatrix<int> MSK = StrictProjectionMask(eTens, 0);
    int TotalNbError = 0;
    int TotalNbErrorMSK = 0;
    std::vector<int> ListNbTime;
    for (auto &eFile : ListFile) {
      //
      // Checking if masks changes as things move on
      //
      Eigen::Tensor<int, 3> eTens_B =
          NC_Read3Dvariable_Mask_file(eFilePrev, "surf_el");
      MyMatrix<int> MSK_B = StrictProjectionMask(eTens_B, 0);
      MyMatrix<int> MSK_diff = MSK - MSK_B;
      int nbRow = MSK_diff.rows();
      int nbCol = MSK_diff.cols();
      int nbErrorMSK = 0;
      for (int iRow = 0; iRow < nbRow; iRow++)
        for (int iCol = 0; iCol < nbCol; iCol++)
          if (MSK_diff(iRow, iCol) != 0)
            nbErrorMSK++;
      //
      // Checking that salinity mask corresponds to the ones of surface eleva
      //
      netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
      netCDF::NcVar data = dataFile.getVar("salinity");
      MyVector<int> StatusFill = NC_ReadVariable_StatusFill_data<int>(data);
      std::vector<size_t> ListDim = NC_ReadVariable_listdim(data);
      int nbTime = ListDim[0];
      int s_vert = ListDim[1];
      int eta = ListDim[2];
      int xi = ListDim[3];
      int idx = 0;
      int nbError = 0;
      ListNbTime.push_back(nbTime);
      for (int iTime = 0; iTime < nbTime; iTime++)
        for (int iS = 0; iS < s_vert; iS++)
          for (int i = 0; i < eta; i++)
            for (int j = 0; j < xi; j++) {
              if (MSK(i, j) == 1) {
                if (StatusFill(idx) != 1) {
                  nbError++;
                }
              }
              idx++;
            }
      //
      std::cerr << "eFile=" << eFile << " nbError=" << nbError
                << " nbErrorMSK=" << nbErrorMSK << " nbTime=" << nbTime << "\n";
      TotalNbError += nbError;
      TotalNbErrorMSK += nbErrorMSK;
    }
    std::cerr << "TotalNbError=" << TotalNbError
              << "  TotalNbErrorMSK=" << TotalNbErrorMSK << "\n";
    CollectedResult<int> eColl = Collected(ListNbTime);
    int nbEnt = eColl.LVal.size();
    for (int iEnt = 0; iEnt < nbEnt; iEnt++) {
      std::cerr << "nbTime=" << eColl.LVal[iEnt] << " attained "
                << eColl.LMult[iEnt] << " times\n";
    }
    std::cerr << "Normal termination of CHECK_HYCOM\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in CHECK_HYCOM\n";
    exit(e.eVal);
  }
  runtime(time1);
}
