// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Model_interpolation.h"

int main(int argc, char *argv[]) {
  srand_random_set();
  SingletonTime time1;
  try {
    double pi = 3.1415926535;
    double g = 9.80665;
    if (argc != 2) {
      std::cerr << "DATA_CreateIdealCurrent [file.nc]\n";
      return -1;
    }
    std::string eFileName = argv[1];
    //
    GridArray GrdArr = WWM_ReadGridFile_netcdf(eFileName);
    MyMatrix<double> LON = GrdArr.GrdArrRho.LON;
    MyMatrix<double> LAT = GrdArr.GrdArrRho.LAT;
    int nx = LON.rows();
    int ny = LON.cols();
    double avgLon = LON.sum() / double(nx * ny);
    double avgLat = LAT.sum() / double(nx * ny);
    double maxLon = LON.maxCoeff();
    double minLon = LON.minCoeff();
    double maxLat = LAT.maxCoeff();
    double minLat = LAT.minCoeff();
    std::cerr << "LON(min/max/avg)=" << minLon << " , " << maxLon << " , "
              << avgLon << "\n";
    std::cerr << "LAT(min/max/avg)=" << minLat << " , " << maxLat << " , "
              << avgLat << "\n";
    //
    std::vector<double> ListTime = NC_ReadTimeFromFile(eFileName, "ocean_time");
    int nbTime = ListTime.size();
    //
    for (int iTime = 0; iTime < nbTime; iTime++) {
      MyMatrix<double> UsurfCurr =
          NETCDF_Get2DvariableSpecEntry(eFileName, GrdArr, "UsurfCurr", iTime);
      MyMatrix<double> VsurfCurr =
          NETCDF_Get2DvariableSpecEntry(eFileName, GrdArr, "VsurfCurr", iTime);
      MyMatrix<double> Hwave =
          NETCDF_Get2DvariableSpecEntry(eFileName, GrdArr, "HS", iTime);
      MyMatrix<double> TM02 =
          NETCDF_Get2DvariableSpecEntry(eFileName, GrdArr, "TM02", iTime);
      MyMatrix<double> MAPSTA =
          NETCDF_Get2DvariableSpecEntry(eFileName, GrdArr, "MAPSTA", iTime);
      MyMatrix<double> Cspeed(nx, ny);
      for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++) {
          double eT = TM02(i, j);
          double sigma = double(2) * pi / eT;
          double k = sigma * sigma / g;
          double lambda = double(2) * pi / k;
          double eCspeed = lambda / eT;
          Cspeed(i, j) = eCspeed;
        }
      MyMatrix<double> Unorm(nx, ny);
      for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++) {
          double eU = UsurfCurr(i, j);
          double eV = VsurfCurr(i, j);
          double eNorm = sqrt(eU * eU + eV * eV);
          Unorm(i, j) = eNorm;
        }
      bool IsFirst = true;
      double MinDist = -1; // just to avoid the warning.
      int iFound = -1, jFound = -1;
      for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
          if (MAPSTA(i, j) == 2) {
            double dLon = LON(i, j) - avgLon;
            double dLat = LAT(i, j) - avgLat;
            double eDist = dLon * dLon + dLat * dLat;
            if (IsFirst) {
              MinDist = eDist;
              iFound = i;
              jFound = j;
            } else {
              if (eDist < MinDist) {
                MinDist = eDist;
                iFound = i;
                jFound = j;
              }
            }
            IsFirst = false;
          }
      double c0 = Cspeed(iFound, jFound);
      double H0 = Hwave(iFound, jFound);
      MyMatrix<double> TM02ret(nx, ny);
      MyMatrix<double> HwaveRet(nx, ny);
      for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++) {
          double U = Unorm(i, j);
          double eQuot = (1 + sqrt(1 + 4 * U / c0)) / 2;
          double c = c0 * eQuot;
          double sigma = g / c;
          double eT = 2 * pi / sigma;
          TM02ret(i, j) = eT;
          double eQuot2 = c0 * c0 / (c * (c + 2 * U));
          double eH = H0 * sqrt(eQuot2);
          HwaveRet(i, j) = eH;
        }
      NETCDF_Write2DvariableSpecEntry(eFileName, "HS", iTime, HwaveRet);
      NETCDF_Write2DvariableSpecEntry(eFileName, "TM02", iTime, TM02ret);
    }
    std::cerr << "Normal termination of DATA_CreateIdealCurrent\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in DATA_CreateIdealCurrent\n";
    exit(e.eVal);
  }
  runtime(time1);
}
