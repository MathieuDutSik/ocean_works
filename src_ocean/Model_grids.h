#ifndef MODEL_GRIDS_INCLUDE
#define MODEL_GRIDS_INCLUDE

#include "Basic_grib.h"
#include "Basic_netcdf.h"
#include "CommonFuncModel.h"
#include "SphericalGeom.h"
#include "Triangulations.h"
#include <limits>
#include <algorithm>
#include <map>
#include <set>
#include <unordered_map>
#include <vector>
#include <string>
#include <utility>

ARVDtyp GetTrivialARrayVerticalDescription() {
  ARVDtyp ARVD;
  ARVD.IsAssigned = false;
  ARVD.Zcoordinate = false;
  ARVD.ModelName = "UNSET";
  return ARVD;
}

ARVDtyp ReadROMSverticalStratification(std::string const &eFile) {
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
  if (dataFile.isNull()) {
    std::cerr << "Error while Netcdf opening of file=" << eFile << "\n";
    throw TerminalException{1};
  }
  //
  netCDF::NcVar data_Vtrans = dataFile.getVar("Vtransform");
  if (data_Vtrans.isNull()) {
    std::cerr << "Error while opening variable Vtransform\n";
    throw TerminalException{1};
  }
  //
  netCDF::NcVar data_Vstret = dataFile.getVar("Vstretching");
  if (data_Vstret.isNull()) {
    std::cerr << "Error while opening variable Vstretching\n";
    throw TerminalException{1};
  }
  //
  netCDF::NcVar data_theta_s = dataFile.getVar("theta_s");
  if (data_theta_s.isNull()) {
    std::cerr << "Error while opening variable theta_s\n";
    throw TerminalException{1};
  }
  //
  netCDF::NcVar data_theta_b = dataFile.getVar("theta_b");
  if (data_theta_b.isNull()) {
    std::cerr << "Error while opening variable theta_b\n";
    throw TerminalException{1};
  }
  //
  netCDF::NcVar data_Tcline = dataFile.getVar("Tcline");
  if (data_Tcline.isNull()) {
    std::cerr << "Error while opening variable Tcline\n";
    throw TerminalException{1};
  }
  //
  netCDF::NcVar data_hc = dataFile.getVar("hc");
  if (data_hc.isNull()) {
    std::cerr << "Error while opening variable hc\n";
    throw TerminalException{1};
  }
  //
  std::vector<int> eValI(1);
  std::vector<double> eValD(1);
  //
  ARVDtyp ARVD;
  ARVD.IsAssigned = true;
  ARVD.ModelName = "ROMS";
  data_Vtrans.getVar(eValI.data());
  ARVD.Vtransform = eValI[0];
  //
  data_Vstret.getVar(eValI.data());
  ARVD.Vstretching = eValI[0];
  //
  data_theta_s.getVar(eValD.data());
  ARVD.theta_s = eValD[0];
  //
  data_theta_b.getVar(eValD.data());
  ARVD.theta_b = eValD[0];
  //
  data_Tcline.getVar(eValD.data());
  ARVD.Tcline = eValD[0];
  //
  data_hc.getVar(eValD.data());
  ARVD.hc = eValD[0];
  //
  netCDF::NcVar data_Cs_r = dataFile.getVar("Cs_r");
  if (data_Cs_r.isNull()) {
    std::cerr << "Error while opening variable Cs_r\n";
    throw TerminalException{1};
  }
  netCDF::NcVar data_Cs_w = dataFile.getVar("Cs_w");
  if (data_Cs_w.isNull()) {
    std::cerr << "Error while opening variable Cs_w\n";
    throw TerminalException{1};
  }
  netCDF::NcVar data_s_r = dataFile.getVar("s_rho");
  if (data_s_r.isNull()) {
    std::cerr << "Error while opening variable s_rho\n";
    throw TerminalException{1};
  }
  netCDF::NcVar data_s_w = dataFile.getVar("s_w");
  if (data_s_w.isNull()) {
    std::cerr << "Error while opening variable s_w\n";
    throw TerminalException{1};
  }
  netCDF::NcDim eDim = data_Cs_r.getDim(0);
  int dim_s_r = eDim.getSize();
  netCDF::NcDim fDim = data_Cs_w.getDim(0);
  int dim_s_w = fDim.getSize();
  std::vector<double> eVarR(dim_s_r), eVarW(dim_s_w);
  data_Cs_r.getVar(eVarR.data());
  MyVector<double> V_r(dim_s_r), V_w(dim_s_w);
  for (int i = 0; i < dim_s_r; i++)
    V_r(i) = eVarR[i];
  ARVD.Cs_r = V_r;
  data_s_r.getVar(eVarR.data());
  for (int i = 0; i < dim_s_r; i++)
    V_r(i) = eVarR[i];
  ARVD.sc_r = V_r;
  data_Cs_w.getVar(eVarW.data());
  for (int i = 0; i < dim_s_w; i++)
    V_w(i) = eVarW[i];
  ARVD.Cs_w = V_w;
  data_s_w.getVar(eVarW.data());
  for (int i = 0; i < dim_s_w; i++)
    V_w(i) = eVarW[i];
  ARVD.sc_w = V_w;
  ARVD.N = dim_s_r;
  ARVD.Zcoordinate = false;
  //
  return ARVD;
}

void WriteROMSverticalStratification(netCDF::NcFile &dataFile,
                                     ARVDtyp const &ARVD) {
  if (dataFile.isNull()) {
    std::cerr << "WriteROMSverticalStratification error, dataFile is null\n";
    throw TerminalException{1};
  }
  //
  // First the scalar values
  //
  std::vector<int> eValI(1);
  std::vector<double> eValD(1);
  std::vector<std::string> EmptyListVar;
  //
  netCDF::NcVar data_Vtrans =
      dataFile.addVar("Vtransform", "int", EmptyListVar);
  if (data_Vtrans.isNull()) {
    std::cerr << "Error while opening variable Vtransform\n";
    throw TerminalException{1};
  }
  eValI[0] = ARVD.Vtransform;
  data_Vtrans.putVar(eValI.data());
  //
  netCDF::NcVar data_Vstret =
      dataFile.addVar("Vstretching", "int", EmptyListVar);
  if (data_Vstret.isNull()) {
    std::cerr << "Error while opening variable Vstretching\n";
    throw TerminalException{1};
  }
  eValI[0] = ARVD.Vstretching;
  data_Vstret.putVar(eValI.data());
  //
  netCDF::NcVar data_theta_s =
      dataFile.addVar("theta_s", "double", EmptyListVar);
  if (data_theta_s.isNull()) {
    std::cerr << "Error while opening variable theta_s\n";
    throw TerminalException{1};
  }
  eValD[0] = ARVD.theta_s;
  data_theta_s.putVar(eValD.data());
  //
  netCDF::NcVar data_theta_b =
      dataFile.addVar("theta_b", "double", EmptyListVar);
  if (data_theta_b.isNull()) {
    std::cerr << "Error while opening variable theta_b\n";
    throw TerminalException{1};
  }
  eValD[0] = ARVD.theta_b;
  data_theta_b.putVar(eValD.data());
  //
  netCDF::NcVar data_Tcline = dataFile.addVar("Tcline", "double", EmptyListVar);
  if (data_Tcline.isNull()) {
    std::cerr << "Error while opening variable Tcline\n";
    throw TerminalException{1};
  }
  eValD[0] = ARVD.Tcline;
  data_Tcline.putVar(eValD.data());
  //
  netCDF::NcVar data_hc = dataFile.addVar("hc", "double", EmptyListVar);
  if (data_hc.isNull()) {
    std::cerr << "Error while opening variable hc\n";
    throw TerminalException{1};
  }
  eValD[0] = ARVD.hc;
  data_hc.putVar(eValD.data());
  //
  // Now the arrays
  //
  int s_rho = ARVD.Cs_r.size();
  int s_w = ARVD.Cs_w.size();
  std::vector<double> eVarR(s_rho), eVarW(s_w);
  std::string strSRho = "s_rho";
  std::string strSW = "s_w";
  netCDF::NcVar data_Cs_r = dataFile.addVar("Cs_r", "double", {strSRho});
  for (int i = 0; i < s_rho; i++)
    eVarR[i] = ARVD.Cs_r(i);
  data_Cs_r.putVar(eVarR.data());
  //
  netCDF::NcVar data_Cs_w = dataFile.addVar("Cs_w", "double", {strSW});
  for (int i = 0; i < s_w; i++)
    eVarW[i] = ARVD.Cs_w(i);
  data_Cs_w.putVar(eVarW.data());
  //
  netCDF::NcVar data_s_r = dataFile.addVar("s_rho", "double", {strSRho});
  for (int i = 0; i < s_rho; i++)
    eVarR[i] = ARVD.sc_r(i);
  data_s_r.putVar(eVarR.data());
  //
  netCDF::NcVar data_s_w = dataFile.addVar("s_w", "double", {strSW});
  for (int i = 0; i < s_w; i++)
    eVarW[i] = ARVD.sc_w(i);
  data_s_w.putVar(eVarW.data());
}

//
// This code is adapted from set_scoord.F of ROMS Rutgers
//
ARVDtyp ROMSgetARrayVerticalDescription(int const &N, int const &Vtransform,
                                        int const &Vstretching,
                                        double const &Tcline, double const &hc,
                                        double const &theta_s,
                                        double const &theta_b) {
  ARVDtyp ARVD;
  ARVD.IsAssigned = true;
  ARVD.ModelName = "ROMS";
  ARVD.Zcoordinate = false;
  ARVD.N = N;
  ARVD.Vtransform = Vtransform;
  ARVD.Vstretching = Vstretching;
  ARVD.Tcline = Tcline;
  ARVD.hc = hc;
  ARVD.theta_s = theta_s;
  ARVD.theta_b = theta_b;
  ARVD.Cs_r = ZeroVector<double>(N);
  ARVD.Cs_w = ZeroVector<double>(N + 1);
  ARVD.sc_r = ZeroVector<double>(N);
  ARVD.sc_w = ZeroVector<double>(N + 1);
  double half = double(1) / double(2);
  double ds = 1 / double(N);
  if (Vstretching == 1) {
    double cff1, cff2;
    if (theta_s > 0) {
      cff1 = 1 / sinh(theta_s);
      cff2 = 1 / (2 * tanh(theta_s / 2));
    } else {
      cff1 = -400; // just to avoid the warning.
      cff2 = -400; // just to avoid the warning.
    }
    ARVD.sc_w(0) = -1;
    ARVD.Cs_w(0) = -1;
    for (int k = 1; k <= N; k++) {
      double eSc_w = ds * double(k - N);
      double eSc_r = ds * (double(k - N) - half);
      ARVD.sc_w(k) = eSc_w;
      ARVD.sc_r(k - 1) = eSc_r;
      if (theta_s > 0) {
        ARVD.Cs_w(k) = (1 - theta_b) * cff1 * sinh(theta_s * eSc_w) +
                       theta_b * (cff2 * tanh(theta_b * (eSc_w + half)) - half);
        ARVD.Cs_r(k - 1) =
            (1 - theta_b) * cff1 * sinh(theta_s * eSc_r) +
            theta_b * (cff2 * tanh(theta_b * (eSc_r + half)) - half);
      } else {
        ARVD.Cs_w(k) = eSc_w;
        ARVD.Cs_r(k - 1) = eSc_r;
      }
    }
  }
  if (Vstretching == 2) {
    double Aweight = 1;
    double Bweight = 1;
    ARVD.sc_w(N) = 0;
    ARVD.Cs_w(N) = 0;
    for (int k = 1; k <= N - 1; k++) {
      double sc_w = ds * double(k - N);
      ARVD.sc_w(k) = sc_w;
      if (theta_s > 0) {
        double Csur = (1 - cosh(theta_s * sc_w)) / (cosh(theta_s) - 1);
        if (theta_b > 0) {
          double Cbot = sinh(theta_b * (sc_w + 1)) / sinh(theta_b) - 1;
          double Cweight =
              pow(sc_w + 1, Aweight) *
              (1 + (Aweight / Bweight) * (1 - pow(sc_w + 1, Bweight)));
          ARVD.Cs_w(k) = Cweight * Csur + (1 - Cweight) * Cbot;
        } else {
          ARVD.Cs_w(k) = Csur;
        }
      } else {
        ARVD.Cs_w(k) = sc_w;
      }
    }
    ARVD.sc_w(0) = -1;
    ARVD.Cs_w(0) = -1;
    for (int k = 1; k <= N; k++) {
      double sc_r = ds * (double(k - N) - half);
      ARVD.sc_r(k - 1) = sc_r;
      if (theta_s > 0) {
        double Csur = (1 - cosh(theta_s * sc_r)) / (cosh(theta_s) - 1);
        if (theta_b > 0) {
          double Cbot = sinh(theta_s * (sc_r + 1)) / sinh(theta_s) - 1;
          double Cweight =
              pow(sc_r + 1, Aweight) *
              (1 + (Aweight / Bweight) * (1 - pow(sc_r + 1, Bweight)));
          ARVD.Cs_r(k - 1) = Cweight * Csur + (1 - Cweight) * Cbot;
        } else {
          ARVD.Cs_r(k - 1) = Csur;
        }
      } else {
        ARVD.Cs_r(k - 1) = sc_r;
      }
    }
  }
  if (Vstretching == 3) {
    double exp_sur = theta_s;
    double exp_bot = theta_b;
    double Hscale = 3;
    ARVD.sc_w(N) = 0;
    ARVD.Cs_w(N) = 0;
    for (int k = 1; k <= N - 1; k++) {
      double sc_w = ds * double(k - N);
      ARVD.sc_w(k) = sc_w;
      double Cbot =
          log(cosh(Hscale * pow(sc_w + 1, exp_bot))) / log(cosh(Hscale)) - 1;
      double Csur =
          -log(cosh(Hscale * pow(abs(sc_w), exp_sur))) / log(cosh(Hscale));
      double Cweight = half * (1 - tanh(Hscale * (sc_w + half)));
      ARVD.Cs_w(k) = Cweight * Cbot + (1 - Cweight) * Csur;
    }
    ARVD.sc_w(0) = -1;
    ARVD.Cs_w(0) = -1;
    for (int k = 1; k <= N; k++) {
      double sc_r = ds * (double(k - N) - half);
      ARVD.sc_r(k - 1) = sc_r;
      double Cbot =
          log(cosh(Hscale * pow(sc_r + 1, exp_bot))) / log(cosh(Hscale)) - 1;
      double Csur =
          -log(cosh(Hscale * pow(abs(sc_r), exp_sur))) / log(cosh(Hscale));
      double Cweight = half * (1 - tanh(Hscale * (sc_r + half)));
      ARVD.Cs_r(k - 1) = Cweight * Cbot + (1 - Cweight) * Csur;
    }
  }
  if (Vstretching == 4) {
    ARVD.sc_w(N) = 0;
    ARVD.Cs_w(N) = 0;
    for (int k = 1; k <= N - 1; k++) {
      double sc_w = ds * double(k - N);
      ARVD.sc_w(k) = sc_w;
      double Csur;
      if (theta_s > 0) {
        Csur = (1 - cosh(theta_s * sc_w)) / (cosh(theta_s) - 1);
      } else {
        Csur = -pow(sc_w, 2);
      }
      if (theta_b > 0) {
        double Cbot = (exp(theta_b * Csur) - 1) / (1 - exp(-theta_b));
        ARVD.Cs_w(k) = Cbot;
      } else {
        ARVD.Cs_w(k) = Csur;
      }
    }
    ARVD.sc_w(0) = -1;
    ARVD.Cs_w(0) = -1;
    for (int k = 1; k <= N; k++) {
      double sc_r = ds * (double(k - N) - half);
      ARVD.sc_r(k - 1) = sc_r;
      double Csur;
      if (theta_s > 0)
        Csur = (1 - cosh(theta_s * sc_r)) / (cosh(theta_s) - 1);
      else
        Csur = -pow(sc_r, 2);
      if (theta_b > 0) {
        double Cbot = (exp(theta_b * Csur) - 1) / (1 - exp(-theta_b));
        ARVD.Cs_r(k - 1) = Cbot;
      } else {
        ARVD.Cs_r(k - 1) = Csur;
      }
    }
  }
  return ARVD;
}

QuadArray GetQuadArray(GridArray const &GrdArr) {
  double MinLon = 0, MaxLon = 0, MinLat = 0, MaxLat = 0;
  if (GrdArr.IsFE == 1) {
    int siz1 = GrdArr.GrdArrRho.LON.rows();
    int siz2 = GrdArr.GrdArrRho.LON.cols();
    if (siz1 == 0 || siz2 == 0) {
      std::cerr << "We need to have a nontrivial matrix\n";
      std::cerr << "siz1=" << siz1 << " siz2=" << siz2 << "\n";
      throw TerminalException{1};
    }
    MinLon = GrdArr.GrdArrRho.LON.minCoeff();
    MaxLon = GrdArr.GrdArrRho.LON.maxCoeff();
    MinLat = GrdArr.GrdArrRho.LAT.minCoeff();
    MaxLat = GrdArr.GrdArrRho.LAT.maxCoeff();
  } else {
    bool IsFirst = true;
    int eta_rho = GrdArr.GrdArrRho.LON.rows();
    int xi_rho = GrdArr.GrdArrRho.LON.cols();
    int eta_rho_msk = GrdArr.GrdArrRho.MSK.rows();
    int xi_rho_msk = GrdArr.GrdArrRho.MSK.cols();
    if (eta_rho_msk != eta_rho || xi_rho_msk != xi_rho) {
      std::cerr << "eta_rho_msk=" << eta_rho_msk << " xi_rho_msk=" << xi_rho_msk
                << "\n";
      std::cerr << "eta_rho    =" << eta_rho << " xi_rho    =" << xi_rho
                << "\n";
      std::cerr << "Dimension error in the arrays\n";
      throw TerminalException{1};
    }
    for (int i = 0; i < eta_rho; i++)
      for (int j = 0; j < xi_rho; j++)
        if (GrdArr.GrdArrRho.MSK(i, j) == 1) {
          double eLon = GrdArr.GrdArrRho.LON(i, j);
          double eLat = GrdArr.GrdArrRho.LAT(i, j);
          if (IsFirst) {
            MinLon = eLon;
            MaxLon = eLon;
            MinLat = eLat;
            MaxLat = eLat;
            IsFirst = false;
          } else {
            if (eLon < MinLon)
              MinLon = eLon;
            if (eLon > MaxLon)
              MaxLon = eLon;
            if (eLat < MinLat)
              MinLat = eLat;
            if (eLat > MaxLat)
              MaxLat = eLat;
          }
        }
  }
  return {MinLon, MaxLon, MinLat, MaxLat};
}

struct DataCFL {
  double MinTimeStep;
  double AvgDist;
  double MinDist;
};

DataCFL ComputeTimeStepCFL(GridArray const &GrdArr) {
  double miss_val = std::numeric_limits<double>::max();
  double ConstantGravity = 9.81;
  auto CompDist = [&](double const &eX, double const &eY, double const &fX,
                      double const &fY) -> double {
    if (GrdArr.IsSpherical) {
      return 1000 * GeodesicDistanceKM(eX, eY, fX, fY);
    }
    double deltaX = eX - fX;
    double deltaY = eY - fY;
    return sqrt(deltaX * deltaX + deltaY * deltaY);
  };
  double MinTimeStep = miss_val;
  double MinDist = miss_val;
  double sumDist = 0;
  size_t nb = 0;
  std::cerr << "IsFE=" << GrdArr.IsFE << "\n";
  const MyMatrix<double> &DEP = GetDEP(GrdArr.GrdArrRho);
  if (GrdArr.IsFE) {
    int mnp = DEP.rows();
    int mne = GrdArr.INE.rows();
    std::vector<double> ListMinDist(mnp, miss_val);
    for (int ie = 0; ie < mne; ie++) {
      for (int i = 0; i < 3; i++) {
        for (int u = 0; u < 2; u++) {
          //          std::cerr << "ie=" << ie << " i=" << i << " u=" << u <<
          //          "\n";
          int ushift = 2 * u + 2;
          int j = (i + ushift) % 3;
          //          std::cerr << "j=" << j << "\n";
          int IP = GrdArr.INE(ie, i);
          int JP = GrdArr.INE(ie, j);
          double eX = GrdArr.GrdArrRho.LON(IP, 0);
          double eY = GrdArr.GrdArrRho.LAT(IP, 0);
          double fX = GrdArr.GrdArrRho.LON(JP, 0);
          double fY = GrdArr.GrdArrRho.LAT(JP, 0);
          double eDist = CompDist(eX, eY, fX, fY);
          if (ListMinDist[IP] > eDist)
            ListMinDist[IP] = eDist;
          sumDist += eDist;
          nb++;
        }
      }
    }
    std::cerr << "We have ListMinDist\n";
    for (int ip = 0; ip < mnp; ip++) {
      double eTimeStep = ListMinDist[ip] / sqrt(ConstantGravity * DEP(ip, 0));
      if (MinTimeStep > eTimeStep)
        MinTimeStep = eTimeStep;
      if (MinDist > ListMinDist[ip])
        MinDist = ListMinDist[ip];
    }
  } else {
    // the other case
    int eta_rho = DEP.rows();
    int xi_rho = DEP.cols();
    std::vector<std::vector<int>> LNeigh{{1, 0}, {0, 1}, {-1, 0}, {0, -1}};
    for (int iEta = 0; iEta < eta_rho; iEta++)
      for (int iXi = 0; iXi < xi_rho; iXi++)
        if (GrdArr.GrdArrRho.MSK(iEta, iXi)) {
          double LocMinDist = miss_val;
          for (auto &eNeigh : LNeigh) {
            int iEtaN = iEta + eNeigh[0];
            int iXiN = iXi + eNeigh[1];
            if (iEtaN >= 0 && iEtaN < eta_rho && iXiN >= 0 && iXiN < xi_rho) {
              double eX = GrdArr.GrdArrRho.LON(iEta, iXi);
              double eY = GrdArr.GrdArrRho.LAT(iEta, iXi);
              double fX = GrdArr.GrdArrRho.LON(iEtaN, iXiN);
              double fY = GrdArr.GrdArrRho.LAT(iEtaN, iXiN);
              double eDist = CompDist(eX, eY, fX, fY);
              sumDist += eDist;
              nb++;
              if (LocMinDist > eDist)
                LocMinDist = eDist;
            }
          }
          //
          if (LocMinDist != miss_val) {
            double eDEP = DEP(iEta, iXi);
            double eTimeStep = MinDist / sqrt(ConstantGravity * eDEP);
            std::cerr << "MinDist=" << MinDist << " eDEP=" << eDEP
                      << " TimeStep=" << eTimeStep << "\n";
            if (MinTimeStep > eTimeStep) {
              MinTimeStep = eTimeStep;
            }
            if (MinDist > LocMinDist) {
              MinDist = LocMinDist;
            }
          }
        }
  }
  double AvgDist = sumDist / nb;
  return {MinTimeStep, AvgDist, MinDist};
}

double ComputeMaxDistance(GridArray const &GrdArr) {
  auto CompDist = [&](double const &eX, double const &eY, double const &fX,
                      double const &fY) -> double {
    if (GrdArr.IsSpherical) {
      return 1000 * GeodesicDistanceKM(eX, eY, fX, fY);
    }
    double deltaX = eX - fX;
    double deltaY = eY - fY;
    return sqrt(deltaX * deltaX + deltaY * deltaY);
  };
  double MaxDist = 0;
  const MyMatrix<double> &DEP = GetDEP(GrdArr.GrdArrRho);
  if (GrdArr.IsFE) {
    int mnp = DEP.rows();
    for (int ip = 0; ip < mnp; ip++) {
      double eX = GrdArr.GrdArrRho.LON(ip, 0);
      double eY = GrdArr.GrdArrRho.LAT(ip, 0);
      for (int jp = ip + 1; jp < mnp; jp++) {
        double fX = GrdArr.GrdArrRho.LON(jp, 0);
        double fY = GrdArr.GrdArrRho.LAT(jp, 0);
        double eDist = CompDist(eX, eY, fX, fY);
        if (eDist > MaxDist)
          MaxDist = eDist;
      }
    }
  } else {
    // the other case
    int eta_rho = DEP.rows();
    int xi_rho = DEP.cols();
    std::vector<std::pair<int, int>> LPair;
    for (int iEta = 0; iEta < eta_rho; iEta++)
      for (int iXi = 0; iXi < xi_rho; iXi++)
        if (GrdArr.GrdArrRho.MSK(iEta, iXi))
          LPair.push_back({iEta, iXi});
    size_t len = LPair.size();
    for (size_t i = 0; i < len; i++) {
      int iEta = LPair[i].first;
      int iXi = LPair[i].second;
      double eX = GrdArr.GrdArrRho.LON(iEta, iXi);
      double eY = GrdArr.GrdArrRho.LAT(iEta, iXi);
      for (size_t j = i + 1; j < len; j++) {
        int jEta = LPair[j].first;
        int jXi = LPair[j].second;
        double fX = GrdArr.GrdArrRho.LON(jEta, jXi);
        double fY = GrdArr.GrdArrRho.LAT(jEta, jXi);
        double eDist = CompDist(eX, eY, fX, fY);
        if (eDist > MaxDist)
          MaxDist = eDist;
      }
    }
  }
  return MaxDist;
}

QuadArray GetQuadArray(MyMatrix<double> const &LON, MyMatrix<double> const &LAT,
                       double const &deltaLL) {
  double MinLon = LON.minCoeff() - deltaLL;
  double MaxLon = LON.maxCoeff() + deltaLL;
  double MinLat = LAT.minCoeff() - deltaLL;
  double MaxLat = LAT.maxCoeff() + deltaLL;
  return {MinLon, MaxLon, MinLat, MaxLat};
}

std::ostream &operator<<(std::ostream &os, QuadArray const &eQ) {
  os << "(Lon(min/max)=" << eQ.MinLon << " / " << eQ.MaxLon
     << " Lat(min/max)=" << eQ.MinLat << " / " << eQ.MaxLat << ")";
  return os;
}

std::vector<std::string> GetAllPossibleModels() {
  std::vector<std::string> vec{"COSMO",        "WAM",
                               "ROMS",         "ROMS_IVICA",
                               "WWM",          "WWM_DAILY",
                               "WW3",          "GRIB_DWD",
                               "GRIB_ALADIN",  "GRIB_ECMWF",
                               "GRIB_GFS",     "GRIB_IFS",
                               "GRIB_COSMO",   "GRIB_WAM_FORT30",
                               "SCHISM_SFLUX", "SCHISM_NETCDF_OUT",
                               "RECTANGULAR",  "WRF",
                               "UNRUNOFF",     "IVICA_UVP",
                               "NEMO",         "HYCOM",
                               "AREG",         "GEOS"};
  return vec;
}

std::string GetKernelModelName(std::string const &eModelName) {
  std::vector<std::string> ListStr = STRING_Split(eModelName, ":");
  return ListStr[0];
}

void CHECK_Model_Allowedness(std::string const &PreModelName) {
  std::string eModelName = GetKernelModelName(PreModelName);
  std::vector<std::string> vec = GetAllPossibleModels();
  bool isPresent = (std::find(vec.begin(), vec.end(), eModelName) != vec.end());
  if (!isPresent) {
    std::cerr << "We did not find the MODEL NAME\n";
    std::cerr << "eModelName = " << eModelName
              << "     PreModelName = " << PreModelName << "\n";
    std::cerr << "List of allowed models\n";
    for (size_t iModel = 0; iModel < vec.size(); iModel++)
      std::cerr << "iModel=" << iModel << " eModel=" << vec[iModel] << "\n";
    throw TerminalException{1};
  }
}

void InitializeIdxJdxWet(CoordGridArrayFD &eCoordGrdArr) {
  int eta = eCoordGrdArr.LON.rows();
  int xi = eCoordGrdArr.LON.cols();
  eCoordGrdArr.Idx.clear();
  eCoordGrdArr.Jdx.clear();
  int nbWet = 0;
  for (int i = 0; i < eta; i++)
    for (int j = 0; j < xi; j++)
      if (eCoordGrdArr.MSK(i, j) == 1) {
        nbWet++;
        eCoordGrdArr.Idx.push_back(i);
        eCoordGrdArr.Jdx.push_back(j);
      }
  eCoordGrdArr.nbWet = nbWet;
}

bool TestEqualityGridArray(GridArray const &GrdArr1, GridArray const &GrdArr2) {
  if (GrdArr1.IsFE != GrdArr2.IsFE)
    return false;
  if (GrdArr1.IsFE) {
    if (GrdArr1.INE.rows() != GrdArr2.INE.rows())
      return false;
    int nbTrig = GrdArr1.INE.rows();
    for (int iTrig = 0; iTrig < nbTrig; iTrig++)
      for (int i = 0; i < 3; i++)
        if (GrdArr1.INE(iTrig, i) != GrdArr2.INE(iTrig, i))
          return false;
  }
  int eta1 = GrdArr1.GrdArrRho.LON.rows();
  int eta2 = GrdArr2.GrdArrRho.LON.rows();
  if (eta1 != eta2)
    return false;
  int eta = eta1;
  int xi1 = GrdArr1.GrdArrRho.LON.cols();
  int xi2 = GrdArr2.GrdArrRho.LON.cols();
  if (xi1 != xi2)
    return false;
  if (GrdArr1.GrdArrRho.DEP.has_value() != GrdArr2.GrdArrRho.DEP.has_value())
    return false;
  //  std::cerr << "GrdArr1.ModelName = " << GrdArr1.ModelName << "\n";
  //  std::cerr << "GrdArr2.ModelName = " << GrdArr2.ModelName << "\n";
  int xi = xi1;
  double err = 0;
  if (GrdArr1.GrdArrRho.DEP) {
    const MyMatrix<double> &DEP1 = GetDEP(GrdArr1.GrdArrRho);
    const MyMatrix<double> &DEP2 = GetDEP(GrdArr2.GrdArrRho);
    for (int i = 0; i < eta; i++)
      for (int j = 0; j < xi; j++) {
        double dep1 = DEP1(i, j);
        double dep2 = DEP2(i, j);
        err += fabs(dep1 - dep2);
      }
  }
  for (int i = 0; i < eta; i++)
    for (int j = 0; j < xi; j++) {
      //      std::cerr << "i=" << i << " j=" << j << "\n";
      double lon1 = GrdArr1.GrdArrRho.LON(i, j);
      //      std::cerr << "step 1\n";
      double lon2 = GrdArr2.GrdArrRho.LON(i, j);
      //      std::cerr << "step 2\n";
      err += fabs(lon1 - lon2);
      double lat1 = GrdArr1.GrdArrRho.LAT(i, j);
      //      std::cerr << "step 3\n";
      double lat2 = GrdArr2.GrdArrRho.LAT(i, j);
      //      std::cerr << "step 4\n";
      err += fabs(lat1 - lat2);
    }
  if (err > double(1))
    return false;
  return true;
}

GridArray NC_ReadGeosGridFile(std::string const &eFile) {
  if (!IsExistingFile(eFile)) {
    std::cerr << "Missing GEOS grid file eFile=" << eFile << "\n";
    throw TerminalException{1};
  }
  std::cerr << "eFile=" << eFile << "\n";
  GridArray GrdArr;
  GrdArr.ARVD.IsAssigned = false;
  GrdArr.ARVD.Zcoordinate = false;
  GrdArr.ModelName = "GEOS";
  GrdArr.IsFE = 0;
  GrdArr.IsSpherical = true;
  //
  // Now reading the RHO part.
  //
  MyVector<double> ListLON = NC_Read1Dvariable(eFile, "lon");
  MyVector<double> ListLAT = NC_Read1Dvariable(eFile, "lat");
  int nbLon = ListLON.size();
  int nbLat = ListLAT.size();
  MyMatrix<double> LON(nbLat, nbLon), LAT(nbLat, nbLon);
  for (int iLat = 0; iLat < nbLat; iLat++)
    for (int iLon = 0; iLon < nbLon; iLon++) {
      LON(iLat, iLon) = ListLON(iLon);
      LAT(iLat, iLon) = ListLAT(iLat);
    }
  int eta_rho = LON.rows();
  int xi_rho = LON.cols();
  MyMatrix<uint8_t> MSK;
  MSK.setConstant(eta_rho, xi_rho, 1);
  GrdArr.GrdArrRho.MSK = MSK;
  GrdArr.GrdArrRho.LON = LON;
  GrdArr.GrdArrRho.LAT = LAT;
  GrdArr.GrdArrRho.ANG = CreateAngleMatrix(LON, LAT);
  return GrdArr;
}

GridArray NC_ReadAregGridFile(std::string const &eFile) {
  if (!IsExistingFile(eFile)) {
    std::cerr << "Missing AREG grid file eFile=" << eFile << "\n";
    throw TerminalException{1};
  }
  std::cerr << "eFile=" << eFile << "\n";
  GridArray GrdArr;
  GrdArr.ARVD.IsAssigned = false;
  GrdArr.ARVD.Zcoordinate = false;
  GrdArr.ModelName = "AREG";
  GrdArr.IsFE = 0;
  GrdArr.IsSpherical = true;
  //
  // Now reading the RHO part.
  //
  GrdArr.GrdArrRho.LON = NC_Read2Dvariable(eFile, "lon");
  GrdArr.GrdArrRho.LAT = NC_Read2Dvariable(eFile, "lat");
  GrdArr.GrdArrRho.DEP = NC_Read2Dvariable(eFile, "depth");
  const MyMatrix<double> &LON = GrdArr.GrdArrRho.LON;
  const MyMatrix<double> &LAT = GrdArr.GrdArrRho.LAT;
  const MyMatrix<double> &DEP = GetDEP(GrdArr.GrdArrRho);
  int ypos = GrdArr.GrdArrRho.LON.rows();
  int xpos = GrdArr.GrdArrRho.LON.cols();
  GrdArr.GrdArrRho.ANG = ZeroMatrix<double>(ypos, xpos);
  const MyMatrix<double> &ANG = GrdArr.GrdArrRho.ANG;
  MyMatrix<double> eMSK_rho_double = NC_Read2Dvariable(eFile, "fsm");
  GrdArr.GrdArrRho.MSK =
      UniversalMatrixConversion<uint8_t, double>(eMSK_rho_double);
  const MyMatrix<uint8_t> &MSK = GrdArr.GrdArrRho.MSK;
  // U
  int eta_u = ypos;
  int xi_u = xpos - 1;
  MyMatrix<uint8_t> MSKu(eta_u, xi_u);
  MyMatrix<double> DEPu(eta_u, xi_u);
  MyMatrix<double> ANGu(eta_u, xi_u);
  MyMatrix<double> LONu(eta_u, xi_u);
  MyMatrix<double> LATu(eta_u, xi_u);
  for (int i = 0; i < eta_u; i++)
    for (int j = 0; j < xi_u; j++) {
      LONu(i, j) = (LON(i, j) + LON(i, j + 1)) / double(2);
      LATu(i, j) = (LAT(i, j) + LAT(i, j + 1)) / double(2);
      DEPu(i, j) = (DEP(i, j) + DEP(i, j + 1)) / double(2);
      ANGu(i, j) = (ANG(i, j) + ANG(i, j + 1)) / double(2);
      MSKu(i, j) = MSK(i, j) * MSK(i, j + 1);
    }
  GrdArr.GrdArrU.MSK = MSKu;
  GrdArr.GrdArrU.DEP = DEPu;
  GrdArr.GrdArrU.ANG = ANGu;
  GrdArr.GrdArrU.LON = LONu;
  GrdArr.GrdArrU.LAT = LATu;
  // V
  int eta_v = ypos - 1;
  int xi_v = xpos;
  MyMatrix<uint8_t> MSKv(eta_v, xi_v);
  MyMatrix<double> DEPv(eta_v, xi_v);
  MyMatrix<double> ANGv(eta_v, xi_v);
  MyMatrix<double> LONv(eta_v, xi_v);
  MyMatrix<double> LATv(eta_v, xi_v);
  for (int i = 0; i < eta_v; i++)
    for (int j = 0; j < xi_v; j++) {
      LONv(i, j) = (LON(i, j) + LON(i + 1, j)) / double(2);
      LATv(i, j) = (LAT(i, j) + LAT(i + 1, j)) / double(2);
      DEPv(i, j) = (DEP(i, j) + DEP(i + 1, j)) / double(2);
      ANGv(i, j) = (ANG(i, j) + ANG(i + 1, j)) / double(2);
      MSKv(i, j) = MSK(i, j) * MSK(i + 1, j);
    }
  GrdArr.GrdArrV.MSK = MSKv;
  GrdArr.GrdArrV.DEP = DEPv;
  GrdArr.GrdArrV.ANG = ANGv;
  GrdArr.GrdArrV.LON = LONv;
  GrdArr.GrdArrV.LAT = LATv;
  // PSI
  int eta_psi = ypos - 1;
  int xi_psi = xpos - 1;
  MyMatrix<uint8_t> MSKp(eta_psi, xi_psi);
  MyMatrix<double> DEPp(eta_psi, xi_psi);
  MyMatrix<double> ANGp(eta_psi, xi_psi);
  MyMatrix<double> LONp(eta_psi, xi_psi);
  MyMatrix<double> LATp(eta_psi, xi_psi);
  for (int i = 0; i < eta_psi; i++)
    for (int j = 0; j < xi_psi; j++) {
      LONp(i, j) =
          (LON(i, j + 1) + LON(i + 1, j + 1) + LON(i, j) + LON(i + 1, j)) /
          double(4);
      LATp(i, j) =
          (LAT(i, j + 1) + LAT(i + 1, j + 1) + LAT(i, j) + LAT(i + 1, j)) /
          double(4);
      DEPp(i, j) =
          (DEP(i, j + 1) + DEP(i + 1, j + 1) + DEP(i, j) + DEP(i + 1, j)) /
          double(4);
      ANGp(i, j) =
          (ANG(i, j + 1) + ANG(i + 1, j + 1) + ANG(i, j) + ANG(i + 1, j)) /
          double(4);
      MSKp(i, j) =
          MSK(i, j + 1) * MSK(i + 1, j + 1) * MSK(i, j) * MSK(i + 1, j);
    }
  GrdArr.GrdArrPsi.MSK = MSKp;
  GrdArr.GrdArrPsi.DEP = DEPp;
  GrdArr.GrdArrPsi.ANG = ANGp;
  return GrdArr;
}

GridArray NC_ReadRomsGridFile(std::string const &eFile) {
  if (!IsExistingFile(eFile)) {
    std::cerr << "Missing roms grid file eFile=" << eFile << "\n";
    throw TerminalException{1};
  }
  std::cerr << "eFile = " << eFile << "\n";
  std::string xName, yName;
  if (NC_IsVar(eFile, "lon_rho")) {
    xName = "lon";
    yName = "lat";
  } else {
    xName = "x";
    yName = "y";
  }
  GridArray GrdArr;
  GrdArr.ARVD.IsAssigned = false;
  GrdArr.ARVD.Zcoordinate = false;
  GrdArr.ModelName = "ROMS";
  GrdArr.IsFE = 0;
  GrdArr.IsSpherical = true;
  // Rho part of the arrays
  GrdArr.GrdArrRho.LON = NC_Read2Dvariable(eFile, xName + "_rho");
  GrdArr.GrdArrRho.LAT = NC_Read2Dvariable(eFile, yName + "_rho");
  GrdArr.GrdArrRho.DEP = NC_Read2Dvariable(eFile, "h");
  GrdArr.GrdArrRho.pm = NC_Read2Dvariable(eFile, "pm");
  GrdArr.GrdArrRho.pn = NC_Read2Dvariable(eFile, "pn");
  GrdArr.GrdArrRho.ANG = NC_Read2Dvariable(eFile, "angle");
  const MyMatrix<double> &ANG = GrdArr.GrdArrRho.ANG;
  const MyMatrix<double> &DEP = GetDEP(GrdArr.GrdArrRho);
  int eta_rho = GrdArr.GrdArrRho.LON.rows();
  int xi_rho = GrdArr.GrdArrRho.LON.cols();
  MyMatrix<double> eMSK_rho_double = NC_Read2Dvariable(eFile, "mask_rho");
  GrdArr.GrdArrRho.MSK =
      UniversalMatrixConversion<uint8_t, double>(eMSK_rho_double);
  const MyMatrix<uint8_t> &MSK = GrdArr.GrdArrRho.MSK;
  InitializeIdxJdxWet(GrdArr.GrdArrRho);
  // U
  GrdArr.GrdArrU.LON = NC_Read2Dvariable(eFile, xName + "_u");
  GrdArr.GrdArrU.LAT = NC_Read2Dvariable(eFile, yName + "_u");
  int eta_u = GrdArr.GrdArrU.LON.rows();
  int xi_u = GrdArr.GrdArrU.LON.cols();
  MyMatrix<uint8_t> MSKu(eta_u, xi_u);
  MyMatrix<double> DEPu(eta_u, xi_u);
  MyMatrix<double> ANGu(eta_u, xi_u);
  for (int i = 0; i < eta_u; i++)
    for (int j = 0; j < xi_u; j++) {
      DEPu(i, j) = (DEP(i, j) + DEP(i, j + 1)) / double(2);
      ANGu(i, j) = (ANG(i, j) + ANG(i, j + 1)) / double(2);
      MSKu(i, j) = MSK(i, j) * MSK(i, j + 1);
    }
  GrdArr.GrdArrU.MSK = MSKu;
  GrdArr.GrdArrU.DEP = DEPu;
  GrdArr.GrdArrU.ANG = ANGu;
  InitializeIdxJdxWet(GrdArr.GrdArrU);
  // V
  GrdArr.GrdArrV.LON = NC_Read2Dvariable(eFile, xName + "_v");
  GrdArr.GrdArrV.LAT = NC_Read2Dvariable(eFile, yName + "_v");
  int eta_v = GrdArr.GrdArrV.LON.rows();
  int xi_v = GrdArr.GrdArrV.LON.cols();
  MyMatrix<uint8_t> MSKv(eta_v, xi_v);
  MyMatrix<double> DEPv(eta_v, xi_v);
  MyMatrix<double> ANGv(eta_v, xi_v);
  for (int i = 0; i < eta_v; i++)
    for (int j = 0; j < xi_v; j++) {
      DEPv(i, j) = (DEP(i, j) + DEP(i + 1, j)) / double(2);
      ANGv(i, j) = (ANG(i, j) + ANG(i + 1, j)) / double(2);
      MSKv(i, j) = MSK(i, j) * MSK(i + 1, j);
    }
  GrdArr.GrdArrV.MSK = MSKv;
  GrdArr.GrdArrV.DEP = DEPv;
  GrdArr.GrdArrV.ANG = ANGv;
  InitializeIdxJdxWet(GrdArr.GrdArrV);
  // PSI
  GrdArr.GrdArrPsi.LON = NC_Read2Dvariable(eFile, xName + "_psi");
  GrdArr.GrdArrPsi.LAT = NC_Read2Dvariable(eFile, yName + "_psi");
  int eta_psi = GrdArr.GrdArrPsi.LON.rows();
  int xi_psi = GrdArr.GrdArrPsi.LON.cols();
  MyMatrix<uint8_t> MSKp(eta_psi, xi_psi);
  MyMatrix<double> DEPp(eta_psi, xi_psi);
  MyMatrix<double> ANGp(eta_psi, xi_psi);
  for (int i = 0; i < eta_psi; i++)
    for (int j = 0; j < xi_psi; j++) {
      DEPp(i, j) =
          (DEP(i, j + 1) + DEP(i + 1, j + 1) + DEP(i, j) + DEP(i + 1, j)) /
          double(4);
      ANGp(i, j) =
          (ANG(i, j + 1) + ANG(i + 1, j + 1) + ANG(i, j) + ANG(i + 1, j)) /
          double(4);
      MSKp(i, j) =
          MSK(i, j + 1) * MSK(i + 1, j + 1) * MSK(i, j) * MSK(i + 1, j);
    }
  GrdArr.GrdArrPsi.MSK = MSKp;
  GrdArr.GrdArrPsi.DEP = DEPp;
  GrdArr.GrdArrPsi.ANG = ANGp;
  InitializeIdxJdxWet(GrdArr.GrdArrPsi);
  std::cerr << "The ROMS grid has been read\n";
  bool PrintKeyInformation = true;
  if (PrintKeyInformation) {
    auto ThePrint = [&](int const &i1, int const &i2) -> void {
      std::cerr << "(i,j) = (" << i1 << "," << i2
                << ") lon=" << GrdArr.GrdArrRho.LON(i1, i2)
                << " lat=" << GrdArr.GrdArrRho.LAT(i1, i2) << "\n";
    };
    ThePrint(0, 0);
    ThePrint(0, xi_rho - 1);
    ThePrint(eta_rho - 1, xi_rho - 1);
    ThePrint(eta_rho - 1, 0);
  }
  return GrdArr;
}

CoordGridArrayFD GRID_ExtendedPsiThi(CoordGridArrayFD const &RecRho,
                                     CoordGridArrayFD const &RecU,
                                     CoordGridArrayFD const &RecV,
                                     CoordGridArrayFD const &RecPsi) {
  int eta_rho = RecRho.LON.rows();
  int xi_rho = RecRho.LON.cols();
  int eta_psi = RecPsi.LON.rows();
  int xi_psi = RecPsi.LON.cols();
  int eta_u = RecU.LON.rows();
  //  int xi_u = RecU.LON.cols();
  //  int eta_v = RecV.LON.rows();
  int xi_v = RecV.LON.cols();
  int eta_psi2 = eta_rho + 1;
  int xi_psi2 = xi_rho + 1;
  // We have eta_u = eta_rho    ;  eta_v = eta_rho-1  ;  eta_psi = eta_rho-1
  // We have xi_u  = xi_rho -1  ;  xi_v  = xi_rho     ;  xi_psi  = xi_rho-1
  auto CompField2 = [&](MyMatrix<double> const &FieldPsi,
                        MyMatrix<double> const &FieldU,
                        MyMatrix<double> const &FieldV) -> MyMatrix<double> {
    MyMatrix<double> FieldPsi2(eta_psi2, xi_psi2);
    for (int i = 0; i < eta_psi; i++)
      for (int j = 0; j < xi_psi; j++)
        FieldPsi2(i + 1, j + 1) = FieldPsi(i, j);
    for (int i = 0; i < eta_psi; i++) {
      FieldPsi2(i + 1, 0) = 2 * FieldV(i, 0) - FieldPsi(i, 0);
      FieldPsi2(i + 1, xi_psi2 - 1) =
          2 * FieldV(i, xi_v - 1) - FieldPsi(i, xi_psi - 1);
    }
    for (int j = 0; j < xi_psi; j++) {
      FieldPsi2(0, j + 1) = 2 * FieldU(0, j) - FieldPsi(0, j);
      FieldPsi2(eta_psi2 - 1, j + 1) =
          2 * FieldU(eta_u - 1, j) - FieldPsi(eta_psi - 1, j);
    }
    FieldPsi2(0, 0) = 2 * FieldPsi2(1, 0) - FieldPsi2(2, 0);
    FieldPsi2(0, xi_psi2 - 1) =
        2 * FieldPsi2(1, xi_psi2 - 1) - FieldPsi2(2, xi_psi2 - 1);
    //
    FieldPsi2(eta_psi2 - 1, 0) =
        2 * FieldPsi2(eta_psi2 - 2, 0) - FieldPsi2(eta_psi2 - 3, 0);
    FieldPsi2(eta_psi2 - 1, xi_psi2 - 1) =
        2 * FieldPsi2(eta_psi2 - 2, xi_psi2 - 1) -
        FieldPsi2(eta_psi2 - 3, xi_psi2 - 1);
    return FieldPsi2;
  };
  MyMatrix<double> LON_psi2 = CompField2(RecPsi.LON, RecU.LON, RecV.LON);
  MyMatrix<double> LAT_psi2 = CompField2(RecPsi.LAT, RecU.LAT, RecV.LAT);
  MyMatrix<double> DEP_psi2 =
      CompField2(GetDEP(RecPsi), GetDEP(RecU), GetDEP(RecV));
  MyMatrix<double> ANG_psi2 = CompField2(RecPsi.ANG, RecU.ANG, RecV.ANG);
  //
  MyMatrix<uint8_t> MSK_psi2 = ZeroMatrix<uint8_t>(eta_psi2, xi_psi2);
  for (int i = 0; i < eta_psi; i++)
    for (int j = 0; j < xi_psi; j++)
      MSK_psi2(i + 1, j + 1) = RecPsi.MSK(i, j);
  for (int i = 0; i < eta_psi; i++) {
    MSK_psi2(i + 1, 0) = MSK_psi2(i + 1, 1);
    MSK_psi2(i + 1, xi_psi2 - 1) = MSK_psi2(i + 1, xi_psi2 - 2);
  }
  for (int j = 0; j < xi_psi; j++) {
    MSK_psi2(0, j + 1) = MSK_psi2(1, j + 1);
    MSK_psi2(eta_psi2 - 1, j + 1) = MSK_psi2(eta_psi2 - 2, j + 1);
  }
  if (MSK_psi2(1, 0) && MSK_psi2(0, 1))
    MSK_psi2(0, 0) = 1;
  if (MSK_psi2(eta_psi2 - 2, 0) && MSK_psi2(eta_psi2 - 1, 1))
    MSK_psi2(eta_psi2 - 1, 0) = 1;
  if (MSK_psi2(1, xi_psi2 - 1) && MSK_psi2(0, xi_psi2 - 2))
    MSK_psi2(0, xi_psi2 - 1) = 1;
  if (MSK_psi2(eta_psi2 - 2, xi_psi2 - 1) &&
      MSK_psi2(eta_psi2 - 1, xi_psi2 - 2))
    MSK_psi2(eta_psi2 - 1, xi_psi2 - 1) = 1;
  //
  CoordGridArrayFD RecPsi2;
  RecPsi2.MSK = MSK_psi2;
  RecPsi2.LON = LON_psi2;
  RecPsi2.LAT = LAT_psi2;
  RecPsi2.DEP = DEP_psi2;
  RecPsi2.ANG = ANG_psi2;
  InitializeIdxJdxWet(RecPsi2);
  return RecPsi2;
}

MyMatrix<double> GetMatrixRadiusROMS(GridArray const &GrdArr) {
  CoordGridArrayFD RecPsi2 = GRID_ExtendedPsiThi(
      GrdArr.GrdArrRho, GrdArr.GrdArrU, GrdArr.GrdArrV, GrdArr.GrdArrPsi);
  MyMatrix<int> MatDir(4, 2);
  MatDir(0, 0) = 0;
  MatDir(0, 1) = 0;
  MatDir(1, 0) = 1;
  MatDir(1, 1) = 0;
  MatDir(2, 0) = 1;
  MatDir(2, 1) = 1;
  MatDir(3, 0) = 0;
  MatDir(3, 1) = 1;
  int eta_rho = GrdArr.GrdArrRho.LON.rows();
  int xi_rho = GrdArr.GrdArrRho.LON.cols();
  MyMatrix<double> MatRadius(eta_rho, xi_rho);
  for (int i = 0; i < eta_rho; i++)
    for (int j = 0; j < xi_rho; j++) {
      double lon1 = GrdArr.GrdArrRho.LON(i, j);
      double lat1 = GrdArr.GrdArrRho.LON(i, j);
      double MinDist = 45555555;
      for (int u = 0; u < 4; u++) {
        int i2 = i + MatDir(u, 0);
        int j2 = j + MatDir(u, 1);
        double lon2 = RecPsi2.LON(i2, j2);
        double lat2 = RecPsi2.LON(i2, j2);
        double dx = lon1 - lon2;
        double dy = lat1 - lat2;
        double dist = sqrt(dx * dx + dy * dy);
        if (dist < MinDist)
          MinDist = dist;
      }
      MatRadius(i, j) = MinDist;
    }
  return MatRadius;
}

GridArray NC_ReadWrfGridFile(std::string const &eFile) {
  GridArray GrdArr;
  GrdArr.ARVD.IsAssigned = false;
  GrdArr.ARVD.Zcoordinate = false;
  GrdArr.ModelName = "WRF";
  GrdArr.IsFE = 0;
  GrdArr.IsSpherical = true;
  // Rho part of the arrays
  MyMatrix<double> LON = NC_Read2Dvariable(eFile, "XLONG");
  MyMatrix<double> LAT = NC_Read2Dvariable(eFile, "XLAT");
  GrdArr.GrdArrRho.LON = LON;
  GrdArr.GrdArrRho.LAT = LAT;
  GrdArr.GrdArrRho.ANG = CreateAngleMatrix(LON, LAT);
  int eta_rho = LON.rows();
  int xi_rho = LON.cols();
  MyMatrix<uint8_t> MSK;
  MSK.setConstant(eta_rho, xi_rho, 1);
  GrdArr.GrdArrRho.MSK = MSK;
  InitializeIdxJdxWet(GrdArr.GrdArrRho);
  return GrdArr;
}

GridArray NC_ReadHycomGridFile(std::string const &eFile) {
  if (!IsExistingFile(eFile)) {
    std::cerr << "Error in NC_ReadHycomGridFile\n";
    std::cerr << "Trying to open non-existing file\n";
    std::cerr << "eFile = " << eFile << "\n";
    throw TerminalException{1};
  }
  std::cerr << "eFile=" << eFile << "\n";
  std::cerr << "NC_ReadHycomGridFile, step 1\n";
  GridArray GrdArr;
  GrdArr.ModelName = "HYCOM";
  GrdArr.IsFE = 0;
  GrdArr.IsSpherical = true;
  // Rho part of the arrays
  MyVector<double> lon1d = NC_Read1Dvariable(eFile, "lon");
  MyVector<double> lat1d = NC_Read1Dvariable(eFile, "lat");
  MyVector<double> dep1d_pre = NC_Read1Dvariable(eFile, "depth");

  Eigen::Tensor<int, 3> eTens = NC_Read3Dvariable_Mask_file(eFile, "surf_el");
  MyMatrix<int> MSK2 = StrictProjectionMask(eTens, 0);
  Eigen::Tensor<int, 3> eTens3 =
      NC_Read3Dvariable_Mask_file(eFile, "water_temp_bottom");
  MyMatrix<int> MSK3 = StrictProjectionMask(eTens3, 0);
  Eigen::Tensor<int, 3> eTens4 =
      NC_Read3Dvariable_Mask_file(eFile, "water_u_bottom");
  MyMatrix<int> MSK4 = StrictProjectionMask(eTens4, 0);
  Eigen::Tensor<int, 3> eTens5 =
      NC_Read3Dvariable_Mask_file(eFile, "water_v_bottom");
  MyMatrix<int> MSK5 = StrictProjectionMask(eTens5, 0);
  Eigen::Tensor<int, 3> eTens6 =
      NC_Read3Dvariable_Mask_file(eFile, "salinity_bottom");
  MyMatrix<int> MSK6 = StrictProjectionMask(eTens6, 0);
  std::cerr << "sum(MSK2)=" << MSK2.sum() << "\n";
  std::cerr << "sum(MSK3)=" << MSK3.sum() << "\n";
  std::cerr << "sum(MSK4)=" << MSK4.sum() << "\n";
  std::cerr << "sum(MSK5)=" << MSK5.sum() << "\n";
  std::cerr << "sum(MSK6)=" << MSK6.sum() << "\n";
  std::cerr << "NC_ReadHycomGridFile, step 2\n";
  int nbLon = lon1d.size();
  int nbLat = lat1d.size();
  int nbDep = dep1d_pre.size();
  std::cerr << "nbLon=" << nbLon << " nbLat=" << nbLat << "\n";
  MyVector<double> dep1d(nbDep);
  // We want index 0 to be deepest and index nbDep-1 to be near surface
  for (int iDep = 0; iDep < nbDep; iDep++)
    dep1d(nbDep - 1 - iDep) = dep1d_pre(iDep);
  for (int iDep = 0; iDep < nbDep; iDep++)
    std::cerr << "iDep=" << iDep << " dep1d=" << dep1d(iDep) << "\n";
  /*
  for (int iDep=0; iDep<nbDep; iDep++)
  std::cerr << "iDep=" << iDep << " dep1d=" << dep1d(iDep) << "\n";*/
  MyMatrix<double> LON(nbLat, nbLon);
  MyMatrix<double> LAT(nbLat, nbLon);
  for (int i = 0; i < nbLat; i++)
    for (int j = 0; j < nbLon; j++) {
      LON(i, j) = lon1d(j);
      LAT(i, j) = lat1d(i);
    }
  std::cerr << "NC_ReadHycomGridFile, step 3\n";
  CheckNetcdfDataArray("NC_ReadHycomGridFile", eFile, "salinity");
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
  netCDF::NcVar data = dataFile.getVar("salinity");
  std::cerr << "NC_ReadHycomGridFile, step 4\n";
  MyVector<uint8_t> StatusFill = NC_ReadVariable_StatusFill_data<uint8_t>(data);
  MyVector<double> VarFill = NC_ReadVariable_data(data);
  size_t TotalSize = StatusFill.size();
  std::cerr << "|StatusFill|=" << TotalSize
            << " min/max=" << int(StatusFill.minCoeff()) << " / "
            << int(StatusFill.maxCoeff()) << " sum=" << int(StatusFill.sum())
            << "\n";
  std::vector<size_t> ListDim = NC_ReadVariable_listdim(data);
  /*
  for (int iTotal=0; iTotal<TotalSize; iTotal++) {
    std::cerr << "iTotal=" << iTotal << " StatusFill=" << StatusFill(iTotal) <<
  " VarFill=" << VarFill(iTotal) << "\n";
    }*/
  int nbTime = ListDim[0];
  int s_vert = ListDim[1];
  int eta = ListDim[2];
  int xi = ListDim[3];
  std::cerr << "nbTime=" << nbTime << " nbDep=" << nbDep << " eta=" << eta
            << " xi=" << xi << "\n";
  if (eta != nbLat || xi != nbLon || s_vert != nbDep) {
    std::cerr << "eta=" << eta << " nbLat=" << nbLat << "\n";
    std::cerr << "xi=" << xi << " nbLon=" << nbLon << "\n";
    std::cerr << "s_vert=" << s_vert << " nbDep=" << nbDep << "\n";
    std::cerr << "Inconsistency between array sizes. Logical error\n";
    throw TerminalException{1};
  }
  //
  // Computing MSK and DEP
  //
  std::cerr << "NC_ReadHycomGridFile, step 5\n";
  MyMatrix<uint8_t> MSK(nbLat, nbLon);
  MyMatrix<double> DEP(nbLat, nbLon);
  Eigen::Tensor<int, 4> StatusTens(nbTime, nbDep, nbLat, nbLon);
  Eigen::Tensor<double, 4> VarTens(nbTime, nbDep, nbLat, nbLon);
  MyMatrix<int> StatusSum = ZeroMatrix<int>(nbLat, nbLon);
  int idx = 0;
  for (int iTime = 0; iTime < nbTime; iTime++)
    for (int iDep = 0; iDep < nbDep; iDep++)
      for (int i = 0; i < nbLat; i++)
        for (int j = 0; j < nbLon; j++) {
          StatusTens(iTime, nbDep - 1 - iDep, i, j) = StatusFill(idx);
          VarTens(iTime, nbDep - 1 - iDep, i, j) = VarFill(idx);
          StatusSum(i, j) += StatusFill(idx);
          idx++;
        }
  //
  std::cerr << "sum(StatusSum)=" << StatusSum.sum()
            << " sum(StatusFill)=" << StatusFill.sum() << "\n";
  std::cerr << "max(StatusSum)=" << StatusSum.maxCoeff() << "\n";
  /*
  std::cerr << "StatusSum:\n";
  for (int i=0; i<nbLat; i++)
    for (int j=0; j<nbLon; j++)
      std::cerr << "i=" << i << " j=" << j << " StatusSum=" << StatusSum(i,j) <<
  " lon=" << LON(i,j) << " lat=" << LAT(i,j) << "\n";
  */
  //
  bool CoherencyCheck = true;
  if (CoherencyCheck) {
    for (int iTime = 0; iTime < nbTime; iTime++)
      for (int i = 0; i < nbLat; i++)
        for (int j = 0; j < nbLon; j++) {
          for (int iDep = 1; iDep < nbDep; iDep++) {
            if (StatusTens(iTime, iDep - 1, i, j) == 0 &&
                StatusTens(iTime, iDep, i, j) == 1) {
              std::cerr << "Found inconsistency at i=" << i << " j=" << j
                        << " iTime = " << iTime << " iDep=" << iDep << "\n";
            }
          }
        }
    std::cerr << "After coherency checks\n";
  }
  std::cerr << "NC_ReadHycomGridFile, step 6\n";
  int ValLand = nbTime * nbDep;
  std::vector<int> VectVal;
  for (int i = 0; i < nbLat; i++)
    for (int j = 0; j < nbLon; j++) {
      VectVal.push_back(StatusSum(i, j));
      int eVal = 1;
      if (StatusSum(i, j) == ValLand)
        eVal = 0;
      MSK(i, j) = eVal;
    }
  CollectedResult<int> eColl = Collected(VectVal);
  for (size_t u = 0; u < eColl.LVal.size(); u++) {
    int eVal = eColl.LVal[u];
    int res = eVal % 245;
    int ldep = (eVal - res) / 245;
    std::cerr << "u=" << u << " eVal=" << eVal << " res=" << res
              << " ldep=" << ldep << " eMult=" << eColl.LMult[u] << "\n";
  }

  int eProd = nbLat * nbLon;
  int nb1_0 = 0;
  int nb0_1 = 0;
  for (int i = 0; i < nbLat; i++)
    for (int j = 0; j < nbLon; j++) {
      if (MSK(i, j) == 0 && MSK2(i, j) == 1)
        nb0_1 = nb0_1 + 1;
      if (MSK(i, j) == 1 && MSK2(i, j) == 0)
        nb1_0 = nb1_0 + 1;
    }
  std::cerr << "nb0_1=" << nb0_1 << " nb1_0=" << nb1_0 << "\n";
  std::cerr << "sum(MSK)=" << MSK.sum() << " sum(MSK2)=" << MSK2.sum()
            << " eProd=" << eProd << "\n";

  size_t minCoeff = MSK.minCoeff();
  size_t maxCoeff = MSK.maxCoeff();
  size_t sumCoeff = MSK.sum();
  std::cerr << "MSK min=" << minCoeff << " / " << maxCoeff
            << " sum=" << sumCoeff << " eProd=" << eProd << "\n";
  std::cerr << "ValLand=" << ValLand << "\n";
  int iTimeRef = 0;
  /*
  int nb48=0;
  for (int i=0; i<nbLat; i++)
    for (int j=0; j<nbLon; j++) {
      if (StatusSum(i,j) == 48) {
        nb48++;
        double lat=lat1d(i);
        double lon=lon1d(j);
        std::cerr << "i=" << i << " j=" << j << " lon/lat=" << lon << " / " <<
  lat << "\n"; for (int iDep=0; iDep<nbDep; iDep++) std::cerr << "iDep=" << iDep
  << " dep=" << dep1d(iDep) << " status=" << StatusTens(iTimeRef,iDep,i,j) <<
  "\n";
      }
    }
    std::cerr << "nb48=" << nb48 << "\n";*/
  std::cerr << "StatusSum  min/max=" << StatusSum.minCoeff() << " / "
            << StatusSum.maxCoeff() << "\n";
  int minCoeffB = MSK.minCoeff();
  int maxCoeffB = MSK.maxCoeff();
  int sumCoeffB = MSK.sum();
  std::cerr << "HYCOM MSK min / max / sum=" << minCoeffB << " / " << maxCoeffB
            << " / " << sumCoeffB << "\n";
  std::cerr << "nbLat=" << nbLat << " nbLon=" << nbLon << "\n";
  for (int i = 0; i < nbLat; i++)
    for (int j = 0; j < nbLon; j++) {
      double eDep = 0;
      if (MSK(i, j) == 1) {
        if (StatusSum(i, j) > 0) {
          for (int iDep = 1; iDep < nbDep; iDep++) {
            if (StatusTens(iTimeRef, iDep - 1, i, j) == 1 &&
                StatusTens(iTimeRef, iDep, i, j) == 0) {
              double dep1 = dep1d(iDep - 1); // rock
              double dep2 = dep1d(iDep);     // sea
              eDep = (dep1 + dep2) / double(2);
            }
          }
        } else {
          eDep = dep1d(0);
        }
        if (eDep == 0) {
          int sumStatus = 0;
          for (int iDep = 0; iDep < nbDep; iDep++)
            sumStatus += StatusTens(iTimeRef, iDep, i, j);
          if (sumStatus != nbDep) {
            for (int iDep = 0; iDep < nbDep; iDep++)
              std::cerr << "iDep=" << iDep
                        << " status=" << StatusTens(iTimeRef, iDep, i, j)
                        << " dep=" << dep1d(iDep) << "\n";
            double lon = lon1d(j);
            double lat = lat1d(i);
            std::cerr << "i=" << i << " j=" << j << " lon=" << lon
                      << " lat=" << lat << "\n";
            std::cerr << "sumStatus=" << sumStatus << " nbDep=" << nbDep
                      << " StatusSum(i,j)=" << StatusSum(i, j)
                      << " MSK=" << MSK(i, j) << "\n";
            std::cerr << "sumStatus is not equal to nbDep. Error in "
                         "NC_ReadHycomGridFile\n";
            throw TerminalException{1};
          }
          eDep = dep1d(nbDep - 1);
        }
      }
      DEP(i, j) = eDep;
    }
  std::cerr << "DEP min/max=" << DEP.minCoeff() << " / " << DEP.maxCoeff()
            << "\n";
  double sumTEM = 0;
  for (int i = 0; i < nbLat; i++)
    for (int j = 0; j < nbLon; j++) {
      if (MSK(i, j) == 1) {
        sumTEM += VarTens(iTimeRef, nbDep - 1, i, j);
      }
    }
  std::cerr << "sumTEM = " << sumTEM << "\n";
  GrdArr.ARVD.N = nbDep;
  GrdArr.ARVD.IsAssigned = true;
  GrdArr.ARVD.ListZ_r = -dep1d;
  GrdArr.ARVD.Zcoordinate = true;
  GrdArr.ARVD.ModelName = "HYCOM";
  //
  // More standard assignation
  //
  GrdArr.GrdArrRho.LON = LON;
  GrdArr.GrdArrRho.LAT = LAT;
  GrdArr.GrdArrRho.DEP = DEP;
  GrdArr.GrdArrRho.ANG = CreateAngleMatrix(LON, LAT);
  GrdArr.GrdArrRho.MSK = MSK;
  InitializeIdxJdxWet(GrdArr.GrdArrRho);
  std::cerr << "Leaving NC_ReadHycomGridFile\n";
  return GrdArr;
}

NEMO_vars ReadNEMO_vars(std::string const &eFile) {
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
  if (dataFile.isNull()) {
    std::cerr << "Unexpected error\n";
    throw TerminalException{1};
  }
  std::vector<std::string> ListVarForbid = {"depth",     "latitude", "lat",
                                            "longitude", "lon",      "time"};
  std::vector<std::string> ListVar = NC_ListVar(eFile);
  std::vector<std::string> List2D_vars;
  std::vector<std::string> List3D_vars;
  for (auto &eVar : ListVar) {
    if (PositionVect(ListVarForbid, eVar) == -1) {
      netCDF::NcVar data = dataFile.getVar(eVar);
      if (data.isNull()) {
        std::cerr << "Unexpected error\n";
        throw TerminalException{1};
      }
      int nbDim = data.getDimCount();
      if (nbDim == 4) { // Worked with 4-dim var only so far, that is time,
                        // vertical, geographic
        List3D_vars.push_back(eVar);
      }
      if (nbDim ==
          3) { // Worked with 3-dim var only so far, that is time, vgeographic
        List2D_vars.push_back(eVar);
      }
    }
  }
  return {List2D_vars, List3D_vars};
}

GridArray NC_ReadNemoGridFile(std::string const &eFile) {
  if (!IsExistingFile(eFile)) {
    std::cerr << "Error in NC_ReadNemoGridFile\n";
    std::cerr << "Trying to open non-existing file\n";
    std::cerr << "eFile = " << eFile << "\n";
    throw TerminalException{1};
  }
  std::cerr << "NC_ReadNemoGridFile with eFile=" << eFile << "\n";
  GridArray GrdArr;
  GrdArr.ModelName = "NEMO";
  GrdArr.IsFE = 0;
  GrdArr.IsSpherical = true;
  //
  // The longitude/latitude part of the grid
  //
  MyVector<double> lon1d, lat1d;
  if (NC_IsVar(eFile, "lon")) {
    lon1d = NC_Read1Dvariable(eFile, "lon");
  } else {
    lon1d = NC_Read1Dvariable(eFile, "longitude");
  }
  if (NC_IsVar(eFile, "lat")) {
    lat1d = NC_Read1Dvariable(eFile, "lat");
  } else {
    lat1d = NC_Read1Dvariable(eFile, "latitude");
  }
  std::cerr << "NC_ReadNemoGridFile, step 1\n";
  int nbLon = lon1d.size();
  int nbLat = lat1d.size();
  /*
  for (int iDep=0; iDep<nbDep; iDep++)
  std::cerr << "iDep=" << iDep << " dep1d=" << dep1d(iDep) << "\n";*/
  MyMatrix<double> LON(nbLat, nbLon);
  MyMatrix<double> LAT(nbLat, nbLon);
  for (int i = 0; i < nbLat; i++)
    for (int j = 0; j < nbLon; j++) {
      LON(i, j) = lon1d(j);
      LAT(i, j) = lat1d(i);
    }
  std::cerr << "We have LON/LAT\n";
  //
  // The depth of the grid
  //
  MyVector<double> dep1d_pre = NC_Read1Dvariable(eFile, "depth");
  int nbDep = dep1d_pre.size();
  MyVector<double> dep1d(nbDep);
  // We want index 0 to be deepest and index nbDep-1 to be near surface
  for (int iDep = 0; iDep < nbDep; iDep++)
    dep1d(nbDep - 1 - iDep) = dep1d_pre(iDep);
  for (int iDep = 0; iDep < nbDep; iDep++)
    std::cerr << "iDep=" << iDep << " z=" << dep1d(iDep) << "\n";
  std::cerr << "We have DEP\n";
  //
  // Now computing the mask, bathymetry and so on.
  //
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
  if (dataFile.isNull()) {
    std::cerr << "Unexpected error\n";
    throw TerminalException{1};
  }
  auto get_var_test = [&]() -> std::string {
    std::vector<std::string> ListVarForbid = {"depth",     "latitude", "lat",
                                              "longitude", "lon",      "time"};
    std::vector<std::string> ListVar = NC_ListVar(eFile);
    for (auto &eVar : ListVar) {
      if (PositionVect(ListVarForbid, eVar) == -1) {
        netCDF::NcVar data = dataFile.getVar(eVar);
        if (data.isNull()) {
          std::cerr << "Unexpected error\n";
          throw TerminalException{1};
        }
        int nbDim = data.getDimCount();
        if (nbDim == 4) { // Worked with 4-dim var only so far, that is time,
                          // vertical, geographic
          return eVar;
        }
      }
    }
    std::cerr << "Failed to find a matching variable\n";
    std::cerr << "Maybe we need to extend the functionality for supporting 3 "
                 "dim vars\n";
    throw TerminalException{1};
  };

  std::cerr << "NC_ReadNemoGridFile, step 1\n";
  std::string eVar = get_var_test();
  std::cerr << "Found variable test eVar=" << eVar << "\n";
  netCDF::NcVar data = dataFile.getVar(eVar);
  if (data.isNull()) {
    std::cerr << "Error while reading thetao\n";
    throw TerminalException{1};
  }
  std::cerr << "We have data\n";
  std::vector<size_t> ListDim = NC_ReadVariable_listdim(data);
  size_t nbTime = ListDim[0];
  size_t s_vert = ListDim[1];
  size_t eta = ListDim[2];
  size_t xi = ListDim[3];
  size_t nbTimeWork, nbTimeCrit = 10;
  if (nbTime > nbTimeCrit) {
    nbTimeWork = nbTimeCrit;
  } else {
    nbTimeWork = nbTime;
  }
  std::cerr << "nbTime=" << nbTime << " nbDep=" << nbDep << " eta=" << eta
            << " xi=" << xi << "\n";

  std::vector<size_t> start{0, 0, 0, 0};
  std::vector<size_t> count{nbTimeWork, s_vert, eta, xi};
  // StatusFill = 1 correspond to missing value.
  MyVector<uint8_t> StatusFill =
      NC_ReadVariable_StatusFill_data_start_count<uint8_t>(data, start, count);
  MyVector<int> StatusFill_i =
      UniversalVectorConversion<int, uint8_t>(StatusFill);
  std::cerr << "We have StatusFill\n";
  MyVector<double> VarFill =
      NC_ReadVariable_data_start_count(data, start, count);
  std::cerr << "|StatusFill|=" << StatusFill.size()
            << " min/max=" << int(StatusFill.minCoeff()) << " / "
            << int(StatusFill.maxCoeff()) << " sum=" << StatusFill_i.sum()
            << "\n";
  if (eta != size_t(nbLat) || xi != size_t(nbLon) || s_vert != size_t(nbDep)) {
    std::cerr << "eta=" << eta << " nbLat=" << nbLat << "\n";
    std::cerr << "xi=" << xi << " nbLon=" << nbLon << "\n";
    std::cerr << "s_vert=" << s_vert << " nbDep=" << nbDep << "\n";
    std::cerr << "Inconsistency between array sizes. Logical error\n";
    throw TerminalException{1};
  }
  size_t tot_siz = VarFill.size();
  for (int iter = 0; iter < 100; iter++) {
    size_t pos = rand() % tot_siz;
    std::cerr << "iter=" << iter << " StatusFill_i=" << StatusFill_i(pos)
              << " VarFill=" << VarFill(pos) << "\n";
  }

  //
  // Computing MSK and DEP
  //
  MyMatrix<uint8_t> MSK(nbLat, nbLon);
  MyMatrix<double> DEP(nbLat, nbLon);
  Eigen::Tensor<uint8_t, 4> StatusTens(nbTimeWork, nbDep, nbLat, nbLon);
  Eigen::Tensor<double, 4> VarTens(nbTime, nbDep, nbLat, nbLon);
  MyMatrix<int> StatusSum = ZeroMatrix<int>(nbLat, nbLon);
  int idx = 0;
  Eigen::Tensor<uint64_t, 3> TensMSKvert_64(nbDep, nbLat, nbLon);
  for (int iDep = 0; iDep < nbDep; iDep++)
    for (int i = 0; i < nbLat; i++)
      for (int j = 0; j < nbLon; j++)
        TensMSKvert_64(iDep, i, j) = 0;
  for (size_t iTime = 0; iTime < nbTimeWork; iTime++)
    for (int iDep = 0; iDep < nbDep; iDep++)
      for (int i = 0; i < nbLat; i++)
        for (int j = 0; j < nbLon; j++) {
          StatusTens(iTime, nbDep - 1 - iDep, i, j) = StatusFill(idx);
          VarTens(iTime, nbDep - 1 - iDep, i, j) = VarFill(idx);
          StatusSum(i, j) += StatusFill(idx);
          TensMSKvert_64(nbDep - 1 - iDep, i, j) += 1 - StatusFill(idx);
          idx++;
        }
  Eigen::Tensor<uint8_t, 3> TensMSKvert(nbDep, nbLat, nbLon);
  for (int iDep = 0; iDep < nbDep; iDep++) {
    int nWet = 0;
    size_t n_error = 0;
    std::map<int, size_t> m_tens;
    for (int i = 0; i < nbLat; i++) {
      for (int j = 0; j < nbLon; j++) {
        int val = TensMSKvert_64(iDep, i, j);
        m_tens[val]++;
        if (val != 0 && val != int(nbTimeWork)) {
          //          std::cerr << "Inconsistency in the TensMSKvert\n";
          n_error++;
        }
        uint8_t val_8 = 0;
        if (val > 0)
          val_8 = 1;
        TensMSKvert(iDep, i, j) = val_8;
        nWet += val_8;
      }
    }
    std::cerr << "iDep=" << iDep << " nWet=" << nWet << " n_error=" << n_error
              << " m_tens=";
    for (auto &kv : m_tens)
      std::cerr << " (" << kv.first << " / " << kv.second << ")";
    std::cerr << "\n";
  }
  bool CoherencyCheck = true;
  if (CoherencyCheck) {
    for (size_t iTime = 0; iTime < nbTimeWork; iTime++)
      for (int i = 0; i < nbLat; i++)
        for (int j = 0; j < nbLon; j++) {
          for (int iDep = 1; iDep < nbDep; iDep++) {
            // StatusTens = 0 corresponds to sea, and StatusTens = 1 to land
            if (StatusTens(iTime, iDep - 1, i, j) == 0 &&
                StatusTens(iTime, iDep, i, j) == 1) {
              std::cerr << "Found inconsistency at i=" << i << " j=" << j
                        << " iTime = " << iTime << " iDep=" << iDep << "\n";
            }
          }
        }
    std::cerr << "After coherency checks\n";
  }
  int ValLand = nbTimeWork * nbDep;
  for (int i = 0; i < nbLat; i++)
    for (int j = 0; j < nbLon; j++) {
      int eVal = 1;
      if (StatusSum(i, j) == ValLand)
        eVal = 0;
      MSK(i, j) = eVal;
    }
  int eProd = nbLat * nbLon;
  MyMatrix<int> MSK_i = UniversalMatrixConversion<int, uint8_t>(MSK);
  std::cerr << "MSK min=" << int(MSK.minCoeff())
            << " max=" << int(MSK.maxCoeff()) << " sum=" << MSK_i.sum()
            << " eProd=" << eProd << "\n";
  std::cerr << "ValLand=" << ValLand << "\n";
  int iTimeRef = 0;
  //  std::cerr << "StatusFill min/max=" << StatusFill.minCoeff() << " / " <<
  //  StatusFill.maxCoeff() << "\n";
  std::cerr << "StatusSum  min/max=" << StatusSum.minCoeff() << " / "
            << StatusSum.maxCoeff() << "\n";
  std::cerr << "NEMO MSK min / max / sum=" << int(MSK.minCoeff()) << " / "
            << int(MSK.maxCoeff()) << " / " << int(MSK.sum()) << "\n";
  std::cerr << "nbLat=" << nbLat << " nbLon=" << nbLon << "\n";
  for (int i = 0; i < nbLat; i++)
    for (int j = 0; j < nbLon; j++) {
      double eDep = 0;
      if (MSK(i, j) == 1) {
        if (StatusSum(i, j) > 0) {
          for (int iDep = 1; iDep < nbDep; iDep++) {
            if (StatusTens(iTimeRef, iDep - 1, i, j) == 1 &&
                StatusTens(iTimeRef, iDep, i, j) == 0) {
              double dep1 = dep1d(iDep - 1); // rock
              double dep2 = dep1d(iDep);     // sea
              eDep = (dep1 + dep2) / double(2);
            }
          }
        } else {
          eDep = dep1d(0);
        }
        if (eDep == 0) {
          int sumStatus = 0;
          for (int iDep = 0; iDep < nbDep; iDep++)
            sumStatus += StatusTens(iTimeRef, iDep, i, j);
          if (sumStatus != nbDep) {
            for (int iDep = 0; iDep < nbDep; iDep++)
              std::cerr << "iDep=" << iDep
                        << " status=" << StatusTens(iTimeRef, iDep, i, j)
                        << " dep=" << dep1d(iDep) << "\n";
            double lon = lon1d(j);
            double lat = lat1d(i);
            std::cerr << "i=" << i << " j=" << j << " lon=" << lon
                      << " lat=" << lat << "\n";
            std::cerr << "sumStatus=" << sumStatus << " nbDep=" << nbDep
                      << " StatusSum(i,j)=" << StatusSum(i, j)
                      << " MSK=" << MSK(i, j) << "\n";
            std::cerr << "sumStatus is not equal to nbDep. Error in "
                         "NC_ReadNemoGridFile\n";
            throw TerminalException{1};
          }
          eDep = dep1d(nbDep - 1);
        }
      }
      DEP(i, j) = eDep;
    }
  std::cerr << "DEP min/max=" << DEP.minCoeff() << " / " << DEP.maxCoeff()
            << "\n";
  double sumTEM = 0;
  for (int i = 0; i < nbLat; i++)
    for (int j = 0; j < nbLon; j++) {
      if (MSK(i, j) == 1) {
        sumTEM += VarTens(iTimeRef, nbDep - 1, i, j);
      }
    }
  std::cerr << "sumTEM = " << sumTEM << "\n";
  GrdArr.ARVD.N = nbDep;
  GrdArr.ARVD.IsAssigned = true;
  GrdArr.ARVD.ListZ_r = -dep1d;
  GrdArr.ARVD.TensMSKvert = TensMSKvert;
  GrdArr.ARVD.Zcoordinate = true;
  GrdArr.ARVD.ModelName = "NEMO";
  //
  // More standard assignation
  //
  GrdArr.GrdArrRho.LON = LON;
  GrdArr.GrdArrRho.LAT = LAT;
  GrdArr.GrdArrRho.DEP = DEP;
  GrdArr.GrdArrRho.ANG = CreateAngleMatrix(LON, LAT);
  GrdArr.GrdArrRho.MSK = MSK;
  InitializeIdxJdxWet(GrdArr.GrdArrRho);
  return GrdArr;
}

GridArray NC_ReadCosmoWamStructGridFile(std::string const &eFile,
                                        std::string const &postfix) {
  GridArray GrdArr;
  GrdArr.ARVD.IsAssigned = false;
  GrdArr.ARVD.Zcoordinate = false;
  if (postfix == "atm")
    GrdArr.ModelName = "COSMO";
  else
    GrdArr.ModelName = "WAM";
  GrdArr.IsFE = 0;
  GrdArr.IsSpherical = true;
  // Rho part of the arrays
  std::string LONstr = "LON_" + postfix;
  //  std::cerr << "Before reading " << LONstr << "\n";
  GrdArr.GrdArrRho.LON = NC_Read2Dvariable(eFile, LONstr);
  std::string LATstr = "LAT_" + postfix;
  //  std::cerr << "Before reading " << LATstr << "\n";
  GrdArr.GrdArrRho.LAT = NC_Read2Dvariable(eFile, LATstr);
  int eta_rho = GrdArr.GrdArrRho.LON.rows();
  int xi_rho = GrdArr.GrdArrRho.LON.cols();
  //  std::cerr << "eta_rho=" << eta_rho << " xi_rho=" << xi_rho << "\n";
  // The bathymetry if available
  std::string DEPstr = "DEP_" + postfix;
  if (NC_IsVar(eFile, DEPstr)) {
    GrdArr.GrdArrRho.DEP = NC_Read2Dvariable(eFile, DEPstr);
  }
  // The mask if available
  std::string MSKstr = "MSK_" + postfix;
  MyMatrix<double> MSK_double(eta_rho, xi_rho);
  if (NC_IsVar(eFile, DEPstr)) {
    MSK_double = NC_Read2Dvariable(eFile, MSKstr);
  } else {
    for (int i = 0; i < eta_rho; i++)
      for (int j = 0; j < xi_rho; j++)
        MSK_double(i, j) = double(1);
  }
  MyMatrix<uint8_t> MSK_int(eta_rho, xi_rho);
  for (int i = 0; i < eta_rho; i++)
    for (int j = 0; j < xi_rho; j++) {
      MSK_int(i, j) = uint8_t(MSK_double(i, j));
    }
  GrdArr.GrdArrRho.MSK = MSK_int;
  // The angle if available
  //  std::cerr << "Before reading angle\n";
  if (postfix == "atm")
    GrdArr.GrdArrRho.ANG = NC_Read2Dvariable(eFile, "ANG_atm");
  else
    GrdArr.GrdArrRho.ANG =
        CreateAngleMatrix(GrdArr.GrdArrRho.LON, GrdArr.GrdArrRho.LAT);
  // U / V / PSI we do not need a priori
  return GrdArr;
}

GridArray NC_ReadSCHISM_sflux_grid(std::string const &eFile) {
  GridArray GrdArr;
  GrdArr.ARVD.IsAssigned = false;
  GrdArr.ARVD.Zcoordinate = false;
  GrdArr.ModelName = "SCHISM_SFLUX";
  GrdArr.IsFE = 0;
  GrdArr.IsSpherical = true;
  // Rho part of the arrays
  std::string LONstr = "lon";
  GrdArr.GrdArrRho.LON = NC_Read2Dvariable(eFile, LONstr);
  std::string LATstr = "lat";
  GrdArr.GrdArrRho.LAT = NC_Read2Dvariable(eFile, LATstr);
  int eta_rho = GrdArr.GrdArrRho.LON.rows();
  int xi_rho = GrdArr.GrdArrRho.LON.cols();
  //
  MyMatrix<uint8_t> MSK_int(eta_rho, xi_rho);
  for (int i = 0; i < eta_rho; i++)
    for (int j = 0; j < xi_rho; j++)
      MSK_int(i, j) = 1;
  GrdArr.GrdArrRho.MSK = MSK_int;
  GrdArr.GrdArrRho.ANG =
      CreateAngleMatrix(GrdArr.GrdArrRho.LON, GrdArr.GrdArrRho.LAT);
  // U / V / PSI we do not need a priori
  return GrdArr;
}

MyMatrix<int> NC_ReadElements(std::string const &eFile,
                              std::string const &eStr) {
  MyMatrix<int> INE = NC_Read2Dvariable_int(eFile, eStr);
  int nbRow = INE.rows();
  int nbCol = INE.cols();
  for (int iRow = 0; iRow < nbRow; iRow++)
    for (int iCol = 0; iCol < nbCol; iCol++) {
      int IP = INE(iRow, iCol);
      INE(iRow, iCol) = IP - 1;
    }
  return INE;
}

GridArray NC_ReadWamGridFile(std::string const &eFile) {
  MyVector<int> eVectLLUNSTR = NC_Read1Dvariable_int(eFile, "LLUNSTR");
  int LLUNSTR = eVectLLUNSTR(0);
  //  std::cerr << "LLUNSTR=" << LLUNSTR << "\n";
  if (LLUNSTR == 0)
    return NC_ReadCosmoWamStructGridFile(eFile, "wav");
  GridArray GrdArr;
  GrdArr.ARVD.IsAssigned = false;
  GrdArr.ARVD.Zcoordinate = false;
  GrdArr.ModelName = "WAM";
  GrdArr.IsFE = 1;
  GrdArr.IsSpherical = true;
  GrdArr.L_IndexSelect = false;
  //
  GrdArr.INE = NC_ReadElements(eFile, "ele");
  //  std::cerr << "NC_ReadWam, step 1\n";
  MyVector<double> LON = NC_Read1Dvariable(eFile, "LON_wav");
  //  std::cerr << "NC_ReadWam, step 2\n";
  MyVector<double> LAT = NC_Read1Dvariable(eFile, "LAT_wav");
  //  std::cerr << "NC_ReadWam, step 3\n";
  MyVector<double> DEP = NC_Read1Dvariable(eFile, "DEP_wav");
  //  std::cerr << "NC_ReadWam, step 4\n";
  int nbPoint = LON.size();
  MyMatrix<double> LONarr(nbPoint, 1);
  MyMatrix<double> LATarr(nbPoint, 1);
  MyMatrix<double> DEParr(nbPoint, 1);
  MyMatrix<double> ANGarr(nbPoint, 1);
  MyMatrix<uint8_t> MSKarr(nbPoint, 1);
  //  std::cerr << "NC_ReadWam, step 5\n";
  for (int iPoint = 0; iPoint < nbPoint; iPoint++) {
    LONarr(iPoint, 0) = LON(iPoint);
    LATarr(iPoint, 0) = LAT(iPoint);
    DEParr(iPoint, 0) = DEP(iPoint);
    ANGarr(iPoint, 0) = 0;
    MSKarr(iPoint, 0) = 1;
  }
  GrdArr.GrdArrRho.LON = LONarr;
  GrdArr.GrdArrRho.LAT = LATarr;
  GrdArr.GrdArrRho.DEP = DEParr;
  GrdArr.GrdArrRho.ANG = ANGarr;
  GrdArr.GrdArrRho.MSK = MSKarr;
  return GrdArr;
}

GridArray WWM_ReadGridFile_netcdf(std::string const &GridFile) {
  GridArray GrdArr;
  GrdArr.ARVD.IsAssigned = false;
  GrdArr.ARVD.Zcoordinate = false;
  GrdArr.IsFE = 1;
  GrdArr.L_IndexSelect = false;
  //
  if (!IsExistingFile(GridFile)) {
    std::cerr << "Error in WWM_ReadGridFile_netcdf\n";
    std::cerr << "GridFile = " << GridFile << "\n";
    std::cerr << "is missing\n";
    throw TerminalException{1};
  }
  GrdArr.INE = NC_ReadElements(GridFile, "ele");
  bool IsVar_X = NC_IsVar(GridFile, "x");
  bool IsVar_Y = NC_IsVar(GridFile, "y");
  bool IsVar_lon = NC_IsVar(GridFile, "lon");
  bool IsVar_lat = NC_IsVar(GridFile, "lat");
  if (IsVar_X != IsVar_Y) {
    std::cerr << "IsVar_X/Y=" << IsVar_X << " / " << IsVar_Y << "\n";
    std::cerr << "They should be identical\n";
    throw TerminalException{1};
  }
  if (IsVar_lon != IsVar_lat) {
    std::cerr << "IsVar_lon/lat=" << IsVar_lon << " / " << IsVar_lat << "\n";
    std::cerr << "They should be identical\n";
    throw TerminalException{1};
  }
  if ((IsVar_X && IsVar_lon) || (!IsVar_X && !IsVar_lon)) {
    std::cerr << "IsVar_X/lon=" << IsVar_X << " / " << IsVar_lon << "\n";
    std::cerr << "One and exactly one option should be selected\n";
    throw TerminalException{1};
  }
  bool LSPHE = IsVar_lon;
  std::string Xname, Yname;
  if (LSPHE) {
    Xname = "lon";
    Yname = "lat";
    GrdArr.IsSpherical = true;
  } else {
    Xname = "x";
    Yname = "y";
    GrdArr.IsSpherical = false;
  }
  MyVector<double> LON = NC_Read1Dvariable(GridFile, Xname);
  MyVector<double> LAT = NC_Read1Dvariable(GridFile, Yname);
  MyVector<double> DEP = NC_Read1Dvariable(GridFile, "depth");
  MyVector<double> IOBP = NC_Read1Dvariable(GridFile, "IOBP");
  int nbPoint = LON.size();
  MyMatrix<double> LONarr(nbPoint, 1);
  MyMatrix<double> LATarr(nbPoint, 1);
  MyMatrix<double> DEParr(nbPoint, 1);
  MyMatrix<double> ANGarr(nbPoint, 1);
  MyMatrix<uint8_t> MSKarr(nbPoint, 1);
  MyVector<int> IOBParr(nbPoint);
  for (int iPoint = 0; iPoint < nbPoint; iPoint++) {
    LONarr(iPoint, 0) = LON(iPoint);
    LATarr(iPoint, 0) = LAT(iPoint);
    DEParr(iPoint, 0) = DEP(iPoint);
    ANGarr(iPoint, 0) = 0;
    MSKarr(iPoint, 0) = 1;
    IOBParr(iPoint) = int(IOBP(iPoint));
  }
  GrdArr.GrdArrRho.LON = LONarr;
  GrdArr.GrdArrRho.LAT = LATarr;
  GrdArr.GrdArrRho.DEP = DEParr;
  GrdArr.GrdArrRho.ANG = ANGarr;
  GrdArr.GrdArrRho.MSK = MSKarr;
  GrdArr.IOBP = IOBParr;
  return GrdArr;
}

// IOBP should be the WWM type of IOBP
// Possible values are:
// --- 0: normal point
// --- 1: Island
// --- 2: Dirichlet
// --- 3: Neumann
// --- 4: New mixed kind
void WriteGridFile_msh(std::string const &GridFile, GridArray const &GrdArr) {
  if (GrdArr.IsFE == 0) {
    std::cerr << "We need the grid to be finite element\n";
    throw TerminalException{1};
  }

  std::ofstream os(GridFile);
  os << std::fixed;
  os << std::setprecision(9);
  int np_total = GrdArr.GrdArrRho.LON.rows();
  std::vector<int> IPbound;
  std::vector<int> IPisland;
  if (GrdArr.IOBP.size() != np_total) {
    std::cerr << "We have |GrdArr.IOBP|=" << GrdArr.IOBP.size() << "\n";
    std::cerr << "when it should be np_total=" << np_total << "\n";
    std::cerr << "Most likely you forgot to provide for boundary file\n";
    throw TerminalException{1};
  }
  for (int i = 0; i < np_total; i++) {
    int eVal = GrdArr.IOBP(i);
    if (eVal == 2) {
      IPbound.push_back(i);
    }
    if (eVal == 1 || eVal == 3 || eVal == 4) {
      IPisland.push_back(i);
    }
  }
  int nbDirichlet = IPbound.size();
  int nbIsland = IPisland.size();
  std::vector<int> ACTIVE(nbDirichlet, 1);
  os << "$MeshFormat\n";
  os << "2 0 8\n";
  os << "$EndMeshFormat\n";
  os << "$Nodes\n";
  os << np_total << "\n";
  const MyMatrix<double> &DEP = GetDEP(GrdArr.GrdArrRho);
  for (int i = 0; i < np_total; i++) {
    int IP = i + 1;
    double eXP = GrdArr.GrdArrRho.LON(i, 0);
    double eYP = GrdArr.GrdArrRho.LAT(i, 0);
    double eDEP = DEP(i, 0);
    os << IP << " " << eXP << " " << eYP << " " << eDEP << "\n";
  }
  //
  // Writing the elements
  //
  int ne_total = GrdArr.INE.rows();
  int nbEle = ne_total + nbDirichlet + nbIsland;
  os << "$EndNodes\n";
  os << "$Elements\n";
  os << nbEle << "\n";
  //
  // Write boundary
  //
  int ie2 = 0;
  for (int i = 0; i < nbDirichlet; i++) {
    ie2++;
    int eAct = ACTIVE[i];
    os << ie2 << " 15 2 " << eAct << " 0 " << IPbound[i] << "\n";
  }
  //
  // Write island
  //
  for (int i = 0; i < nbIsland; i++) {
    ie2++;
    int ip = i + 1;
    int eIPisland = IPisland[i] + 1;
    os << ie2 << " 15 2 0 " << ip << " " << eIPisland << "\n";
  }
  //
  // Write gmsh elements
  //
  for (int ie = 0; ie < ne_total; ie++) {
    ie2++;
    int ieRel = ie + 1;
    int ip1 = GrdArr.INE(ie, 0) + 1;
    int ip2 = GrdArr.INE(ie, 1) + 1;
    int ip3 = GrdArr.INE(ie, 2) + 1;
    os << ie2 << " 2 3 0 " << ieRel << " 0 " << ip1 << " " << ip2 << " " << ip3
       << "\n";
  }
  os << "$EndElements\n";
}

void WriteWWMboundaryGR3(std::string const &BndFile, GridArray const &GrdArr) {
  int nbVert = GrdArr.IOBP.size();
  int nbNode = GrdArr.GrdArrRho.LON.size();
  int nbTrig = GrdArr.INE.rows();
  if (nbVert != nbNode) {
    std::cerr << "We have nbVert not equal to nbNode\n";
    std::cerr << "nbNode = " << nbNode << "\n";
    std::cerr << "nbVert = " << nbVert << "\n";
    throw TerminalException{1};
  }
  std::ofstream os(BndFile);
  os << std::fixed;
  os << std::setprecision(9);
  os << "boundary file\n";
  os << nbTrig << " " << nbNode << "\n";
  for (int iVert = 0; iVert < nbVert; iVert++) {
    double eLon = GrdArr.GrdArrRho.LON(iVert, 0);
    double eLat = GrdArr.GrdArrRho.LAT(iVert, 0);
    int iPos = iVert + 1;
    int eIOBP = GrdArr.IOBP(iVert);
    os << iPos << " " << eLon << " " << eLat << " " << eIOBP << "\n";
  }
}

MyVector<int> WWM_ReadBoundFile_gr3(std::string const &BoundFile) {
  if (!IsExistingFile(BoundFile)) {
    std::cerr << "Error in WWM_ReadBoundFile_gr3\n";
    std::cerr << "Missing BoundFile=" << BoundFile << "\n";
    throw TerminalException{1};
  }
  std::ifstream IN(BoundFile);
  std::string line;
  std::getline(IN, line);
  int mne, mnp;
  IN >> mne;
  IN >> mnp;
  MyVector<int> eVect(mnp);
  for (int i = 0; i < mnp; i++) {
    int KTMP;
    double XPDTMP, YPDTMP, ZPDTMP;
    IN >> KTMP >> XPDTMP >> YPDTMP >> ZPDTMP;
    if (KTMP != i + 1) {
      std::cerr << "Inconsistency at this level\n";
      throw TerminalException{1};
    }
    int eIOBP = int(ZPDTMP);
    eVect(i) = eIOBP;
  }
  return eVect;
}

/* This strategy fails if the domain has a width of less than 90m */
bool GuessIsSpherical(GridArray const &GrdArr) {
  bool IsSpherical = true;
  double LONmax = GrdArr.GrdArrRho.LON.maxCoeff();
  double LONmin = GrdArr.GrdArrRho.LON.minCoeff();
  double deltaLON = LONmax - LONmin;
  if (deltaLON > 360)
    IsSpherical = false;
  double LATmax = GrdArr.GrdArrRho.LAT.maxCoeff();
  double LATmin = GrdArr.GrdArrRho.LAT.minCoeff();
  if (LATmax > 90 || LATmin < -90)
    IsSpherical = false;
  return IsSpherical;
}

GridArray WWM_ReadGridFile_msh(std::string const &GridFile) {
  GridArray GrdArr;
  std::ifstream is(GridFile);
  std::string line;
  std::getline(is, line);
  if (line != "$MeshFormat") {
    std::cerr << "line=" << line << "\n";
    std::cerr << "First line should be $MeshFormat\n";
    throw TerminalException{1};
  }
  std::getline(is, line);
  std::cerr << "GMSH version number / fileType / dataSize = " << line << "\n";
  //
  std::getline(is, line);
  if (line != "$EndMeshFormat") {
    std::cerr << "line=" << line << "\n";
    std::cerr << "Line should be $EndMeshFormat\n";
    throw TerminalException{1};
  }
  //
  auto WaitForString = [&](std::string const &strSearch) -> void {
    while (true) {
      try {
        std::string strRead;
        std::getline(is, strRead);
        if (strRead == strSearch)
          return;
      } catch (...) {
        std::cerr << "Error in data reading\n";
        throw TerminalException{1};
      }
    }
  };
  //
  WaitForString("$Nodes");
  //
  std::getline(is, line);
  int mnp;
  std::istringstream(line) >> mnp;
  std::cerr << "mnp=" << mnp << "\n";
  MyMatrix<double> LON(mnp, 1), LAT(mnp, 1), DEP(mnp, 1);
  for (int ip = 0; ip < mnp; ip++) {
    int idx;
    double eLon, eLat, eDep;
    std::getline(is, line);
    std::vector<std::string> LStr = STRING_Split(line, " ");
    if (LStr.size() < 4) {
      std::cerr << "|LStr|=" << LStr.size() << "\n";
      std::cerr << "line=" << line << "\n";
      std::cerr << "Unfortunately LStr is too small\n";
      throw TerminalException{1};
    }
    std::istringstream(LStr[0]) >> idx;
    std::istringstream(LStr[1]) >> eLon;
    std::istringstream(LStr[2]) >> eLat;
    std::istringstream(LStr[3]) >> eDep;
    if (idx != ip + 1) {
      std::cerr << "idx=" << idx << " eLon=" << eLon << " eLat=" << eLat
                << " eDep=" << eDep << "\n";
      std::cerr << "idx=" << idx << " ip=" << ip << "\n";
      std::cerr << "Error in the indices\n";
      throw TerminalException{1};
    }
    LON(ip, 0) = eLon;
    LAT(ip, 0) = eLat;
    DEP(ip, 0) = eDep;
  }
  GrdArr.GrdArrRho.LON = LON;
  GrdArr.GrdArrRho.LAT = LAT;
  GrdArr.GrdArrRho.DEP = DEP;
  //
  std::getline(is, line);
  if (line != "$EndNodes") {
    std::cerr << "line=" << line << "\n";
    std::cerr << "Line should be $EndNodes\n";
    throw TerminalException{1};
  }
  //
  WaitForString("$Elements");
  //
  std::getline(is, line);
  int nbElemTot;
  std::istringstream(line) >> nbElemTot;
  int nb15 = 0;
  int nb1 = 0;
  int nb2 = 0;
  int nb3 = 0;
  std::cerr << "nbElemTot=" << nbElemTot << "\n";
  std::vector<int> List15, List1, List2, List3;
  for (int ie = 0; ie < nbElemTot; ie++) {
    std::getline(is, line);
    std::vector<int> LInt = STRING_Split_Int(line, " ");
    int len = LInt.size();
    if (len < 2) {
      std::cerr << "len is not large enough. len=" << len << "\n";
      throw TerminalException{1};
    }
    int idx = LInt[0];
    int elm_type = LInt[1];
    if (idx != ie + 1) {
      std::cerr << "idx=" << idx << " ie=" << ie << "\n";
      std::cerr << "Error in the indices\n";
      throw TerminalException{1};
    }
    if (elm_type != 15 && elm_type != 1 && elm_type != 2 && elm_type != 3) {
      std::cerr << "elm_type=" << elm_type << "\n";
      std::cerr << "Only values allowed are 15, 1, 2 and 3\n";
      throw TerminalException{1};
    }
    if (elm_type == 15) {
      int a4 = LInt[len - 1];
      List15.push_back(a4 - 1);
      nb15++;
    }
    if (elm_type == 1) {
      int a1 = LInt[len - 2];
      int a2 = LInt[len - 1];
      List1.push_back(a1 - 1);
      List1.push_back(a2 - 1);
      nb1++;
    }
    if (elm_type == 2) {
      int a1 = LInt[len - 3];
      int a2 = LInt[len - 2];
      int a3 = LInt[len - 1];
      List2.push_back(a1 - 1);
      List2.push_back(a2 - 1);
      List2.push_back(a3 - 1);
      nb2++;
    }
    if (elm_type == 3) {
      int a1 = LInt[len - 4];
      int a2 = LInt[len - 3];
      int a3 = LInt[len - 2];
      int a4 = LInt[len - 1];
      List3.push_back(a1 - 1);
      List3.push_back(a2 - 1);
      List3.push_back(a3 - 1);
      List3.push_back(a4 - 1);
      nb3++;
    }
  }
  if (nb2 > 0 && nb3 > 0) {
    std::cerr << "nb2=" << nb2 << " nb3=" << nb3 << "\n";
    std::cerr << "Right now we cannot mix triangles and quadrangles\n";
    throw TerminalException{1};
  }
  std::getline(is, line);
  if (line != "$EndElements") {
    std::cerr << "line=" << line << "\n";
    std::cerr << "We reach an error here. Should be $EndElements\n";
    throw TerminalException{1};
  }
  if (nb2 > 0) {
    MyMatrix<int> INE(nb2, 3);
    int idx = 0;
    for (int ie = 0; ie < nb2; ie++) {
      INE(ie, 0) = List2[idx];
      idx++;
      INE(ie, 1) = List2[idx];
      idx++;
      INE(ie, 2) = List2[idx];
      idx++;
    }
    GrdArr.INE = INE;
  }
  if (nb3 > 0) {
    MyMatrix<int> INE(nb3, 4);
    int idx = 0;
    for (int ie = 0; ie < nb3; ie++) {
      INE(ie, 0) = List3[idx];
      idx++;
      INE(ie, 1) = List3[idx];
      idx++;
      INE(ie, 2) = List3[idx];
      idx++;
      INE(ie, 3) = List3[idx];
      idx++;
    }
    GrdArr.INE = INE;
  }
  GrdArr.ARVD.IsAssigned = false;
  GrdArr.ARVD.Zcoordinate = false;
  GrdArr.IsSpherical = true;
  return GrdArr;
}

GridArray WWM_ReadGridFile_gr3(std::string const &GridFile) {
  GridArray GrdArr;
  GrdArr.ARVD.IsAssigned = false;
  GrdArr.ARVD.Zcoordinate = false;
  GrdArr.IsFE = 1;
  GrdArr.L_IndexSelect = false;
  //
  if (!IsExistingFile(GridFile)) {
    std::cerr << "Error in WWM_ReadGridFile_gr3\n";
    std::cerr << "GridFile = " << GridFile << "\n";
    std::cerr << "is missing\n";
    throw TerminalException{1};
  }
  std::ifstream IN(GridFile);
  // read first line
  std::string line;
  std::getline(IN, line);
  std::cerr << "line=" << line << "\n";
  //
  int mne, mnp;
  IN >> mne;
  IN >> mnp;
  std::cerr << "mne=" << mne << " mnp=" << mnp << "\n";
  GrdArr.INE = MyMatrix<int>(mne, 3);
  GrdArr.GrdArrRho.LON = MyMatrix<double>(mnp, 1);
  GrdArr.GrdArrRho.LAT = MyMatrix<double>(mnp, 1);
  MyMatrix<double> DEP(mnp, 1);
  GrdArr.GrdArrRho.ANG = MyMatrix<double>(mnp, 1);
  GrdArr.GrdArrRho.MSK = MyMatrix<uint8_t>(mnp, 1);
  for (int iP = 0; iP < mnp; iP++) {
    int KTMP;
    double XPDTMP, YPDTMP, ZPDTMP;
    IN >> KTMP >> XPDTMP >> YPDTMP >> ZPDTMP;
    //    std::cerr << "iP=" << iP << " XYZ=" << XPDTMP << " " << YPDTMP << " "
    //    << ZPDTMP << "\n";
    GrdArr.GrdArrRho.LON(iP, 0) = XPDTMP;
    GrdArr.GrdArrRho.LAT(iP, 0) = YPDTMP;
    DEP(iP, 0) = ZPDTMP;
    GrdArr.GrdArrRho.ANG(iP, 0) = 0;
    GrdArr.GrdArrRho.MSK(iP, 0) = 1;
  }
  GrdArr.GrdArrRho.DEP = DEP;
  auto check_ip = [&](int const &ip) -> void {
    if (ip < 1 || ip > mnp) {
      std::cerr << "We have ip=" << ip
                << " but it should be in the interval [1..mnp] with mnp=" << mnp
                << "\n";
      throw TerminalException{1};
    }
  };
  std::vector<int> Nmatch(mnp, 0);
  int KTMP, LTMP, ip, idx;
  for (int iE = 0; iE < mne; iE++) {
    IN >> KTMP >> LTMP;
    for (int i = 0; i < 3; i++) {
      IN >> ip;
      check_ip(ip);
      idx = ip - 1;
      GrdArr.INE(iE, i) = idx;
      Nmatch[idx]++;
    }
  }
  for (int ip = 0; ip < mnp; ip++) {
    if (Nmatch[ip] == 0) {
      std::cerr << "Some vertices are contained in ZERO triangles. Clear bug\n";
      throw TerminalException{1};
    }
  }

  GrdArr.IsSpherical = GuessIsSpherical(GrdArr);
  return GrdArr;
}

MyVector<int> WWM_ReadBoundFile_DAT(std::string const &BoundFile) {
  if (!IsExistingFile(BoundFile)) {
    std::cerr << "Error in WWM_ReadBoundFile_DAT\n";
    std::cerr << "Missing BoundFile=" << BoundFile << "\n";
    throw TerminalException{1};
  }
  std::ifstream IN(BoundFile);
  std::string line;
  for (int i = 0; i < 2; i++)
    std::getline(IN, line);
  int ITMP, JTMP;
  IN >> ITMP;
  std::getline(IN, line);
  IN >> JTMP;
  int mnp = ITMP + JTMP;
  for (int i = 0; i < 7; i++)
    std::getline(IN, line);
  MyVector<int> eVect(mnp);
  for (int i = 0; i < mnp; i++) {
    int KTMP;
    double XPDTMP, YPDTMP, ZPDTMP;
    IN >> KTMP >> XPDTMP >> YPDTMP >> ZPDTMP;
    if (KTMP != i + 1) {
      std::cerr << "Inconsistency error\n";
      throw TerminalException{1};
    }
    int eIOBP = int(ZPDTMP);
    eVect(i) = eIOBP;
  }
  return eVect;
}

GridArray WWM_ReadGridFile_obj(std::string const &GridFile) {
  GridArray GrdArr;
  GrdArr.ARVD.IsAssigned = false;
  GrdArr.ARVD.Zcoordinate = false;
  GrdArr.IsFE = 1;
  GrdArr.L_IndexSelect = false;
  //
  if (!IsExistingFile(GridFile)) {
    std::cerr << "Error in WWM_ReadGridFile_obj\n";
    std::cerr << "GridFile = " << GridFile << "\n";
    std::cerr << "is missing\n";
    throw TerminalException{1};
  }
  std::string line;
  std::ifstream IN(GridFile);
  std::getline(IN, line);
  //  std::cerr << "line=" << line << "\n";
  std::vector<std::vector<double>> ListVect;
  std::vector<std::vector<int>> ListTrig;
  while (true) {
    if (IN.eof() == 1)
      break;
    std::getline(IN, line);
    if (line.size() == 0)
      break;
    //    std::cerr << "line=" << line << " |line|=" << line.size() << "\n";
    std::istringstream is(line);
    std::string eChar;
    is >> eChar;
    bool IsDone = false;
    if (eChar == "v") {
      IsDone = true;
      double lon, lat, dep;
      is >> lon >> lat >> dep;
      ListVect.push_back({lon, lat, dep});
    }
    if (eChar == "f") {
      IsDone = true;
      int i1, i2, i3;
      is >> i1 >> i2 >> i3;
      ListTrig.push_back({i1 - 1, i2 - 1, i3 - 1});
    }
    if (!IsDone) {
      std::cerr << "Error while reading data\n";
      std::cerr << "eChar=" << eChar << "\n";
      throw TerminalException{1};
    }
  }
  int mnp = ListVect.size();
  GrdArr.GrdArrRho.LON = MyMatrix<double>(mnp, 1);
  GrdArr.GrdArrRho.LAT = MyMatrix<double>(mnp, 1);
  MyMatrix<double> DEP(mnp, 1);
  GrdArr.GrdArrRho.ANG = MyMatrix<double>(mnp, 1);
  GrdArr.GrdArrRho.MSK = MyMatrix<uint8_t>(mnp, 1);
  for (int i = 0; i < mnp; i++) {
    GrdArr.GrdArrRho.LON(i, 0) = ListVect[i][0];
    GrdArr.GrdArrRho.LAT(i, 0) = ListVect[i][1];
    DEP(i, 0) = ListVect[i][2];
    GrdArr.GrdArrRho.ANG(i, 0) = 0;
    GrdArr.GrdArrRho.MSK(i, 0) = 1;
  }
  GrdArr.GrdArrRho.DEP = DEP;
  int mne = ListTrig.size();
  GrdArr.INE = MyMatrix<int>(mne, 3);
  for (int i = 0; i < mne; i++)
    for (int j = 0; j < 3; j++)
      GrdArr.INE(i, j) = ListTrig[i][j];
  GrdArr.IsSpherical = GuessIsSpherical(GrdArr);
  return GrdArr;
}

GridArray WWM_ReadGridFile_DAT(std::string const &GridFile) {
  GridArray GrdArr;
  GrdArr.ARVD.IsAssigned = false;
  GrdArr.ARVD.Zcoordinate = false;
  GrdArr.IsFE = 1;
  GrdArr.L_IndexSelect = false;
  //
  if (!IsExistingFile(GridFile)) {
    std::cerr << "Error in WWM_ReadGridFile_DAT\n";
    std::cerr << "GridFile = " << GridFile << "\n";
    std::cerr << "is missing\n";
    throw TerminalException{1};
  }
  std::string line;
  int ITMP, JTMP;
  std::ifstream IN(GridFile);
  for (int i = 0; i < 2; i++)
    std::getline(IN, line);
  std::getline(IN, line);
  std::istringstream(line) >> ITMP;
  std::cerr << "ITMP=" << ITMP << "\n";
  std::getline(IN, line);
  std::getline(IN, line);
  std::istringstream(line) >> JTMP;
  std::cerr << "JTMP=" << JTMP << "\n";
  int mnp = ITMP + JTMP;
  std::cerr << "mnp=" << mnp << "\n";
  for (int i = 0; i < 7; i++)
    std::getline(IN, line);
  GrdArr.GrdArrRho.LON = MyMatrix<double>(mnp, 1);
  GrdArr.GrdArrRho.LAT = MyMatrix<double>(mnp, 1);
  MyMatrix<double> DEP(mnp, 1);
  GrdArr.GrdArrRho.ANG = MyMatrix<double>(mnp, 1);
  GrdArr.GrdArrRho.MSK = MyMatrix<uint8_t>(mnp, 1);
  for (int iP = 0; iP < mnp; iP++) {
    int KTMP;
    double XPDTMP, YPDTMP, ZPDTMP;
    std::getline(IN, line);
    std::istringstream(line) >> KTMP >> XPDTMP >> YPDTMP >> ZPDTMP;
    if (KTMP != iP) {
      std::cerr << "KTMP=" << KTMP << " iP=" << iP << "\n";
      std::cerr << "Inconsistency in the values\n";
      throw TerminalException{1};
    }
    GrdArr.GrdArrRho.LON(iP, 0) = XPDTMP;
    GrdArr.GrdArrRho.LAT(iP, 0) = YPDTMP;
    DEP(iP, 0) = ZPDTMP;
    GrdArr.GrdArrRho.ANG(iP, 0) = 0;
    GrdArr.GrdArrRho.MSK(iP, 0) = 1;
  }
  GrdArr.GrdArrRho.DEP = DEP;
  for (int i = 0; i < 2; i++)
    std::getline(IN, line);
  std::getline(IN, line);
  int mne;
  std::istringstream(line) >> mne;
  GrdArr.INE = MyMatrix<int>(mne, 3);
  for (int i = 0; i < 3; i++)
    std::getline(IN, line);
  for (int iE = 0; iE < mne; iE++) {
    int KTMP, LTMP, ip1, ip2, ip3;
    std::getline(IN, line);
    std::istringstream(line) >> ip1 >> ip2 >> ip3 >> KTMP >> LTMP;
    if (LTMP != iE || KTMP != 0) {
      std::cerr << "LTMP=" << LTMP << " iE=" << iE << " (should be equal)\n";
      std::cerr << "KTMP=" << KTMP << " (should be 0)\n";
      std::cerr << "Inconsistent values in the .dat file\n";
      throw TerminalException{1};
    }
    GrdArr.INE(iE, 0) = ip1;
    GrdArr.INE(iE, 1) = ip2;
    GrdArr.INE(iE, 2) = ip3;
  }
  GrdArr.IsSpherical = GuessIsSpherical(GrdArr);
  return GrdArr;
}

GridArray NC_ReadWW3_GridFile(std::string const &eFile) {
  GridArray GrdArr;
  GrdArr.ARVD.IsAssigned = false;
  GrdArr.ARVD.Zcoordinate = false;
  GrdArr.ModelName = "WWM";
  GrdArr.IsFE = 1;
  GrdArr.L_IndexSelect = false;
  std::cerr << "NC_ReadWW3_GRidFile\n";
  //
  GrdArr.INE = NC_ReadElements(eFile, "tri");
  MyVector<double> LON = NC_Read1Dvariable(eFile, "longitude");
  MyVector<double> LAT = NC_Read1Dvariable(eFile, "latitude");
  int nbPoint = LON.size();
  MyMatrix<double> LONarr(nbPoint, 1);
  MyMatrix<double> LATarr(nbPoint, 1);
  MyMatrix<double> DEParr(nbPoint, 1);
  MyMatrix<double> ANGarr(nbPoint, 1);
  MyMatrix<uint8_t> MSKarr(nbPoint, 1);
  for (int iPoint = 0; iPoint < nbPoint; iPoint++) {
    LONarr(iPoint, 0) = LON(iPoint);
    LATarr(iPoint, 0) = LAT(iPoint);
    DEParr(iPoint, 0) = 0;
    ANGarr(iPoint, 0) = 0;
    MSKarr(iPoint, 0) = 1;
  }
  GrdArr.GrdArrRho.LON = LONarr;
  GrdArr.GrdArrRho.LAT = LATarr;
  GrdArr.GrdArrRho.DEP = DEParr;
  GrdArr.GrdArrRho.ANG = ANGarr;
  GrdArr.GrdArrRho.MSK = MSKarr;
  GrdArr.ModelName = "WW3";
  GrdArr.IsSpherical = GuessIsSpherical(GrdArr);
  return GrdArr;
}

int TheSignFct(double const &eVal) {
  if (eVal > 0)
    return 1;
  if (eVal < 0)
    return -1;
  return 0;
}

MyMatrix<double> get_angle_corr_rho(MyMatrix<double> const &LON_rho,
                                    MyMatrix<double> const &LAT_rho) {
  std::string spheroid = "wgs84";
  double A = -1, B = -1, E = -1;
  if (spheroid == "sph") {
    A = 6371000.0;
    B = A;
    E = sqrt(A * A - B * B) / A;
  }
  if (spheroid == "cla") {
    A = 6378206.4E0;
    B = 6356583.8E0;
    E = sqrt(A * A - B * B) / A;
  }
  if (spheroid == "iau") {
    A = 6378160.e0;
    B = 6356774.516E0;
    E = sqrt(A * A - B * B) / A;
  }
  if (spheroid == "wgs84") {
    A = 6378137.;
    E = 0.081819191;
    double alpha = A * E;
    B = sqrt(A * A - alpha * alpha);
  }
  double eps = E * E / (1 - E * E);
  int eta_rho = LON_rho.rows();
  int xi_rho = LON_rho.cols();
  int eta_u = eta_rho;
  int xi_u = xi_rho - 1;
  MyMatrix<double> LONrad_u(eta_u, xi_u);
  MyMatrix<double> LATrad_u(eta_u, xi_u);
  double pi = 3.1415926535;
  double eFact = pi / double(360);
  for (int i = 0; i < eta_u; i++)
    for (int j = 0; j < xi_u; j++) {
      double eLON = (LON_rho(i, j) + LON_rho(i, j + 1)) * eFact;
      double eLAT = (LAT_rho(i, j) + LAT_rho(i, j + 1)) * eFact;
      if (eLAT == 0)
        eLAT = eps;
      LONrad_u(i, j) = eLON;
      LATrad_u(i, j) = eLAT;
    }
  MyMatrix<double> azim(eta_u, xi_u - 1);
  for (int i = 0; i < eta_u; i++)
    for (int j = 0; j < xi_u - 1; j++) {
      double PHI1 = LATrad_u(i, j);
      double XLAM1 = LONrad_u(i, j);
      double PHI2 = LATrad_u(i, j + 1);
      double XLAM2 = LONrad_u(i, j + 1);
      if (PHI1 == PHI2)
        PHI2 = PHI2 + 1e-14;
      if (XLAM1 == XLAM2)
        XLAM2 = XLAM2 + 1e-14;
      //
      double EsPHI1 = E * sin(PHI1);
      double EsPHI2 = E * sin(PHI2);
      double xnu1 = A / sqrt(1 - EsPHI1 * EsPHI1);
      double xnu2 = A / sqrt(1 - EsPHI2 * EsPHI2);
      double TPSI2 = (1 - E * E) * tan(PHI2) +
                     E * E * xnu1 * sin(PHI1) / (xnu2 * cos(PHI2));
      double DLAM = XLAM2 - XLAM1;
      double CTA12 = (cos(PHI1) * TPSI2 - sin(PHI1) * cos(DLAM)) / sin(DLAM);
      double DLAM2 = DLAM;
      if (DLAM2 >= pi)
        DLAM2 = DLAM2 - 2 * pi;
      if (DLAM2 <= -pi)
        DLAM2 = DLAM2 + 2 * pi;
      double eAzim = atan(1 / CTA12);
      if (eAzim < -pi)
        eAzim += 2 * pi;
      if (eAzim > pi)
        eAzim += -2 * pi;
      if (TheSignFct(eAzim) != TheSignFct(DLAM2)) {
        eAzim += pi * double(TheSignFct(-eAzim));
      }
      azim(i, j) = eAzim;
    }
  MyMatrix<double> angle(eta_rho, xi_rho);
  for (int i = 0; i < eta_rho; i++)
    for (int j = 1; j < xi_u; j++) {
      double eAzim = azim(i, j - 1);
      double eAngle = (pi / double(2)) - eAzim;
      angle(i, j) = eAngle;
    }
  for (int i = 0; i < eta_rho; i++) {
    double eAngle = angle(i, 1);
    angle(i, 0) = eAngle;
    eAngle = angle(i, xi_u - 1);
    angle(i, xi_u) = eAngle;
  }
  return angle;
}

void DifferenceLonRenormalize(double &Lon) {
  if (Lon > 180)
    Lon = Lon - 360;
  if (Lon < -180)
    Lon = Lon + 360;
}

GridArray CFONE_GRID_ARRAY(std::string const &GridFile) {
  MyVector<double> LonArr = NC_Read1Dvariable(GridFile, "lon");
  MyVector<double> LatArr = NC_Read1Dvariable(GridFile, "lat");
  int nbLon = LonArr.size();
  int nbLat = LatArr.size();
  MyMatrix<double> LON(nbLon, nbLat);
  MyMatrix<double> LAT(nbLon, nbLat);
  MyMatrix<uint8_t> MSK(nbLon, nbLat);
  MyMatrix<double> DEP(nbLon, nbLat);
  MyMatrix<double> ANG(nbLon, nbLat);
  for (int iLon = 0; iLon < nbLon; iLon++)
    for (int iLat = 0; iLat < nbLat; iLat++) {
      LON(iLon, iLat) = LonArr(iLon);
      LAT(iLon, iLat) = LatArr(iLat);
      MSK(iLon, iLat) = 1;
      DEP(iLon, iLat) = 0;
      ANG(iLon, iLat) = 0;
    }
  CoordGridArrayFD GrdArrRho;
  GrdArrRho.LON = LON;
  GrdArrRho.LAT = LAT;
  GrdArrRho.MSK = MSK;
  GrdArrRho.DEP = DEP;
  GrdArrRho.ANG = ANG;
  //
  GridArray GrdArr;
  GrdArr.ARVD.IsAssigned = false;
  GrdArr.ARVD.Zcoordinate = false;
  GrdArr.ModelName = "CFconvention";
  GrdArr.IsFE = 0;
  GrdArr.IsSpherical = true;
  GrdArr.GrdArrRho = GrdArrRho;
  return GrdArr;
}

MyMatrix<double> MatrixSubsample(MyMatrix<double> const &F, int const &splitRow,
                                 int const &splitCol) {
  int nbRow = F.rows();
  int nbCol = F.cols();
  auto GetListIdx = [](int const &nbPos,
                       int const &splitPos) -> std::vector<int> {
    int nbPosRed = nbPos / splitPos;
    double multCoef = double(nbPos - 1) / double(nbPosRed - 1);
    std::vector<int> ListIdx(nbPosRed);
    for (int i = 0; i < nbPosRed; i++) {
      double xPos = double(i) * multCoef;
      int iPos = int(round(xPos));
      int ePos = std::max(0, std::min(nbPos - 1, iPos));
      ListIdx[i] = ePos;
    }
    ListIdx[0] = 0;
    ListIdx[nbPosRed - 1] = nbPos - 1;
    return ListIdx;
  };
  std::vector<int> ListIdxRow = GetListIdx(nbRow, splitRow);
  std::vector<int> ListIdxCol = GetListIdx(nbCol, splitCol);
  int nbRowRed = ListIdxRow.size();
  int nbColRed = ListIdxCol.size();
  MyMatrix<double> Fred(nbRowRed, nbColRed);
  for (int iRowRed = 0; iRowRed < nbRowRed; iRowRed++)
    for (int iColRed = 0; iColRed < nbColRed; iColRed++) {
      int iRow = ListIdxRow[iRowRed];
      int iCol = ListIdxCol[iColRed];
      Fred(iRowRed, iColRed) = F(iRow, iCol);
    }
  return Fred;
}

GridArray CURVILINEAR_GRID_ARRAY(MyMatrix<double> const &LON,
                                 MyMatrix<double> const &LAT) {
  if (LON.rows() != LAT.rows() || LON.cols() != LAT.cols()) {
    std::cerr << "LON and LAT should have same size\n";
    throw TerminalException{1};
  }
  int nbRow = LON.rows();
  int nbCol = LON.cols();
  MyMatrix<uint8_t> MSK(nbRow, nbCol);
  MyMatrix<double> DEP(nbRow, nbCol);
  MyMatrix<double> ANG(nbRow, nbCol);
  for (int iRow = 0; iRow < nbRow; iRow++)
    for (int iCol = 0; iCol < nbCol; iCol++) {
      MSK(iRow, iCol) = 1;
      DEP(iRow, iCol) = 0;
      ANG(iRow, iCol) = 0;
    }
  CoordGridArrayFD GrdArrRho;
  GrdArrRho.LON = LON;
  GrdArrRho.LAT = LAT;
  GrdArrRho.MSK = MSK;
  GrdArrRho.DEP = DEP;
  GrdArrRho.ANG = ANG;
  GrdArrRho.nbWet = nbRow * nbCol;
  //
  GridArray GrdArr;
  GrdArr.ARVD.IsAssigned = false;
  GrdArr.ARVD.Zcoordinate = false;
  GrdArr.ModelName = "CURVILINEAR";
  GrdArr.IsFE = 0;
  GrdArr.IsSpherical = true;
  GrdArr.GrdArrRho = GrdArrRho;
  return GrdArr;
}

GridArray RECTANGULAR_GRID_ARRAY(QuadArray const &eQuad, int const &nbSplitLon,
                                 int const &nbSplitLat) {
  double MinLon = eQuad.MinLon;
  double MinLat = eQuad.MinLat;
  double MaxLon = eQuad.MaxLon;
  double MaxLat = eQuad.MaxLat;
  double deltaLon = (MaxLon - MinLon) / double(nbSplitLon - 1);
  double deltaLat = (MaxLat - MinLat) / double(nbSplitLat - 1);
  std::cerr << "nbSplitLon=" << nbSplitLon << " nbSplitLat=" << nbSplitLat
            << "\n";
  MyMatrix<double> LON(nbSplitLon, nbSplitLat);
  MyMatrix<double> LAT(nbSplitLon, nbSplitLat);
  for (int iLon = 0; iLon < nbSplitLon; iLon++)
    for (int iLat = 0; iLat < nbSplitLat; iLat++) {
      LON(iLon, iLat) = MinLon + deltaLon * iLon;
      LAT(iLon, iLat) = MinLat + deltaLat * iLat;
    }
  GridArray GrdArr = CURVILINEAR_GRID_ARRAY(LON, LAT);
  GrdArr.ARVD.IsAssigned = false;
  GrdArr.ARVD.Zcoordinate = false;
  GrdArr.ModelName = "RECTANGULAR";
  return GrdArr;
}

void CutWorldMap(GridArray &GrdArr) {
  // We cut at -180 - eps
  double eps = 1e-8;
  int nbPoint = GrdArr.GrdArrRho.LON.rows();
  double LonSplit = 0;
  while (true) {
    double MinDist = 2400;
    for (int iPoint = 0; iPoint < nbPoint; iPoint++) {
      double eLon = GrdArr.GrdArrRho.LON(iPoint, 0);
      if (eLon > 0)
        eLon -= 360;
      LonSplit = -180 - eps;
      double dist = fabs(eLon - LonSplit);
      if (dist < MinDist)
        MinDist = dist;
    }
    std::cerr << "eps=" << eps << " MinDist=" << MinDist << "\n";
    if (MinDist > eps / 2)
      break;
    eps *= 2;
  }
  int nbTrig = GrdArr.INE.rows();
  std::vector<int> ListStatus(nbTrig);
  int SumStatus = 0;
  for (int iTrig = 0; iTrig < nbTrig; iTrig++) {
    int i1 = GrdArr.INE(iTrig, 0);
    int i2 = GrdArr.INE(iTrig, 1);
    int i3 = GrdArr.INE(iTrig, 2);
    double eLon1 = GrdArr.GrdArrRho.LON(i1, 0);
    double eLon2 = GrdArr.GrdArrRho.LON(i2, 0);
    double eLon3 = GrdArr.GrdArrRho.LON(i3, 0);
    eLon1 -= LonSplit;
    eLon2 -= LonSplit;
    eLon3 -= LonSplit;
    DifferenceLonRenormalize(eLon1);
    DifferenceLonRenormalize(eLon2);
    DifferenceLonRenormalize(eLon3);
    int eStatus = 1;
    double UpperLimit = 90;
    if (fabs(eLon1) < UpperLimit && fabs(eLon2) < UpperLimit &&
        fabs(eLon3) < UpperLimit) {
      double eProd12 = eLon1 * eLon2;
      double eProd23 = eLon2 * eLon3;
      double eProd31 = eLon3 * eLon1;
      if (eProd12 < 0 || eProd23 < 0 || eProd31 < 0)
        eStatus = 0;
    }
    ListStatus[iTrig] = eStatus;
    SumStatus += eStatus;
  }
  std::cerr << "SumStatus = " << SumStatus << "   nbTrig = " << nbTrig << "\n";
  int nbTrigNew = 0;
  for (int iTrig = 0; iTrig < nbTrig; iTrig++)
    if (ListStatus[iTrig] == 1)
      nbTrigNew++;
  MyMatrix<int> INEnew(nbTrigNew, 3);
  int iTrigNew = 0;
  for (int iTrig = 0; iTrig < nbTrig; iTrig++)
    if (ListStatus[iTrig] == 1) {
      int i1 = GrdArr.INE(iTrig, 0);
      int i2 = GrdArr.INE(iTrig, 1);
      int i3 = GrdArr.INE(iTrig, 2);
      INEnew(iTrigNew, 0) = i1;
      INEnew(iTrigNew, 1) = i2;
      INEnew(iTrigNew, 2) = i3;
      iTrigNew++;
    }
  GrdArr.INE = INEnew;
}

void CUT_HigherLatitude(GridArray &GrdArr, double MinLatCut, double MaxLatCut) {
  int mnp = GrdArr.GrdArrRho.LON.rows();
  int mne = GrdArr.INE.rows();
  std::vector<int> ListStatus(mnp);
  std::vector<int> Index(mnp);
  std::vector<int> RevIndex(mnp);
  int iNodeNew = 0;
  std::vector<int> I_IndexSelectOld;
  if (GrdArr.L_IndexSelect) {
    I_IndexSelectOld = GrdArr.I_IndexSelect;
  } else {
    for (int i = 0; i < mnp; i++)
      I_IndexSelectOld.push_back(i);
  }
  std::vector<int> I_IndexSelect;
  for (int iNode = 0; iNode < mnp; iNode++) {
    double eLat = GrdArr.GrdArrRho.LAT(iNode, 0);
    if (MinLatCut < eLat && eLat < MaxLatCut) {
      ListStatus[iNode] = 1;
      Index[iNode] = iNodeNew;
      RevIndex[iNodeNew] = iNode;
      int iNodeMain = I_IndexSelectOld[iNode];
      I_IndexSelect.push_back(iNodeMain);
      iNodeNew++;
    } else {
      ListStatus[iNode] = 0;
    }
  }
  int nbNodeNew = iNodeNew;
  MyMatrix<double> LONnew(nbNodeNew, 1);
  MyMatrix<double> LATnew(nbNodeNew, 1);
  MyMatrix<double> DEPnew(nbNodeNew, 1);
  MyMatrix<double> ANGnew(nbNodeNew, 1);
  MyMatrix<uint8_t> MSKnew(nbNodeNew, 1);
  const MyMatrix<double> &DEP = GetDEP(GrdArr.GrdArrRho);
  for (int iNodeNewB = 0; iNodeNewB < nbNodeNew; iNodeNewB++) {
    int iNode = RevIndex[iNodeNewB];
    double eLon = GrdArr.GrdArrRho.LON(iNode, 0);
    double eLat = GrdArr.GrdArrRho.LAT(iNode, 0);
    double eDep = DEP(iNode, 0);
    double eAng = GrdArr.GrdArrRho.ANG(iNode, 0);
    int eMsk = GrdArr.GrdArrRho.MSK(iNode, 0);
    LONnew(iNodeNewB, 0) = eLon;
    LATnew(iNodeNewB, 0) = eLat;
    DEPnew(iNodeNewB, 0) = eDep;
    ANGnew(iNodeNewB, 0) = eAng;
    MSKnew(iNodeNewB, 0) = eMsk;
  }
  GrdArr.GrdArrRho.LON = LONnew;
  GrdArr.GrdArrRho.LAT = LATnew;
  GrdArr.GrdArrRho.DEP = DEPnew;
  GrdArr.GrdArrRho.ANG = ANGnew;
  GrdArr.GrdArrRho.MSK = MSKnew;
  int nbTrigNew = 0;
  for (int iTrig = 0; iTrig < mne; iTrig++) {
    int i1 = GrdArr.INE(iTrig, 0);
    int i2 = GrdArr.INE(iTrig, 1);
    int i3 = GrdArr.INE(iTrig, 2);
    if (ListStatus[i1] == 1 && ListStatus[i2] == 1 && ListStatus[i3] == 1)
      nbTrigNew++;
  }
  MyMatrix<int> INEnew(nbTrigNew, 3);
  int iTrigNew = 0;
  for (int iTrig = 0; iTrig < mne; iTrig++) {
    int i1 = GrdArr.INE(iTrig, 0);
    int i2 = GrdArr.INE(iTrig, 1);
    int i3 = GrdArr.INE(iTrig, 2);
    if (ListStatus[i1] == 1 && ListStatus[i2] == 1 && ListStatus[i3] == 1) {
      INEnew(iTrigNew, 0) = Index[i1];
      INEnew(iTrigNew, 1) = Index[i2];
      INEnew(iTrigNew, 2) = Index[i3];
      iTrigNew++;
    }
  }
  GrdArr.INE = INEnew;
  GrdArr.L_IndexSelect = true;
  GrdArr.I_IndexSelect = I_IndexSelect;
}

double GetGridSpacing(GridArray const &GrdArr) {
  int IsFE = GrdArr.IsFE;
  double SumDistKM = 0;
  int SumNb = 0;
  if (IsFE == 1) {
    int nbEle = GrdArr.INE.rows();
    for (int iEle = 0; iEle < nbEle; iEle++)
      for (int i = 0; i < 3; i++) {
        int j = NextIdx(3, i);
        int iNode1 = GrdArr.INE(iEle, i);
        int iNode2 = GrdArr.INE(iEle, j);
        double eLon1 = GrdArr.GrdArrRho.LON(iNode1, 0);
        double eLat1 = GrdArr.GrdArrRho.LAT(iNode1, 0);
        double eLon2 = GrdArr.GrdArrRho.LON(iNode2, 0);
        double eLat2 = GrdArr.GrdArrRho.LAT(iNode2, 0);
        double DistKM = GeodesicDistanceKM(eLon1, eLat1, eLon2, eLat2);
        SumDistKM += DistKM;
        SumNb += 1;
      }
  } else {
    int nbRow = GrdArr.GrdArrRho.LON.rows();
    int nbCol = GrdArr.GrdArrRho.LON.cols();
    for (int iRow = 0; iRow < nbRow; iRow++) {
      for (int iCol = 1; iCol < nbCol; iCol++) {
        double eLon1 = GrdArr.GrdArrRho.LON(iRow, iCol);
        double eLat1 = GrdArr.GrdArrRho.LAT(iRow, iCol);
        double eLon2 = GrdArr.GrdArrRho.LON(iRow, iCol - 1);
        double eLat2 = GrdArr.GrdArrRho.LAT(iRow, iCol - 1);
        double DistKM = GeodesicDistanceKM(eLon1, eLat1, eLon2, eLat2);
        SumDistKM += DistKM;
        SumNb += 1;
      }
    }
    for (int iRow = 1; iRow < nbRow; iRow++) {
      for (int iCol = 0; iCol < nbCol; iCol++) {
        double eLon1 = GrdArr.GrdArrRho.LON(iRow, iCol);
        double eLat1 = GrdArr.GrdArrRho.LAT(iRow, iCol);
        double eLon2 = GrdArr.GrdArrRho.LON(iRow - 1, iCol);
        double eLat2 = GrdArr.GrdArrRho.LAT(iRow - 1, iCol);
        double DistKM = GeodesicDistanceKM(eLon1, eLat1, eLon2, eLat2);
        SumDistKM += DistKM;
        SumNb++;
      }
    }
  }
  double avgDistKM = SumDistKM / double(SumNb);
  return avgDistKM;
}

struct GridSymbolic {
public:
  GridSymbolic() {
    Sphericity = "unset";
    CutWorldMap = false;
    HigherLatitudeCut = false;
    MinLatCut = 0;
    MaxLatCut = 0;
    MinLat = 0;
    MaxLat = 0;
    MinLon = 0;
    MaxLon = 0;
    deltaKM = 0;
  }
  GridSymbolic(std::string _Sphericity, bool _CutWorldMap,
               bool _HigherLatitudeCut, double _MinLatCut, double _MaxLatCut,
               double _MinLat, double _MaxLat, double _MinLon, double _MaxLon,
               double _deltaKM) {
    Sphericity = _Sphericity;
    CutWorldMap = _CutWorldMap;
    HigherLatitudeCut = _HigherLatitudeCut;
    MinLatCut = _MinLatCut;
    MaxLatCut = _MaxLatCut;
    MinLat = _MinLat;
    MaxLat = _MaxLat;
    MinLon = _MinLon;
    MaxLon = _MaxLon;
    deltaKM = _deltaKM;
  }
  std::string Sphericity;
  bool CutWorldMap;
  bool HigherLatitudeCut;
  double MinLatCut;
  double MaxLatCut;
  double MinLat;
  double MaxLat;
  double MinLon;
  double MaxLon;
  double deltaKM;
};

struct TripleModelDesc {
  std::string ModelName;
  std::string GridFile;
  std::string BoundFile;
  std::string HisPrefix;
  GridSymbolic RecGridSymb;
};

std::string GET_GRID_FILE(TripleModelDesc const &eTriple) {
  std::string PreModelName = eTriple.ModelName;
  std::string eModelName = GetKernelModelName(PreModelName);
  std::string HisPrefix = eTriple.HisPrefix;
  if (eModelName == "RECTANGULAR")
    return "irrelevant";
  if (eModelName == "COSMO")
    return HisPrefix + "0001.nc";
  if (eModelName == "WAM")
    return HisPrefix + "0001.nc";
  if (eModelName == "ROMS" || eModelName == "ROMS_IVICA")
    return eTriple.GridFile;
  if (eModelName == "WWM" || eModelName == "WWM_DAILY" ||
      eModelName == "UNRUNOFF")
    return eTriple.GridFile;
  if (eModelName == "SCHISM_SFLUX")
    return eTriple.GridFile;
  if (eModelName == "WRF")
    return HisPrefix + "0001.nc"; // but maybe should be something else
  if (eModelName == "SCHISM_NETCDF_OUT")
    return eTriple.GridFile;
  if (eModelName == "AREG") {
    std::vector<std::string> ListFile =
        FILE_DirectoryFilesSpecificExtension(HisPrefix, "nc");
    if (ListFile.size() == 0) {
      std::cerr << "No netcdf .nc files found in HisPrefix = " << HisPrefix
                << "\n";
      throw TerminalException{1};
    }
    return ListFile[0];
  }
  if (eModelName == "GEOS") {
    std::cerr << "Before FILE_DirectoryFilesSpecificExtension 1\n";
    std::string HisPrefixNaked = FILE_GetDirectoryOfFileName(HisPrefix);
    std::cerr << "HisPrefixaked=" << HisPrefixNaked << "\n";
    std::vector<std::string> ListFile =
        FILE_DirectoryFilesSpecificExtension(HisPrefixNaked, "nc4");
    std::cerr << " After FILE_DirectoryFilesSpecificExtension 1\n";
    if (ListFile.size() == 0) {
      std::cerr << "No netcdf .nc4 files found in HisPrefix = " << HisPrefix
                << "\n";
      throw TerminalException{1};
    }
    return ListFile[0];
  }
  if (eModelName == "WW3") {
    std::string ThePrefix = HisPrefix + "*";
    std::vector<std::string> ListFile = ls_operation(ThePrefix);
    return ListFile[0];
  }
  if (eModelName == "NEMO") {
    std::vector<std::string> ListFile =
        FILE_DirectoryMatchingPrefixExtension(HisPrefix, "nc");
    std::cerr << "NEMO : |ListFile|=" << ListFile.size() << "\n";
    for (auto &eFile : ListFile) {
      NEMO_vars nemo_vars = ReadNEMO_vars(eFile);
      if (nemo_vars.List3D_vars.size() > 0) {
        std::cerr << "Finding NEMO grid file eFile=" << eFile << "\n";
        return eFile;
      }
    }
    std::cerr
        << "We failed to find the matching file with a 3D variable for NEMO\n";
    throw TerminalException{1};
  }
  if (eModelName == "HYCOM") {
    std::cerr << "HisPrefix=" << HisPrefix << "\n";
    std::vector<std::string> ListFile =
        FILE_DirectoryFilesSpecificExtension_Gen(HisPrefix, "nc");
    if (ListFile.size() > 0) {
      std::cerr << "ListFile[0]=" << ListFile[0] << "\n";
      return ListFile[0];
    }
    std::cerr << "We failed to find the matching file for HYCOM\n";
    throw TerminalException{1};
  }
  if (eModelName == "GRIB_DWD" || eModelName == "GRIB_GFS" ||
      eModelName == "GRIB_ECMWF" || eModelName == "GRIB_COSMO" ||
      eModelName == "GRIB_ALADIN" || eModelName == "GRIB_IFS") {
    std::vector<std::string> ListFile =
        FILE_DirectoryFilesSpecificExtension(HisPrefix, "grb");
    if (ListFile.size() == 0) {
      std::cerr << "The list of files is empty\n";
      throw TerminalException{1};
    }
    return ListFile[0];
  }
  if (eModelName == "GRIB_WAM_FORT30") {
    std::string eFile = HisPrefix;
    if (!IsExistingFile(eFile)) {
      std::cerr << "The file eFile = " << eFile << " is missing\n";
      std::cerr << "It serves as grid and should be put in HisPrefix\n";
      throw TerminalException{1};
    }
    return eFile;
  }
  std::cerr << "Error in GET_GRID_FILE\n";
  std::cerr << "Did not find the matching model for the grid\n";
  std::cerr << "Please correct\n";
  throw TerminalException{1};
}

void WriteUnstructuredGrid_UGRID_CF_NC(std::string const &GridFile,
                                       GridArray const &GrdArr) {
  if (!GrdArr.IsSpherical) {
    std::cerr << "We shuld have a spherical grid\n";
    throw TerminalException{1};
  }
  if (!FILE_IsFileMakeable(GridFile)) {
    std::cerr << "Request to create file GridFile=" << GridFile << "\n";
    std::cerr << "but the directory does not exist\n";
    throw TerminalException{1};
  }
  netCDF::NcFile dataFile(GridFile, netCDF::NcFile::replace,
                          netCDF::NcFile::nc4);
  int mnp = GrdArr.GrdArrRho.LON.rows();
  int nbFace = GrdArr.INE.rows();
  MyMatrix<int> LEdge = GetEdgeSet(GrdArr.INE, mnp);
  int nbEdge = LEdge.rows();
  int nMaxMesh2_face_nodes = 3;
  netCDF::NcDim eDimTwo = dataFile.addDim("two", 2);
  netCDF::NcDim eDimMnp = dataFile.addDim("nMesh2_node", mnp);
  netCDF::NcDim eDimMne = dataFile.addDim("nMesh2_face", nbFace);
  netCDF::NcDim eDimNbedge = dataFile.addDim("nMesh2_edge", nbEdge);
  netCDF::NcDim eDimMax =
      dataFile.addDim("nMaxMesh2_face_nodes", nMaxMesh2_face_nodes);
  std::vector<std::string> ListDim_Nodes{"nMesh2_node"};
  std::vector<std::string> ListDim_Edges{"nMesh2_edge"};
  std::vector<std::string> ListDim_Faces{"nMesh2_face"};
  std::vector<std::string> ListDim_Conn1{"nMesh2_edge", "two"};
  std::vector<std::string> ListDim_Conn2{"nMesh2_face", "nMaxMesh2_face_nodes"};
  //
  // Nodes coordinates
  //
  int Aint_1653[1], Aint_1652[1], Aint_zero[1], Aint_m999[1];
  Aint_1653[0] = 1653;
  Aint_1652[0] = 1652;
  Aint_zero[0] = 0;
  Aint_m999[0] = -999;
  {
    std::vector<double> Anode(mnp);
    //
    netCDF::NcVar eVar_lon_node =
        dataFile.addVar("Mesh2_node_lon", "double", ListDim_Nodes);
    eVar_lon_node.putAtt("long_name", "longitude of 2D mesh nodes");
    eVar_lon_node.putAtt("units", "degrees_east");
    eVar_lon_node.putAtt("name_id", netCDF::NcType::nc_INT, 1, Aint_1653);
    eVar_lon_node.putAtt("standard_name", "longitude");
    for (int ip = 0; ip < mnp; ip++)
      Anode[ip] = GrdArr.GrdArrRho.LON(ip, 0);
    eVar_lon_node.putVar(Anode.data());
    //
    netCDF::NcVar eVar_lat_node =
        dataFile.addVar("Mesh2_node_lat", "double", ListDim_Nodes);
    eVar_lat_node.putAtt("long_name", "latitude of 2D mesh nodes");
    eVar_lat_node.putAtt("units", "degrees_north");
    eVar_lat_node.putAtt("name_id", netCDF::NcType::nc_INT, 1, Aint_1652);
    eVar_lat_node.putAtt("standard_name", "latitude");
    for (int ip = 0; ip < mnp; ip++)
      Anode[ip] = GrdArr.GrdArrRho.LAT(ip, 0);
    eVar_lat_node.putVar(Anode.data());
  }
  //
  // Edges coordinates
  //
  {
    std::vector<double> Aedge(nbEdge);
    //
    MyMatrix<double> LonLatEdge(nbEdge, 2);
    for (int iEdge = 0; iEdge < nbEdge; iEdge++) {
      int iNode1 = LEdge(iEdge, 0);
      int iNode2 = LEdge(iEdge, 1);
      double eLon =
          (GrdArr.GrdArrRho.LON(iNode1, 0) + GrdArr.GrdArrRho.LON(iNode2, 0)) /
          double(2);
      double eLat =
          (GrdArr.GrdArrRho.LAT(iNode1, 0) + GrdArr.GrdArrRho.LAT(iNode2, 0)) /
          double(2);
      LonLatEdge(iEdge, 0) = eLon;
      LonLatEdge(iEdge, 1) = eLat;
    }
    //
    netCDF::NcVar eVar_lon_edge =
        dataFile.addVar("Mesh2_edge_lon", "double", ListDim_Edges);
    eVar_lon_edge.putAtt("long_name", "longitude of 2D mesh edges (center)");
    eVar_lon_edge.putAtt("units", "degrees_east");
    eVar_lon_edge.putAtt("name_id", netCDF::NcType::nc_INT, 1, Aint_1653);
    eVar_lon_edge.putAtt("bounds", "Mesh2_edge_lon_bnd");
    eVar_lon_edge.putAtt("standard_name", "longitude");
    for (int iEdge = 0; iEdge < nbEdge; iEdge++)
      Aedge[iEdge] = LonLatEdge(iEdge, 0);
    eVar_lon_edge.putVar(Aedge.data());
    //
    netCDF::NcVar eVar_lat_edge =
        dataFile.addVar("Mesh2_edge_lat", "double", ListDim_Edges);
    eVar_lat_edge.putAtt("long_name", "latitude of 2D mesh edges (center)");
    eVar_lat_edge.putAtt("units", "degrees_north");
    eVar_lat_edge.putAtt("name_id", netCDF::NcType::nc_INT, 1, Aint_1652);
    eVar_lat_edge.putAtt("bounds", "Mesh2_edge_lat_bnd");
    eVar_lat_edge.putAtt("standard_name", "latitude");
    for (int iEdge = 0; iEdge < nbEdge; iEdge++)
      Aedge[iEdge] = LonLatEdge(iEdge, 1);
    eVar_lat_edge.putVar(Aedge.data());
  }
  //
  // Faces coordinates
  //
  {
    std::vector<double> Aface(nbFace);
    //
    MyMatrix<double> LonLatFace(nbFace, 4);
    for (int iFace = 0; iFace < nbFace; iFace++) {
      int iNode1 = GrdArr.INE(iFace, 0);
      int iNode2 = GrdArr.INE(iFace, 1);
      int iNode3 = GrdArr.INE(iFace, 2);
      double eLonCG =
          (GrdArr.GrdArrRho.LON(iNode1, 0) + GrdArr.GrdArrRho.LON(iNode2, 0) +
           GrdArr.GrdArrRho.LON(iNode3, 0)) /
          double(2);
      double eLatCG =
          (GrdArr.GrdArrRho.LAT(iNode1, 0) + GrdArr.GrdArrRho.LAT(iNode2, 0) +
           GrdArr.GrdArrRho.LAT(iNode3, 0)) /
          double(2);
      LonLatFace(iFace, 0) = eLonCG; // Center of gravity
      LonLatFace(iFace, 1) = eLatCG; // Center of gravity
      LonLatFace(iFace, 2) = eLonCG; // Circumcenter
      LonLatFace(iFace, 3) = eLatCG; // Circumcenter
    }
    //
    netCDF::NcVar eVar_lon_faceCG =
        dataFile.addVar("Mesh2_face_lon", "double", ListDim_Faces);
    eVar_lon_faceCG.putAtt("long_name",
                           "longitude of 2D mesh nodes (center of gravity)");
    eVar_lon_faceCG.putAtt("units", "degrees_east");
    eVar_lon_faceCG.putAtt("name_id", netCDF::NcType::nc_INT, 1, Aint_1653);
    eVar_lon_faceCG.putAtt("bounds", "Mesh2_face_lon_bnd");
    eVar_lon_faceCG.putAtt("standard_name", "longitude");
    for (int iFace = 0; iFace < nbFace; iFace++)
      Aface[iFace] = LonLatFace(iFace, 0);
    eVar_lon_faceCG.putVar(Aface.data());
    //
    netCDF::NcVar eVar_lat_faceCG =
        dataFile.addVar("Mesh2_face_lon", "double", ListDim_Faces);
    eVar_lat_faceCG.putAtt("long_name",
                           "latitude of 2D mesh nodes (center of gravity)");
    eVar_lat_faceCG.putAtt("units", "degrees_north");
    eVar_lat_faceCG.putAtt("name_id", netCDF::NcType::nc_INT, 1, Aint_1652);
    eVar_lat_faceCG.putAtt("bounds", "Mesh2_face_lon_bnd");
    eVar_lat_faceCG.putAtt("standard_name", "latitude");
    for (int iFace = 0; iFace < nbFace; iFace++)
      Aface[iFace] = LonLatFace(iFace, 1);
    eVar_lat_faceCG.putVar(Aface.data());
    //
    netCDF::NcVar eVar_lon_faceCC =
        dataFile.addVar("Mesh2_face_lon", "double", ListDim_Faces);
    eVar_lon_faceCC.putAtt("long_name",
                           "longitude of 2D mesh nodes (circumcenter)");
    eVar_lon_faceCC.putAtt("units", "degrees_east");
    eVar_lon_faceCC.putAtt("name_id", netCDF::NcType::nc_INT, 1, Aint_1653);
    eVar_lon_faceCC.putAtt("standard_name", "longitude");
    for (int iFace = 0; iFace < nbFace; iFace++)
      Aface[iFace] = LonLatFace(iFace, 2);
    eVar_lon_faceCC.putVar(Aface.data());
    //
    netCDF::NcVar eVar_lat_faceCC =
        dataFile.addVar("Mesh2_face_lon", "double", ListDim_Faces);
    eVar_lat_faceCC.putAtt("long_name",
                           "latitude of 2D mesh nodes (center of gravity)");
    eVar_lat_faceCC.putAtt("units", "degrees_north");
    eVar_lat_faceCC.putAtt("name_id", netCDF::NcType::nc_INT, 1, Aint_1652);
    eVar_lat_faceCC.putAtt("standard_name", "latitude");
    for (int iFace = 0; iFace < nbFace; iFace++)
      Aface[iFace] = LonLatFace(iFace, 3);
    eVar_lat_faceCC.putVar(Aface.data());
  }
  //
  // Edge-node connectivity
  //
  {
    std::vector<int> Aconn1(nbEdge * 2);
    //
    int idx = 0;
    for (int iEdge = 0; iEdge < nbEdge; iEdge++)
      for (int i = 0; i < 2; i++) {
        Aconn1[idx] = LEdge(iEdge, i);
        idx++;
      }
    //
    netCDF::NcVar eVar_edge_node =
        dataFile.addVar("Mesh2_edge_nodes", "int", ListDim_Conn1);
    eVar_edge_node.putAtt("long_name",
                          "list of nodes for all edges, start node - end node");
    eVar_edge_node.putAtt("cf_role", "edge_node_connectivity");
    eVar_edge_node.putAtt("start_index", netCDF::NcType::nc_INT, 1, Aint_zero);
    eVar_edge_node.putVar(Aconn1.data());
  }
  //
  // Edge-face connectivity
  //
  {
    std::vector<int> Aconn1(nbEdge * 2);
    //
    MyMatrix<int> FaceEdgeConn =
        GetFaceEdgeConnectivity(mnp, LEdge, GrdArr.INE);
    int idx = 0;
    for (int iEdge = 0; iEdge < nbEdge; iEdge++)
      for (int i = 0; i < 2; i++) {
        Aconn1[idx] = LEdge(iEdge, i);
        idx++;
      }
    //
    netCDF::NcVar eVar_edge_face =
        dataFile.addVar("Mesh2_edge_faces", "int", ListDim_Conn1);
    eVar_edge_face.putAtt("long_name",
                          "list of (adjacent) faces (polygons) for all edges - "
                          "left and right neigbour");
    eVar_edge_face.putAtt("cf_role", "edge_face_connectivity");
    eVar_edge_face.putAtt("start_index", netCDF::NcType::nc_INT, 1, Aint_zero);
    eVar_edge_face.putAtt("_FillValue", netCDF::NcType::nc_INT, 1, Aint_m999);
    eVar_edge_face.putVar(Aconn1.data());
  }
  //
  // Face-node connectivity
  //
  {
    std::vector<int> Aconn2(nbFace * 2);
    //
    int idx = 0;
    for (int iFace = 0; iFace < nbFace; iFace++)
      for (int i = 0; i < 2; i++) {
        Aconn2[idx] = GrdArr.INE(iFace, i);
        idx++;
      }
    //
    netCDF::NcVar eVar_face_node =
        dataFile.addVar("Mesh2_edge_faces", "int", ListDim_Conn2);
    eVar_face_node.putAtt(
        "long_name",
        "list of nodes for all faces (polygons), counterclockwise");
    eVar_face_node.putAtt("cf_role", "face_node_connectivity");
    eVar_face_node.putAtt("start_index", netCDF::NcType::nc_INT, 1, Aint_zero);
    eVar_face_node.putAtt("_FillValue", netCDF::NcType::nc_INT, 1, Aint_m999);
    eVar_face_node.putVar(Aconn2.data());
  }
}

void WriteUnstructuredGrid_GR3(std::string const &GridFile,
                               GridArray const &GrdArr) {
  std::ofstream os(GridFile);
  os << std::fixed;
  os << std::setprecision(9);
  os << GridFile << "\n";
  int mnp = GrdArr.GrdArrRho.LON.rows();
  int mne = GrdArr.INE.rows();
  const MyMatrix<double> &DEP = GetDEP(GrdArr.GrdArrRho);
  os << mne << " " << mnp << "\n";
  for (int ip = 0; ip < mnp; ip++) {
    double lon = GrdArr.GrdArrRho.LON(ip);
    double lat = GrdArr.GrdArrRho.LAT(ip);
    double dep = DEP(ip);
    int idx = ip + 1;
    os << idx << " " << lon << " " << lat << " " << dep << "\n";
  }
  for (int ie = 0; ie < mne; ie++) {
    int ip1 = GrdArr.INE(ie, 0) + 1;
    int ip2 = GrdArr.INE(ie, 1) + 1;
    int ip3 = GrdArr.INE(ie, 2) + 1;
    int idx = ie + 1;
    os << idx << " 3 " << ip1 << " " << ip2 << " " << ip3 << "\n";
  }
}

GridArray WWM_ReadGridFile_Ricchiuto_grd(std::string const &GridFile) {
  std::ifstream IN(GridFile);
  GridArray GrdArr;
  GrdArr.ARVD.IsAssigned = false;
  GrdArr.ARVD.Zcoordinate = false;
  //
  std::string FirstLine;
  std::getline(IN, FirstLine);
  //
  std::string lineDATA;
  std::getline(IN, lineDATA);
  std::vector<std::string> LStr = STRING_Split(lineDATA, " ");
  int mne, mnp;
  std::istringstream(LStr[1]) >> mne;
  std::istringstream(LStr[2]) >> mnp;
  //
  MyMatrix<int> INE(mne, 3);
  for (int ie = 0; ie < mne; ie++) {
    std::string lineINE;
    std::getline(IN, lineINE);
    std::vector<std::string> LStrB = STRING_Split(lineINE, " ");
    for (int i = 0; i < 3; i++) {
      int IP;
      std::istringstream(LStrB[i]) >> IP;
      INE(ie, i) = IP;
    }
  }
  GrdArr.INE = INE;
  GrdArr.IsFE = true;
  //
  MyMatrix<double> LON(mnp, 1), LAT(mnp, 1), DEP(mnp, 1);
  for (int ip = 0; ip < mnp; ip++) {
    std::string lineXYD;
    std::getline(IN, lineXYD);
    std::vector<std::string> LStrC = STRING_Split(lineXYD, " ");
    double eLon, eLat, eDep;
    std::istringstream(LStrC[0]) >> eLon;
    std::istringstream(LStrC[1]) >> eLat;
    std::istringstream(LStrC[2]) >> eDep;
    LON(ip, 0) = eLon;
    LAT(ip, 0) = eLat;
    DEP(ip, 0) = eDep;
  }
  GrdArr.GrdArrRho.LON = LON;
  GrdArr.GrdArrRho.LAT = LAT;
  GrdArr.GrdArrRho.DEP = DEP;
  GrdArr.IsSpherical = true;
  return GrdArr;
}

void WriteUnstructuredGrid_Ricchiuto_GRD(std::string const &GridFile,
                                         GridArray const &GrdArr) {
  std::ofstream os(GridFile);
  os << std::fixed;
  os << std::setprecision(9);
  os << GridFile << "\n";
  int NN = GrdArr.GrdArrRho.LON.rows();
  int NE = GrdArr.INE.rows();
  int NBF = 0;
  const MyMatrix<double> &DEP = GetDEP(GrdArr.GrdArrRho);
  os << " 2 " << NE << " " << NN << " " << NBF << "\n";
  for (int ie = 0; ie < NE; ie++) {
    int ip1 = GrdArr.INE(ie, 0);
    int ip2 = GrdArr.INE(ie, 1);
    int ip3 = GrdArr.INE(ie, 2);
    os << ip1 << " " << ip2 << " " << ip3 << "\n";
  }
  for (int ip = 0; ip < NN; ip++) {
    double lon = GrdArr.GrdArrRho.LON(ip);
    double lat = GrdArr.GrdArrRho.LAT(ip);
    double dep = DEP(ip);
    os << lon << " " << lat << " " << dep << "\n";
  }
  for (int ib = 0; ib < NBF; ib++) {
    /* Write something here */
  }
}

void WriteUnstructuredGrid_DAT(std::string const &GridFile,
                               GridArray const &GrdArr) {
  std::ofstream os(GridFile);
  os << std::fixed;
  os << std::setprecision(9);
  int mnp = GrdArr.GrdArrRho.LON.rows();
  int mne = GrdArr.INE.rows();
  os << "C system.dat, erzeugt von xf am sometime\n";
  os << "C Anzahl der Randknoten:\n";
  os << 0 << "\n";
  os << "C Anzahl der Gebietsknoten:\n";
  os << mnp << "\n";
  os << "C Koordinaten und Skalarwerte der Knoten\n";
  os << "C --------------------------------------\n";
  os << "C Zuerst die Randknoten  (Anzahl s.o.),\n";
  os << "C dann die Gebietsknoten (Anzahl s.o.).\n";
  os << "C ------------+-------------+-------------+---------------\n";
  os << "C     Nr.     |  x-Koord.   |   y-Koord.  |   Skalarwert\n";
  os << "C ------------+-------------+-------------+---------------\n";
  const MyMatrix<double> &DEP = GetDEP(GrdArr.GrdArrRho);
  for (int ip = 0; ip < mnp; ip++) {
    double lon = GrdArr.GrdArrRho.LON(ip);
    double lat = GrdArr.GrdArrRho.LAT(ip);
    double dep = DEP(ip);
    int idx = ip;
    os << idx << " " << lon << " " << lat << " " << dep << "\n";
  }
  os << "C ------------------------------------------------------------\n";
  os << "C Anzahl der Elemente:\n";
  os << mne << "\n";
  os << "C Elementverzeichnis\n";
  os << "C ------------------------------------------------------------\n";
  os << "C    Knoten i  Knoten j  Knoten k   Kennung     Nr.\n";
  for (int ie = 0; ie < mne; ie++) {
    int ip1 = GrdArr.INE(ie, 0);
    int ip2 = GrdArr.INE(ie, 1);
    int ip3 = GrdArr.INE(ie, 2);
    int idx = ie;
    os << ip1 << " " << ip2 << " " << ip3 << " 0 " << idx << "\n";
  }
  os << "C ------------------------------------------------------------\n";
}

void WriteUnstructuredGrid_NC(std::string const &GridFile,
                              GridArray const &GrdArr) {
  if (!FILE_IsFileMakeable(GridFile)) {
    std::cerr << "Request to create file GridFile=" << GridFile << "\n";
    std::cerr << "but the directory does not exist\n";
    throw TerminalException{1};
  }
  netCDF::NcFile dataFile(GridFile, netCDF::NcFile::replace,
                          netCDF::NcFile::nc4);
  int mnp = GrdArr.GrdArrRho.LON.rows();
  int mne = GrdArr.INE.rows();
  netCDF::NcDim eDimOne = dataFile.addDim("one", 1);
  netCDF::NcDim eDimTwo = dataFile.addDim("two", 2);
  netCDF::NcDim eDimThree = dataFile.addDim("three", 3);
  netCDF::NcDim eDimMnp = dataFile.addDim("mnp", mnp);
  netCDF::NcDim eDimMne = dataFile.addDim("mne", mne);
  std::vector<std::string> ListDimLL{"mnp"};
  std::vector<std::string> ListDimINE{"mne", "three"};
  //
  std::string strLON, strLAT;
  if (GrdArr.IsSpherical) {
    strLON = "lon";
    strLAT = "lat";
  } else {
    strLON = "x";
    strLAT = "y";
  }
  auto putVARmnp = [&](std::string const &name,
                       MyMatrix<double> const &VAR) -> void {
    netCDF::NcVar eVar = dataFile.addVar(name, "double", ListDimLL);
    std::vector<double> A(mnp);
    for (int ip = 0; ip < mnp; ip++)
      A[ip] = VAR(ip, 0);
    eVar.putVar(A.data());
  };
  putVARmnp(strLON, GrdArr.GrdArrRho.LON);
  putVARmnp(strLAT, GrdArr.GrdArrRho.LAT);
  putVARmnp("depth", GetDEP(GrdArr.GrdArrRho));
  //
  if (GrdArr.IOBP.size() != mnp) {
    std::cerr << "We have |GrdArr.IOBP.size()|=" << GrdArr.IOBP.size() << "\n";
    std::cerr << "Which is distinct from mnp=" << mnp << "\n";
    std::cerr << "Most likely you forgot to put the input file\n";
    throw TerminalException{1};
  }
  netCDF::NcVar eVar = dataFile.addVar("IOBP", "int", ListDimLL);
  std::vector<int> A(mnp);
  for (int ip = 0; ip < mnp; ip++)
    A[ip] = GrdArr.IOBP(ip, 0);
  eVar.putVar(A.data());
  //
  netCDF::NcVar eVarIne = dataFile.addVar("ele", "int", ListDimINE);
  std::vector<int> Aine(mne * 3);
  int idx = 0;
  for (int ie = 0; ie < mne; ie++)
    for (int j = 0; j < 3; j++) {
      Aine[idx] = GrdArr.INE(ie, j) + 1;
      idx++;
    }
  eVarIne.putVar(Aine.data());
  //
  MyMatrix<int> LEdge = GetEdgeSet(GrdArr.INE, mnp);
  int nbEdge = LEdge.rows();
  netCDF::NcDim eDimNbEdge = dataFile.addDim("nbedge", nbEdge);
  std::vector<std::string> ListDimEDGE{"nbedge", "two"};
  std::vector<int> Aedge(nbEdge * 2);
  netCDF::NcVar eVarEdge = dataFile.addVar("edges", "int", ListDimEDGE);
  int idx2 = 0;
  for (int iedge = 0; iedge < nbEdge; iedge++)
    for (int j = 0; j < 2; j++) {
      Aedge[idx2] = LEdge(iedge, j) + 1;
      idx2++;
    }
  eVarEdge.putVar(Aedge.data());
}

void WriteUnstructuredGrid(std::string const &GridFile,
                           GridArray const &GrdArr) {
  std::string eExtension = FILE_GetExtension(GridFile);
  std::cerr << "WriteUnstructuredGrid  with  eExtension=" << eExtension << "\n";
  if (eExtension == "gr3") {
    WriteUnstructuredGrid_GR3(GridFile, GrdArr);
    return;
  }
  if (eExtension == "dat") {
    WriteUnstructuredGrid_DAT(GridFile, GrdArr);
    return;
  }
  if (eExtension == "grd") {
    WriteUnstructuredGrid_Ricchiuto_GRD(GridFile, GrdArr);
    return;
  }
  if (eExtension == "nc") {
    WriteUnstructuredGrid_NC(GridFile, GrdArr);
    return;
  }
  if (eExtension == "msh") {
    std::cerr << "Before call to WriteGridFile_msh\n";
    WriteGridFile_msh(GridFile, GrdArr);
    std::cerr << " After call to WriteGridFile_msh\n";
    return;
  }
  std::cerr << "eExtension = " << eExtension << "\n";
  std::cerr << "Programmed extensions are (gr3 , dat , grd , nc , msh)\n";
  throw TerminalException{1};
}

void WriteGrid(std::string const &GridFile, GridArray const &GrdArr) {
  if (GrdArr.IsFE == 1)
    return WriteUnstructuredGrid(GridFile, GrdArr);
  //
  if (GrdArr.ModelName == "ROMS") {
    NETCDF_Write2Dvariable(GridFile, "h", GetDEP(GrdArr.GrdArrRho));
    return;
  }
  //
  std::cerr << "Missing code for the model in the WriteGrid\n";
  throw TerminalException{1};
}

void WriteUnstructuredGridTot(std::string const &GridFile,
                              std::string const &BoundFile,
                              GridArray const &GrdArr) {
  std::string eExtension = FILE_GetExtension(GridFile);
  GridArray GrdArrCopy = GrdArr;
  int nbNode = GrdArr.IOBP.size();
  MyMatrix<double> IOBPmat(nbNode, 1);
  for (int ip = 0; ip < nbNode; ip++) {
    double eVal = double(GrdArr.IOBP(ip));
    IOBPmat(ip, 0) = eVal;
  }
  GrdArrCopy.GrdArrRho.DEP = IOBPmat;
  if (eExtension == "gr3") {
    WriteUnstructuredGrid_GR3(GridFile, GrdArr);
    WriteUnstructuredGrid_GR3(BoundFile, GrdArrCopy);
    return;
  }
  if (eExtension == "dat") {
    WriteUnstructuredGrid_DAT(GridFile, GrdArr);
    WriteUnstructuredGrid_DAT(BoundFile, GrdArrCopy);
    return;
  }
  if (eExtension == "nc") {
    WriteUnstructuredGrid_NC(GridFile, GrdArr);
    return;
  }
  if (eExtension == "msh") {
    WriteGridFile_msh(GridFile, GrdArr);
    return;
  }
  std::cerr << "eExtension = " << eExtension << "\n";
  std::cerr << "Programmed extensions are .gr3 , .dat and .nc only\n";
  throw TerminalException{1};
}

GridArray ReadUnstructuredGrid(std::string const &GridFile,
                               std::string const &BoundFile) {
  std::string eExtension = FILE_GetExtension(GridFile);
  if (!IsExistingFile(GridFile)) {
    std::cerr << "The file GridFile = " << GridFile << " is not existent\n";
    throw TerminalException{1};
  }
  //    std::cerr << "eExtension=" << eExtension << "\n";
  if (eExtension == "gr3" || eExtension == "ll") {
    GridArray GrdArr = WWM_ReadGridFile_gr3(GridFile);
    if (BoundFile != "unset") {
      std::string fExtension = FILE_GetExtension(BoundFile);
      if (fExtension != "gr3" && fExtension != "ll") {
        std::cerr << "Error in ReasUnstructuredGrid\n";
        std::cerr << "GridFile = " << GridFile << " so we are in gr3 case\n";
        std::cerr << "But BoundFile = " << BoundFile << "\n";
        std::cerr << "which is not ending by gr3 or ll\n";
        throw TerminalException{1};
      }
      std::cerr << "BoundFile=" << BoundFile << "\n";
      MyVector<int> eVect = WWM_ReadBoundFile_gr3(BoundFile);
      if (eVect.size() != GrdArr.GrdArrRho.LON.size()) {
        std::cerr << "not same number of vertices between grid file and "
                     "boundary file\n";
        std::cerr << "nbVert(grid)=" << GrdArr.GrdArrRho.LON.size() << "\n";
        std::cerr << "nbVert(bound)=" << eVect.size() << "\n";
        throw TerminalException{1};
      }
      GrdArr.IOBP = eVect;
    }
    return GrdArr;
  }
  if (eExtension == "dat") {
    std::cerr << "Before WWM_ReadGridFile_DAT\n";
    GridArray GrdArr = WWM_ReadGridFile_DAT(GridFile);
    std::cerr << "After  WWM_ReadGridFile_DAT\n";
    if (BoundFile != "unset") {
      std::string fExtension = FILE_GetExtension(BoundFile);
      if (fExtension != "dat") {
        std::cerr << "Error in ReasUnstructuredGrid\n";
        std::cerr << "GridFile = " << GridFile
                  << " so we are in dat/xfn case\n";
        std::cerr << "But BoundFile = " << BoundFile << "\n";
        std::cerr << "which is not ending by dat\n";
        throw TerminalException{1};
      }
      MyVector<int> eVect = WWM_ReadBoundFile_DAT(BoundFile);
      if (eVect.size() != GrdArr.GrdArrRho.LON.size()) {
        std::cerr << "not same number of vertices between grid file and "
                     "boundary file\n";
        std::cerr << "nbVert(grid)=" << GrdArr.GrdArrRho.LON.size() << "\n";
        std::cerr << "nbVert(bound)=" << eVect.size() << "\n";
        throw TerminalException{1};
      }
      GrdArr.IOBP = eVect;
    }
    return GrdArr;
  }
  if (eExtension == "grd") {
    GridArray GrdArr = WWM_ReadGridFile_Ricchiuto_grd(GridFile);
    return GrdArr;
  }
  if (eExtension == "obj") {
    std::cerr << "Before WWM_ReadGridFile_obj\n";
    GridArray GrdArr = WWM_ReadGridFile_obj(GridFile);
    std::cerr << "After  WWM_ReadGridFile_obj\n";
    return GrdArr;
  }
  if (eExtension == "nc")
    return WWM_ReadGridFile_netcdf(GridFile);
  if (eExtension == "msh")
    return WWM_ReadGridFile_msh(GridFile);
  std::cerr << "Error in reading grid for WWM\n";
  std::cerr << "Supported formats: gr3, ll, dat, grd, obj, nc, msh\n";
  std::cerr << "We did not find the right kind\n";
  throw TerminalException{1};
}

GridArray PRE_RETRIEVE_GRID_ARRAY(TripleModelDesc const &eTriple) {
  std::string PreModelName = eTriple.ModelName;
  std::string eModelName = GetKernelModelName(PreModelName);
  CHECK_Model_Allowedness(eModelName);
  std::string GridFile = GET_GRID_FILE(eTriple);
  std::cerr << "PRE_RETRIEVE_GRID_ARRAY : eModelName=" << eModelName << "\n";
  if (eModelName == "RECTANGULAR") {
    QuadArray eQuad;
    eQuad.MaxLat = eTriple.RecGridSymb.MaxLat;
    eQuad.MinLat = eTriple.RecGridSymb.MinLat;
    eQuad.MaxLon = eTriple.RecGridSymb.MaxLon;
    eQuad.MinLon = eTriple.RecGridSymb.MinLon;
    std::cerr << "Lat (min/max)=" << eQuad.MinLat << " / " << eQuad.MaxLat
              << "\n";
    std::cerr << "Lon (min/max)=" << eQuad.MinLon << " / " << eQuad.MaxLon
              << "\n";
    double deltaKM = eTriple.RecGridSymb.deltaKM;
    std::cerr << "deltaKM=" << deltaKM << "\n";
    double distLON = GeodesicDistanceKM(eQuad.MinLon, eQuad.MinLat,
                                        eQuad.MaxLon, eQuad.MinLat);
    double distLAT = GeodesicDistanceKM(eQuad.MinLon, eQuad.MinLat,
                                        eQuad.MinLon, eQuad.MaxLat);
    double nbLON = int(distLON / deltaKM);
    double nbLAT = int(distLAT / deltaKM);
    return RECTANGULAR_GRID_ARRAY(eQuad, nbLON, nbLAT);
  }
  if (eModelName == "COSMO")
    return NC_ReadCosmoWamStructGridFile(GridFile, "atm");
  if (eModelName == "WAM")
    return NC_ReadWamGridFile(GridFile);
  if (eModelName == "ROMS" || eModelName == "ROMS_IVICA")
    return NC_ReadRomsGridFile(GridFile);
  if (eModelName == "WRF")
    return NC_ReadWrfGridFile(GridFile);
  if (eModelName == "NEMO")
    return NC_ReadNemoGridFile(GridFile);
  if (eModelName == "HYCOM")
    return NC_ReadHycomGridFile(GridFile);
  if (eModelName == "WWM" || eModelName == "WWM_DAILY") {
    std::string BoundFile = eTriple.BoundFile;
    std::cerr << "GridFile=" << GridFile << "\n";
    std::cerr << "BoundFile=" << BoundFile << "\n";
    GridArray GrdArr = ReadUnstructuredGrid(GridFile, BoundFile);
    GrdArr.ModelName = "WWM";
    return GrdArr;
  }
  if (eModelName == "UNRUNOFF") {
    std::string BoundFile = eTriple.BoundFile;
    GridArray GrdArr = ReadUnstructuredGrid(GridFile, BoundFile);
    GrdArr.ModelName = "UNRUNOFF";
    return GrdArr;
  }
  if (eModelName == "WW3")
    return NC_ReadWW3_GridFile(GridFile);
  if (eModelName == "SCHISM_SFLUX")
    return NC_ReadSCHISM_sflux_grid(GridFile);
  if (eModelName == "SCHISM_NETCDF_OUT") {
    std::string BoundFile = eTriple.BoundFile;
    GridArray GrdArr = ReadUnstructuredGrid(GridFile, BoundFile);
    GrdArr.ModelName = "SCHISM_NETCDF_OUT";
    return GrdArr;
  }
  if (eModelName == "AREG") {
    std::string HisPrefix = eTriple.HisPrefix;
    std::vector<std::string> ListFile =
        FILE_DirectoryFilesSpecificExtension(HisPrefix, "grb");
    if (ListFile.size() == 0) {
      std::cerr << "The list of files for AREG retrieval of grid is empty\n";
      throw TerminalException{1};
    }
    std::string eFileName = ListFile[0];
    GridArray GrdArr = NC_ReadAregGridFile(eFileName);
    std::cerr << "We have GrdArr in AREG case\n";
    return GrdArr;
  }
  if (eModelName == "GEOS") {
    std::string HisPrefix = eTriple.HisPrefix;
    std::cerr << "Before FILE_DirectoryFilesSpecificExtension 2\n";
    std::string HisPrefixNaked = FILE_GetDirectoryOfFileName(HisPrefix);
    std::cerr << "HisPrefixaked=" << HisPrefixNaked << "\n";
    std::vector<std::string> ListFile =
        FILE_DirectoryFilesSpecificExtension(HisPrefixNaked, "nc4");
    std::cerr << " After FILE_DirectoryFilesSpecificExtension 2\n";
    if (ListFile.size() == 0) {
      std::cerr << "The list of files for GEOS retrieval of grid is empty\n";
      throw TerminalException{1};
    }
    std::string eFileName = ListFile[0];
    GridArray GrdArr = NC_ReadGeosGridFile(eFileName);
    std::cerr << "We have GrdArr in GEOS case\n";
    return GrdArr;
  }
  if (eModelName == "GRIB_DWD" || eModelName == "GRIB_GFS" ||
      eModelName == "GRIB_ECMWF" || eModelName == "GRIB_COSMO" ||
      eModelName == "GRIB_ALADIN" || eModelName == "GRIB_IFS") {
    std::string HisPrefix = eTriple.HisPrefix;
    std::vector<std::string> ListFile =
        FILE_DirectoryFilesSpecificExtension(HisPrefix, "grb");
    if (ListFile.size() == 0) {
      std::cerr << "The list of files is empty\n";
      throw TerminalException{1};
    }
    std::string eFileName = ListFile[0];
    GridArray GrdArr = GRIB_ReadGridArray(eFileName, eModelName);
    std::cerr << "We have GrdArr in GRIB case\n";
    return GrdArr;
  }
  if (eModelName == "GRIB_WAM_FORT30") {
    std::string eFileName = eTriple.HisPrefix;
    if (!IsExistingFile(eFileName)) {
      std::cerr << "The file eFileName = " << eFileName << " is missing\n";
      std::cerr << "This is set by HisPrefx and serves for the data storage\n";
      throw TerminalException{1};
    }
    return GRIB_ReadGridArray(eFileName, eModelName);
  }
  std::cerr << "Error in PRE_RETRIEVE_GRID_ARRAY\n";
  std::cerr << "Did not find the matching model for the grid\n";
  std::cerr << "Please correct\n";
  throw TerminalException{1};
}

void PrintGridArray(std::ostream &os, GridArray const &GrdArr) {
  os << "IsFE=" << GrdArr.IsFE << "\n";
  QuadArray eArr = GetQuadArray(GrdArr);
  os << "Lon (min/max)=" << eArr.MinLon << " / " << eArr.MaxLon << "\n";
  os << "Lat (min/max)=" << eArr.MinLat << " / " << eArr.MaxLat << "\n";
}

GridArray RETRIEVE_GRID_ARRAY(TripleModelDesc const &eTriple) {
  GridArray GrdArr = PRE_RETRIEVE_GRID_ARRAY(eTriple);
  std::string strSphericity = eTriple.RecGridSymb.Sphericity;
  if (strSphericity != "unset") {
    if (strSphericity != "Spherical" && strSphericity != "Cartesian") {
      std::cerr << "Error, we need the Sphericity option to be set to either\n";
      std::cerr << "unset: then the plotting software is guessing the right "
                   "value or reading directly from netcdf file\n";
      std::cerr << "Spherical: if the grid is spherical (needed only if the "
                   "guess is wrong)\n";
      std::cerr << "Cartesian: if the grid is cartesian (needed only if the "
                   "guess is wrong)\n";
      throw TerminalException{1};
    }
    if (strSphericity == "Spherical")
      GrdArr.IsSpherical = true;
    if (strSphericity == "Cartesian")
      GrdArr.IsSpherical = false;
  }
  std::cerr << "IsSpherical=" << GrdArr.IsSpherical << "\n";
  if (GrdArr.IsFE == 0)
    return GrdArr;
  if (eTriple.RecGridSymb.CutWorldMap) {
    CutWorldMap(GrdArr);
    std::cerr << "After CutWorldMap\n";
  }
  if (eTriple.RecGridSymb.HigherLatitudeCut) {
    double MinLatCut = eTriple.RecGridSymb.MinLatCut;
    double MaxLatCut = eTriple.RecGridSymb.MaxLatCut;
    CUT_HigherLatitude(GrdArr, MinLatCut, MaxLatCut);
    std::cerr << "After CUT_HigherLatitude\n";
  }
  CHECK_UnstructuredGrid(GrdArr);
  CHECK_CombinatorialGrid(GrdArr);
  CHECK_COORDINATE_ORIENTATION(GrdArr);
  std::cerr << "Before returning GrdArr in RETRIEVE_GRID_ARRAY\n";
  return GrdArr;
}

//
// NEMO code. After difficult mess, we decided to do the following structure.
// ---A single HisPrefix in the input
// ---An additional prefix "nut", "tem", etc. specifying the kind of files to
// get
//     We have that from the website of the marine.copernicus.eu
// ---For each files, a sequencing _0001.nc , _0002.nc , etc.
//    Each files has to be sequenced in the same way.
// ---A priori it is allowed to have files of differernt sizes. Just they have
// to be synchrone
//    over all variables.
// ---If the files are available as a
//

ArrayHistory NC_ReadArrayHistory_NEMO(TripleModelDesc const &eTriple) {
  std::string StringTime = "ocean_time";
  std::string PreModelName = eTriple.ModelName;
  std::string eModelName = GetKernelModelName(PreModelName);
  std::string HisPrefix = eTriple.HisPrefix;

  std::vector<std::string> ListFile =
      FILE_DirectoryMatchingPrefixExtension(HisPrefix, "nc");
  if (ListFile.size() == 0) {
    std::cerr << "We did not found any matching files in ListFile\n";
    throw TerminalException{1};
  }
  std::string eFile = ListFile[0];
  std::cerr << "eFile=" << eFile << "\n";
  std::vector<size_t> LPos;
  for (size_t i = 0; i < eFile.size(); i++) {
    if (eFile.substr(i, 1) == "_")
      LPos.push_back(i);
  }
  size_t len = LPos.size();
  if (len < 2) {
    std::cerr << "We need at least two _ in the filename\n";
    std::cerr << "For example a model is Reset_file/med_nut_0001.nc\n";
    std::cerr << "But the file we have is named eFile=" << eFile << "\n";
    throw TerminalException{1};
  }
  std::string HisPrefixRed = eFile.substr(0, LPos[len - 1] + 1);
  std::string HisPrefixNew = eFile.substr(0, LPos[len - 2] + 1);
  std::cerr << "HisPrefix=" << HisPrefix << "\n";
  std::cerr << "HisPrefixRed=" << HisPrefixRed << "\n";
  std::cerr << "HisPrefixNew=" << HisPrefixNew << "\n";
  ArrayHistory eArr = Sequential_ReadArrayHistory(HisPrefixRed, "time");
  eArr.nbDigit = 4;
  eArr.HisPrefix = HisPrefixNew;
  std::cerr << "NEMO : |ListFile|=" << ListFile.size() << "\n";
  for (auto &eFile : ListFile) {
    std::vector<std::string> LStr = STRING_Split(eFile, "_");
    std::cerr << "eFile=" << eFile << "\n";
    if (LStr.size() >= 3) { // Corresponding for example to data/med_nut_0001.nc
      std::string postfix = LStr[LStr.size() - 2];
      NEMO_vars nemo_vars = ReadNEMO_vars(eFile);
      std::cerr << "NEMO |List2D_vars|=" << nemo_vars.List2D_vars.size()
                << " |List3D_vars|=" << nemo_vars.List3D_vars.size() << "\n";
      for (auto &eVar : nemo_vars.List2D_vars) {
        std::cerr << "2D : eFile=" << eFile << " postfix=" << postfix
                  << " eVar=" << eVar << "\n";
        eArr.NEMO_vars_to_postfix[eVar] = postfix;
      }
      for (auto &eVar : nemo_vars.List3D_vars) {
        std::cerr << "3D : eFile=" << eFile << " postfix=" << postfix
                  << " eVar=" << eVar << "\n";
        eArr.NEMO_vars_to_postfix[eVar] = postfix;
      }
    }
  }
  std::cerr << "Returning eArr\n";
  return eArr;
}

ArrayHistory NC_ReadArrayHistory(TripleModelDesc const &eTriple) {
  std::string StringTime = "ocean_time";
  std::string PreModelName = eTriple.ModelName;
  std::string eModelName = GetKernelModelName(PreModelName);
  std::string HisPrefix = eTriple.HisPrefix;
  std::cerr << "Debug: NC_ReadArrayHistory\n";
  std::cerr << "NC_ReadArrayHistory : eModelName=" << eModelName << "\n";
  // special models first
  if (eModelName == "WW3") {
    std::string HisFile = GET_GRID_FILE(eTriple);
    return WW3_ReadArrayHistory(HisFile, HisPrefix);
  }
  if (eModelName == "ROMS_IVICA" || eModelName == "WWM_DAILY")
    return Sequential_ReadArrayHistory(HisPrefix, "ocean_time");
  if (eModelName == "SCHISM_SFLUX")
    return NC_ReadArrayHistory_Kernel(HisPrefix, "time", 3);
  if (eModelName == "NEMO")
    return NC_ReadArrayHistory_NEMO(eTriple);
  if (eModelName == "AREG")
    return Sequential_ReadArrayHistory(HisPrefix, "time");
  if (eModelName == "GEOS")
    return Sequential_ReadArrayHistory(HisPrefix, "time");
  if (eModelName == "HYCOM")
    return Sequential_ReadArrayHistory(HisPrefix, "time");
  return NC_ReadArrayHistory_Kernel(HisPrefix, StringTime, 4);
}

/*
   When we have several runs (typical in operational runs)
   then we need to select the best one.
 */
double GetOptimalTimeShiftLength(std::string const &eModelName) {
  std::vector<std::string> LStr = STRING_Split(eModelName, ":");
  for (int i = 1; i < int(LStr.size()); i++) {
    std::vector<std::string> LStrB = STRING_Split(LStr[i], "_");
    if (LStrB[0] == "optimaltime") {
      double eValHour;
      std::istringstream(LStrB[1]) >> eValHour;
      double eValDay = eValHour / double(24);
      return eValDay;
    }
  }
  return 0;
}

bool RetrieveAllStates(std::string const &eModelName) {
  std::vector<std::string> LStr = STRING_Split(eModelName, ":");
  for (int i = 1; i < int(LStr.size()); i++)
    if (LStr[i] == "retrieveallstates")
      return true;
  return false;
}

ArrayHistory
GRIB_ReadArrayHistory_Kernel(std::vector<std::string> const &ListFile,
                             std::string const &eModelName) {
  int nbFile = ListFile.size();
  std::cerr << "Beginning of GRIB_ReadArrayHistory_Kernel nbFile=" << nbFile
            << "\n";
  //
  // Determining parameters of the search
  //
  GRIB_CheckAllowedKeywords(eModelName);
  double OptimalShiftDay = GetOptimalTimeShiftLength(eModelName);
  std::cerr << "OptimalShiftDay = " << OptimalShiftDay << "\n";
  double MaxErrorTime = 0.01;
  bool RetAllStates = RetrieveAllStates(eModelName);
  //
  // Determining the list of messages
  //
  std::vector<GRIB_MessageInfo> ListAllMessage;
  std::set<std::string> TotalListShortName;
  bool PrintGribFileInfo = false;
  for (int iFile = 0; iFile < nbFile; iFile++) {
    std::string eFile = ListFile[iFile];
    std::cerr << "iFile=" << iFile << " / " << nbFile << " eFile=" << eFile
              << "\n";
    std::vector<GRIB_MessageInfo> ListMessage =
        GRIB_GetAllMessagesFromFile(eFile, eModelName);
    int nbMessage = ListMessage.size();
    if (PrintGribFileInfo) {
      std::cerr << "eFile = " << eFile << " nbMessage=" << nbMessage << "\n";
      PrintVectorGRIBmessageInfo(std::cerr, ListMessage);
      for (auto &eMesg : ListMessage)
        std::cerr << "  shortName=" << eMesg.shortName << "\n";
    }
    if (nbMessage == 0) {
      std::cerr << "Remark: eFile = " << eFile << " has zero messages\n";
    }
    for (auto &eMesg : ListMessage) {
      ListAllMessage.push_back(eMesg);
      TotalListShortName.insert(eMesg.shortName);
    }
  }
  int TotalNbMessage = ListAllMessage.size();
  std::cerr << "TotalNbMessage=" << TotalNbMessage << "\n";
  if (TotalNbMessage == 0) {
    std::cerr << "TotalNbMessage=" << TotalNbMessage << "\n";
    std::cerr << "|ListFile|=" << ListFile.size() << "\n";
    std::cerr << "We have zero messages. No work can be done\n";
    throw TerminalException{1};
  }
  //  std::cerr << "1: |ListAllMessage|=" << ListAllMessage.size() << "\n";
  //
  // Now reordering the messages
  //
  sort(ListAllMessage.begin(), ListAllMessage.end(),
       [&](GRIB_MessageInfo const &a, GRIB_MessageInfo const &b) -> bool {
         if (a.time < b.time)
           return true;
         return false;
       });
  std::cerr << "After sorting by time\n";
  //
  // Determining the list of times
  //
  double TimePrev = ListAllMessage[0].time;
  std::vector<double> ListTime{TimePrev};
  std::vector<int> ListITime(TotalNbMessage);
  std::vector<int> ListIFile(TotalNbMessage, -1);
  std::vector<int> ListIRec(TotalNbMessage, -1);
  int posTime = 0;
  for (int iMesg = 0; iMesg < TotalNbMessage; iMesg++) {
    GRIB_MessageInfo eMesg = ListAllMessage[iMesg];
    double eTime = eMesg.time;
    double TimeDiff = fabs(eTime - TimePrev);
    if (TimeDiff > MaxErrorTime) {
      ListTime.push_back(eTime);
      TimePrev = eTime;
      posTime++;
    }
    ListITime[iMesg] = posTime;
  }
  //
  // Determination of List of starttime from the list of messages.
  //
  std::vector<double> ListStartTime;
  std::vector<int> ListIStartTime(TotalNbMessage);
  double tolDay = double(1) / double(10000);
  auto GetIStartTime = [&](double const &timeStart) -> int {
    int nbTimeStart = ListStartTime.size();
    for (int iTimeStart = 0; iTimeStart < nbTimeStart; iTimeStart++)
      if (fabs(timeStart - ListStartTime[iTimeStart]) < tolDay)
        return iTimeStart;
    ListStartTime.push_back(timeStart);
    return nbTimeStart;
  };
  for (int iMesg = 0; iMesg < TotalNbMessage; iMesg++) {
    double eStartTime = ListAllMessage[iMesg].timeStart;
    int iTimeStart = GetIStartTime(eStartTime);
    ListIStartTime[iMesg] = iTimeStart;
  }
  std::cerr << "ListIStartTime determined\n";
  int nbTimeStart = ListStartTime.size();
  std::vector<double> ListEndTime(nbTimeStart, 0);
  std::vector<std::vector<int>> ListListIMesg(nbTimeStart);
  for (int iMesg = 0; iMesg < TotalNbMessage; iMesg++) {
    int iTimeStart = ListIStartTime[iMesg];
    ListListIMesg[iTimeStart].push_back(iMesg);
    double time = ListAllMessage[iMesg].time;
    if (time > ListEndTime[iTimeStart])
      ListEndTime[iTimeStart] = time;
  }
  std::cerr << "ListEndTime determined\n";
  bool PrintEssentialInfo = false;
  if (PrintEssentialInfo) {
    for (int iMesg = 0; iMesg < TotalNbMessage; iMesg++) {
      GRIB_MessageInfo eMesg = ListAllMessage[iMesg];
      double time = eMesg.time;
      double timestart = eMesg.timeStart;
      double deltaTimeTimeStart = timestart - time;
      std::cerr << "iMesg=" << iMesg << " shortName=" << eMesg.shortName
                << " time=" << eMesg.time << " deltaTT=" << deltaTimeTimeStart
                << " file=" << eMesg.FileName << " p=" << ListITime[iMesg]
                << "\n";
    }
  }
  int nbTime = ListTime.size();
  std::cerr << "nbTime=" << nbTime << "\n";
  //
  // Determining the list of messages according to time
  //
  std::vector<std::vector<GRIB_MessageInfo>> ListListMessages(nbTime);
  for (int iMesg = 0; iMesg < TotalNbMessage; iMesg++) {
    GRIB_MessageInfo eMesg = ListAllMessage[iMesg];
    //    std::cerr << "iMesg=" << iMesg << " eMesg.shortName=" <<
    //    eMesg.shortName << "\n";
    int iTime = ListITime[iMesg];
    //    std::cerr << "iTime=" << iTime << "\n";
    if (iTime < 0 || iTime >= nbTime) {
      std::cerr << "iTime=" << iTime << " but nbTime=" << nbTime << "\n";
      throw TerminalException{1};
    }
    ListListMessages[iTime].push_back(eMesg);
  }
  std::cerr << "After determination of ListListMessages\n";
  //
  // Now cleaning the entries by prefering the entries with timeStart +
  // OptimalTimeDay as small as possible.
  //
  std::map<std::string,
           std::vector<std::pair<double, std::vector<GRIB_MessageInfo>>>>
      FullOrganizedInfo;
  if (RetAllStates) {
    for (auto &eShortName : TotalListShortName) {
      FullOrganizedInfo[eShortName] = {};
    }
  }
  std::set<std::string> SetRawNames;
  bool PrintInfo = false;
  for (int iTime = 0; iTime < nbTime; iTime++) {
    double eTime = ListTime[iTime];
    std::string strPres = DATE_ConvertMjd2mystringPres(eTime);
    std::vector<GRIB_MessageInfo> ListMessages = ListListMessages[iTime];
    if (PrintInfo) {
      std::cerr << "iTime=" << iTime << " / " << nbTime << " eTime=" << eTime
                << " date=" << strPres << "\n";
      std::cerr << "1: |ListMessages|=" << ListMessages.size() << "\n";
    }
    std::set<std::string> ListShortName;
    for (auto &eMesg : ListMessages) {
      ListShortName.insert(eMesg.shortName);
      SetRawNames.insert(eMesg.shortName);
    }
    std::vector<GRIB_MessageInfo> NewListMessages;
    std::vector<GRIB_MessageInfo> ListMatchingMessages;
    for (auto &eShortName : ListShortName) {
      GRIB_MessageInfo NewMesg;
      double minPenaltyFct = 10 ^ (30);
      int nbMatch = 0;
      for (auto &eMesg : ListMessages) {
        double ePenaltyFct =
            fabs(eMesg.timeStart + OptimalShiftDay - eMesg.time);
        if (eMesg.shortName == eShortName) {
          nbMatch++;
          if (RetAllStates) {
            ListMatchingMessages.push_back(eMesg);
          }
          if (PrintInfo)
            std::cerr << " " << eMesg.FileName << "\n";
          if (ePenaltyFct < minPenaltyFct) {
            NewMesg = eMesg;
            minPenaltyFct = ePenaltyFct;
          }
        }
      }
      if (PrintInfo) {
        std::cerr << "1: eShortName=" << eShortName << " nbMatch=" << nbMatch
                  << " Selected=" << NewMesg.FileName << "\n";
      }
      NewListMessages.push_back(NewMesg);
      if (RetAllStates) {
        std::pair<double, std::vector<GRIB_MessageInfo>> ePair{
            eTime, ListMatchingMessages};
        FullOrganizedInfo[eShortName].push_back(ePair);
      }
    }
    ListListMessages[iTime] = NewListMessages;
  }
  std::cerr << "After determination of ListListMessages (second operation)\n";
  std::vector<std::string> RawVarNames;
  std::unordered_map<std::string, std::vector<int>> MatchingByVariable;
  for (auto &eName : SetRawNames) {
    RawVarNames.push_back(eName);
    MatchingByVariable[eName] = {};
  }
  for (int iTime = 0; iTime < nbTime; iTime++) {
    for (auto &eMesg : ListListMessages[iTime]) {
      std::string eName = eMesg.shortName;
      MatchingByVariable[eName].push_back(iTime);
    }
  }
  std::cerr << "After to MatchingByVariable\n";
  double FirstTime = ListTime[0];
  double LastTime = ListTime[nbTime - 1];
  ArrayHistory eArr;
  eArr.nbFile = nbFile;
  eArr.nbTime = nbTime;
  eArr.ListListMessages = ListListMessages;
  eArr.ListAllMessage = ListAllMessage;
  eArr.ListStartTime = ListStartTime;
  eArr.ListEndTime = ListEndTime;
  eArr.ListListIMesg = ListListIMesg;
  eArr.ListIStartTime = ListIStartTime;
  eArr.RawVarNames = RawVarNames;
  eArr.MatchingByVariable = MatchingByVariable;
  eArr.FullOrganizedInfo = FullOrganizedInfo;
  eArr.ListITime = ListITime;
  eArr.ListIFile = ListIFile;
  eArr.ListIRec = ListIRec;
  eArr.ListTime = ListTime;
  eArr.FirstTime = FirstTime;
  eArr.LastTime = LastTime;
  eArr.KindArchive = "GRIB";
  eArr.TimeSteppingInfo = "classic";
  eArr.AppendVarName = false;
  std::cerr << "Ending of GRIB_ReadArrayHistory\n";
  return eArr;
}

ArrayHistory GRIB_ReadArrayHistory(std::string const &HisPrefix,
                                   std::string const &eModelName) {
  //
  // Determining the list of files.
  //
  std::cerr << "Begin of GRIB_ReadArrayHistory\n";
  std::vector<std::string> ListFile;
  if (IsExistingFile(HisPrefix) && FILE_IsRegularFile(HisPrefix)) {
    ListFile = {HisPrefix};
  } else {
    ListFile = FILE_DirectoryFilesSpecificExtension(HisPrefix, "grb");
  }
  return GRIB_ReadArrayHistory_Kernel(ListFile, eModelName);
}

ArrayHistory ReadArrayHistory(TripleModelDesc const &eTriple) {
  ArrayHistory eArr;
  std::string PreModelName = eTriple.ModelName;
  std::string eModelName = GetKernelModelName(PreModelName);
  CHECK_Model_Allowedness(eModelName);
  std::vector<std::string> ListModelGrib{
      "GRIB_DWD",    "GRIB_GFS", "GRIB_COSMO",     "GRIB_ECMWF",
      "GRIB_ALADIN", "GRIB_IFS", "GRIB_WAM_FORT30"};
  if (PositionVect(ListModelGrib, eModelName) != -1) {
    ArrayHistory eArr = GRIB_ReadArrayHistory(eTriple.HisPrefix, PreModelName);
    eArr.eModelName = eModelName;
    return eArr;
  }
  std::vector<std::string> ListModelNetcdf{
      "COSMO",       "WAM",          "ROMS",
      "ROMS_IVICA",  "WWM",          "WWM_DAILY",
      "WW3",         "SCHISM_SFLUX", "SCHISM_NETCDF_OUT",
      "RECTANGULAR", "WRF",          "UNRUNOFF",
      "IVICA_UVP",   "NEMO",         "HYCOM",
      "AREG",        "GEOS"};
  if (PositionVect(ListModelNetcdf, eModelName) != -1) {
    ArrayHistory eArr = NC_ReadArrayHistory(eTriple);
    eArr.eModelName = eModelName;
    return eArr;
  }
  std::cerr << "Error in ReadArrayHistory. Could not find a matching method "
               "for creating the array history\n";
  throw TerminalException{1};
}

VerticalInfo GetVerticalInfo(int const &N) {
  MyVector<double> Hz(N);
  MyVector<double> z_w(N + 1);
  MyVector<double> z_r(N);
  return {std::move(Hz), std::move(z_w), std::move(z_r)};
}

void ComputeHz(ARVDtyp const &ARVD, double const &hwater, double const &eZeta,
               VerticalInfo &eVert) {
  int Vtrans = ARVD.Vtransform;
  int N = ARVD.N;
  if (Vtrans == 1) {
    double hinv = 1 / hwater;
    eVert.z_w(0) = -hwater;
    for (int k = 1; k <= N; k++) {
      double cff_r = ARVD.hc * (ARVD.sc_r(k - 1) - ARVD.Cs_r(k - 1));
      double cff_w = ARVD.hc * (ARVD.sc_w(k) - ARVD.Cs_w(k));
      double cff1_r = ARVD.Cs_r(k - 1);
      double cff1_w = ARVD.Cs_w(k);
      double z_w0 = cff_w + cff1_w * hwater;
      eVert.z_w(k) = z_w0 + eZeta * (1 + z_w0 * hinv);
      double z_r0 = cff_r + cff1_r * hwater;
      eVert.z_r(k - 1) = z_r0 + eZeta * (1 + z_r0 * hinv);
      eVert.Hz(k - 1) = eVert.z_w(k) - eVert.z_w(k - 1);
    }
    return;
  }
  if (Vtrans == 2) {
    eVert.z_w(0) = -hwater;
    double hinv = 1 / (ARVD.hc + hwater);
    for (int k = 1; k <= N; k++) {
      double cff_r = ARVD.hc * ARVD.sc_r(k - 1);
      double cff_w = ARVD.hc * ARVD.sc_w(k);
      double cff1_r = ARVD.Cs_r(k - 1);
      double cff1_w = ARVD.Cs_w(k);
      double cff2_r = (cff_r + cff1_r * hwater) * hinv;
      double cff2_w = (cff_w + cff1_w * hwater) * hinv;
      eVert.z_w(k) = eZeta + (eZeta + hwater) * cff2_w;
      eVert.z_r(k - 1) = eZeta + (eZeta + hwater) * cff2_r;
      eVert.Hz(k - 1) = eVert.z_w(k) - eVert.z_w(k - 1);
    }
    return;
  }
  std::cerr << "Failed to find matching entry of Vtrans Vtrans=" << Vtrans
            << "\n";
  throw TerminalException{1};
}

Eigen::Tensor<double, 3>
ROMS_ComputeVerticalGlobalCoordinate_r(GridArray const &GrdArr,
                                       MyMatrix<double> const &zeta) {
  int N = GrdArr.ARVD.N;
  int eta_rho = zeta.rows();
  int xi_rho = zeta.cols();
  VerticalInfo eVert = GetVerticalInfo(N);
  Eigen::Tensor<double, 3> Zmat(N, eta_rho, xi_rho);
  const MyMatrix<double> &DEP = GetDEP(GrdArr.GrdArrRho);
  for (int i = 0; i < eta_rho; i++)
    for (int j = 0; j < xi_rho; j++) {
      for (int k = 0; k < N; k++)
        eVert.z_r(k) = 0;
      if (GrdArr.GrdArrRho.MSK(i, j) == 1)
        ComputeHz(GrdArr.ARVD, DEP(i, j), zeta(i, j), eVert);
      for (int k = 0; k < N; k++)
        Zmat(k, i, j) = eVert.z_r(k);
    }
  return Zmat;
}

Eigen::Tensor<double, 3>
ROMS_ComputeVerticalGlobalCoordinate_w(GridArray const &GrdArr,
                                       MyMatrix<double> const &zeta) {
  int N = GrdArr.ARVD.N;
  int eta_rho = zeta.rows();
  int xi_rho = zeta.cols();
  VerticalInfo eVert = GetVerticalInfo(N);
  Eigen::Tensor<double, 3> Zmat(N + 1, eta_rho, xi_rho);
  const MyMatrix<double> &DEP = GetDEP(GrdArr.GrdArrRho);
  for (int i = 0; i < eta_rho; i++)
    for (int j = 0; j < xi_rho; j++) {
      for (int k = 0; k < N + 1; k++)
        eVert.z_w(k) = 0;
      if (GrdArr.GrdArrRho.MSK(i, j) == 1)
        ComputeHz(GrdArr.ARVD, DEP(i, j), zeta(i, j), eVert);
      for (int k = 0; k < N + 1; k++)
        Zmat(k, i, j) = eVert.z_w(k);
    }
  return Zmat;
}

struct PairMSKfield {
  MyMatrix<int> MSK;
  MyMatrix<double> field;
};

PairMSKfield
VerticalInterpolation_P1_W(ARVDtyp const &ARVD, MyMatrix<double> const &h,
                           MyMatrix<double> const &zeta,
                           MyMatrix<int> const &MSK, double const &dep,
                           Eigen::Tensor<double, 3> const &VertField_W) {
  int eta = h.rows();
  int xi = h.cols();
  MyMatrix<int> MSKret(eta, xi);
  MyMatrix<double> FieldRet(eta, xi);
  if (!IsEqualSizeMatrices(h, zeta)) {
    std::cerr << "   |h|=" << h.rows() << " / " << h.cols() << "\n";
    std::cerr << "|zeta|=" << zeta.rows() << " / " << zeta.cols() << "\n";
    std::cerr << "Most likely the grid array does not match the history file\n";
    throw TerminalException{1};
  }
  int N = ARVD.N;
  VerticalInfo eVert = GetVerticalInfo(N);
  for (int i = 0; i < eta; i++)
    for (int j = 0; j < xi; j++) {
      int eMSK = MSK(i, j);
      double eField = 0;
      if (eMSK == 1) {
        if (dep < -h(i, j)) {
          eMSK = 0;
        } else {
          ComputeHz(ARVD, h(i, j), zeta(i, j), eVert);
          for (int iVert = 0; iVert < N; iVert++) {
            double dep1 = eVert.z_w(iVert);
            double dep2 = eVert.z_w(iVert + 1);
            if (dep1 >= dep && dep <= dep2) {
              double alpha1 = (dep2 - dep) / (dep2 - dep1);
              double alpha2 = (dep - dep1) / (dep2 - dep1);
              eField = alpha1 * VertField_W(iVert, i, j) +
                       alpha2 * VertField_W(iVert + 1, i, j);
            }
          }
        }
      }
      MSKret(i, j) = eMSK;
      FieldRet(i, j) = eField;
    }
  return {std::move(MSKret), std::move(FieldRet)};
}

MyMatrix<double>
VerticalInterpolation_P2_W(ARVDtyp const &ARVD, MyMatrix<double> const &h,
                           MyMatrix<double> const &zeta,
                           MyMatrix<int> const &MSK, double const &dep,
                           Eigen::Tensor<double, 3> const &VertField_W) {
  if (!IsEqualSizeMatrices(h, zeta)) {
    std::cerr << "   |h|=" << h.rows() << " / " << h.cols() << "\n";
    std::cerr << "|zeta|=" << zeta.rows() << " / " << zeta.cols() << "\n";
    std::cerr << "Most likely the grid array does not match the history file\n";
    throw TerminalException{1};
  }
  PairMSKfield ePair =
      VerticalInterpolation_P1_W(ARVD, h, zeta, MSK, dep, VertField_W);
  int eta = h.rows();
  int xi = h.cols();
  MyMatrix<double> FieldRet(eta, xi);
  for (int i = 0; i < eta; i++)
    for (int j = 0; j < xi; j++) {
      double eVal;
      if (ePair.MSK(i, j) == 0) {
        eVal = 0;
      } else {
        eVal = ePair.field(i, j);
      }
      FieldRet(i, j) = eVal;
    }
  return FieldRet;
}

MyVector<double> GetVertCoord_R(ARVDtyp const &ARVD, double const &eDep,
                                double const &eZeta) {
  if (ARVD.Zcoordinate)
    return ARVD.ListZ_r;
  if (ARVD.ModelName == "ROMS") {
    VerticalInfo eVert = GetVerticalInfo(ARVD.N);
    ComputeHz(ARVD, eDep, eZeta, eVert);
    return eVert.z_r;
  }
  std::cerr << "Did not find any matching entry\n";
  throw TerminalException{1};
}

MyVector<double> GetVertCoord_W(ARVDtyp const &ARVD, double const &eDep,
                                double const &eZeta) {
  if (ARVD.Zcoordinate)
    return ARVD.ListZ_w;
  if (ARVD.ModelName == "ROMS") {
    VerticalInfo eVert = GetVerticalInfo(ARVD.N);
    ComputeHz(ARVD, eDep, eZeta, eVert);
    return eVert.z_w;
  }
  std::cerr << "Did not find any matching entry\n";
  throw TerminalException{1};
}

Eigen::Tensor<double, 3>
VerticalInterpolationTensor_R(GridArray const &GrdArrOut, ARVDtyp const &ARVDin,
                              MyMatrix<double> const &DEPin,
                              Eigen::Tensor<double, 3> const &TensIn) {
  int eta_rho = GrdArrOut.GrdArrRho.LON.rows();
  int xi_rho = GrdArrOut.GrdArrRho.LON.cols();
  int NvertOut = GrdArrOut.ARVD.N;
  int NvertIn = ARVDin.N;
  Eigen::Tensor<double, 3> TensOut(NvertOut, eta_rho, xi_rho);
  const MyMatrix<double> &DEP = GetDEP(GrdArrOut.GrdArrRho);
  for (int i = 0; i < eta_rho; i++)
    for (int j = 0; j < xi_rho; j++) {
      double eDepOut = DEP(i, j);
      double eDepIn = DEPin(i, j);
      double eZeta = 0;
      MyVector<double> Zr_out = GetVertCoord_R(GrdArrOut.ARVD, eDepOut, eZeta);
      MyVector<double> Zr_in = GetVertCoord_R(ARVDin, eDepIn, eZeta);
      //      std::cerr << "Zr_out=" << Zr_out << "\n";
      //      std::cerr << "Zr_in=" << Zr_in << "\n";
      //      std::cerr << "NvertIn=" << NvertIn << "\n";
      //      std::cerr << "|Zr_out|=" << Zr_out.size() << " |Zr_in|=" <<
      //      Zr_in.size() << " ARVDin.N=" << ARVDin.N << "\n"; std::cerr <<
      //      "NvertOut=" << NvertOut << " NvertIn=" << NvertIn << "\n";
      for (int k = 0; k < NvertOut; k++) {
        //	std::cerr << "k=" << k << "\n";
        double depW = Zr_out(k);
        //	std::cerr << "depW=" << depW << "\n";
        double eValOut = 0;
        //        std::cerr << "|TensIn|=" << TensIn.size() << " i=" << i << "
        //        j=" << j << " k=" << k << " depW=" << depW << "\n";
        bool IsAssigned = false;
        if (depW < Zr_in(0)) {
          //          std::cerr << "Case 1\n";
          eValOut = TensIn(0, i, j);
          IsAssigned = true;
        } else {
          if (depW > Zr_in(NvertIn - 1)) {
            //            std::cerr << "Case 2\n";
            eValOut = TensIn(NvertIn - 1, i, j);
            IsAssigned = true;
          } else {
            //            std::cerr << "Case 3\n";
            for (int u = 1; u < NvertIn; u++) {
              double dep1 = Zr_in(u - 1);
              double dep2 = Zr_in(u);
              double alpha1 = (dep2 - depW) / (dep2 - dep1);
              double alpha2 = (depW - dep1) / (dep2 - dep1);
              if (dep1 <= depW && depW <= dep2) {
                eValOut =
                    alpha1 * TensIn(u - 1, i, j) + alpha2 * TensIn(u, i, j);
                IsAssigned = true;
              }
            }
          }
        }
        if (!IsAssigned) {
          std::cerr << "Failed to do the assignation of the variable\n";
          throw TerminalException{1};
        }
        TensOut(k, i, j) = eValOut;
      }
    }
  return TensOut;
}

MyMatrix<double> ConvertBaroclinic_to_Barotropic_ARVD_Coord(
    Eigen::Tensor<double, 3> const &F3, MyMatrix<double> const &zeta,
    ARVDtyp const &ARVD, CoordGridArrayFD const &GrdArrRho) {
  auto LDim = F3.dimensions();
  int N = LDim[0];
  int eta_rho = LDim[1];
  int xi_rho = LDim[2];
  int eta_rho_grid = GrdArrRho.LON.rows();
  int xi_rho_grid = GrdArrRho.LON.cols();
  if (N != ARVD.N) {
    std::cerr << "First dimension of F1 is N=" << N << "\n";
    std::cerr << "But GrdArr.ARVD.N=" << ARVD.N << "\n";
    throw TerminalException{1};
  }
  if (eta_rho != eta_rho_grid || xi_rho != xi_rho_grid) {
    std::cerr << "eta_rho      = " << eta_rho << " xi_rho      = " << xi_rho
              << "\n";
    std::cerr << "eta_rho_grid = " << eta_rho_grid << " xi_rho_grid = "
              << "\n";
    std::cerr << "They should be equal\n";
    throw TerminalException{1};
  }
  MyMatrix<double> F(eta_rho, xi_rho);
  VerticalInfo eVert = GetVerticalInfo(N);
  const MyMatrix<double> &DEP = GetDEP(GrdArrRho);
  for (int i = 0; i < eta_rho; i++)
    for (int j = 0; j < xi_rho; j++) {
      double eVal = 0;
      if (GrdArrRho.MSK(i, j) == 1) {
        ComputeHz(ARVD, DEP(i, j), zeta(i, j), eVert);
        double VertInt = 0;
        for (int k = 0; k < N; k++) {
          double dep = eVert.z_w(k + 1) - eVert.z_w(k);
          VertInt += dep * F3(k, i, j);
        }
        double TotalDep = zeta(i, j) + DEP(i, j);
        eVal = VertInt / TotalDep;
      }
      F(i, j) = eVal;
    }
  return F;
}

MyMatrix<double>
ConvertBaroclinic_to_Barotropic(Eigen::Tensor<double, 3> const &F3,
                                MyMatrix<double> const &zeta,
                                GridArray const &GrdArr) {
  return ConvertBaroclinic_to_Barotropic_ARVD_Coord(F3, zeta, GrdArr.ARVD,
                                                    GrdArr.GrdArrRho);
}

PairMSKfield
VerticalInterpolation_P1_R(ARVDtyp const &ARVD, MyMatrix<double> const &h,
                           MyMatrix<double> const &zeta,
                           MyMatrix<uint8_t> const &MSK, double const &dep,
                           Eigen::Tensor<double, 3> const &VertField_R,
                           int const &Choice) {
  if (!IsEqualSizeMatrices(h, zeta)) {
    std::cerr << "   |h|=" << h.rows() << " / " << h.cols() << "\n";
    std::cerr << "|zeta|=" << zeta.rows() << " / " << zeta.cols() << "\n";
    std::cerr << "Most likely the grid array does not match the history file\n";
    throw TerminalException{1};
  }
  if (dep > 0) {
    std::cerr << "dep=" << dep << "\n";
    std::cerr << "dep should be negative because we are looking after oceanic "
                 "interpolation\n";
    throw TerminalException{1};
  }
  int eta = h.rows();
  int xi = h.cols();
  MyMatrix<int> MSKret(eta, xi);
  MyMatrix<double> FieldRet(eta, xi);
  int N = ARVD.N;
  VerticalInfo eVert = GetVerticalInfo(N);
  std::cerr << "Choice=" << Choice << "\n";
  for (int i = 0; i < eta; i++)
    for (int j = 0; j < xi; j++) {
      int eMSK = MSK(i, j);
      double depW;
      if (Choice == 1)
        depW = dep;
      else
        depW = dep + zeta(i, j);
      double eField = 0;
      if (eMSK == 1) {
        if (depW < -h(i, j)) {
          eMSK = 0;
        } else {
          ComputeHz(ARVD, h(i, j), zeta(i, j), eVert);
          bool WeMatch = false;
          if (depW <= eVert.z_r(0)) {
            eField = VertField_R(0, i, j);
            WeMatch = true;
          }
          if (eVert.z_r(N - 1) <= depW) {
            eField = VertField_R(N - 1, i, j);
            WeMatch = true;
          } else {
            for (int iVert = 0; iVert < N - 1; iVert++) {
              double dep1 = eVert.z_r(iVert);
              double dep2 = eVert.z_r(iVert + 1);
              if (dep1 <= depW && depW <= dep2) {
                double alpha1 = (dep2 - depW) / (dep2 - dep1);
                double alpha2 = (depW - dep1) / (dep2 - dep1);
                eField = alpha1 * VertField_R(iVert, i, j) +
                         alpha2 * VertField_R(iVert + 1, i, j);
                WeMatch = true;
              }
            }
          }
          if (!WeMatch) {
            std::cerr << "No assignation for the vertical interpolation\n";
            std::cerr << "i=" << i << " j=" << j << "\n";
            std::cerr << "zeta=" << zeta(i, j) << " h=" << h(i, j) << "\n";
            std::cerr << "dep=" << dep << " N=" << N << "\n";
            std::cerr << "depW=" << depW << "\n";
            for (int iVert = 0; iVert < N; iVert++)
              std::cerr << "  iVert=" << iVert << " z_z=" << eVert.z_r(iVert)
                        << "\n";
            throw TerminalException{1};
          }
        }
      }
      MSKret(i, j) = eMSK;
      FieldRet(i, j) = eField;
    }
  return {std::move(MSKret), std::move(FieldRet)};
}

MyMatrix<double> VerticalInterpolation_SCHISM_ZNL(
    Eigen::Tensor<double, 3> const &znl, MyMatrix<double> const &zeta,
    double const &dep, Eigen::Tensor<double, 3> const &VertField_R,
    int const &Choice) {
  if (Choice != 1 && Choice != 2) {
    std::cerr << "We should have Choice=1 or 2\n";
    throw TerminalException{1};
  }
  int eta = zeta.rows();
  int xi = zeta.cols();
  auto LDim = znl.dimensions();
  int N = LDim[0];
  MyMatrix<double> FieldRet(eta, xi);
  for (int i = 0; i < eta; i++)
    for (int j = 0; j < xi; j++) {
      double depSearch;
      if (Choice == 1)
        depSearch = dep;
      else
        depSearch = dep + zeta(i, j);
      double eField = 0;
      if (depSearch <= znl(0, i, j)) {
        eField = VertField_R(0, i, j);
      }
      if (znl(N - 1, i, j) <= depSearch) {
        eField = VertField_R(N - 1, i, j);
      } else {
        for (int iVert = 0; iVert < N - 1; iVert++) {
          double dep1 = znl(iVert, i, j);
          double dep2 = znl(iVert + 1, i, j);
          if (dep1 <= depSearch && depSearch <= dep2) {
            double alpha1 = (dep2 - depSearch) / (dep2 - dep1);
            double alpha2 = (depSearch - dep1) / (dep2 - dep1);
            eField = alpha1 * VertField_R(iVert, i, j) +
                     alpha2 * VertField_R(iVert + 1, i, j);
          }
        }
      }
      FieldRet(i, j) = eField;
    }
  return FieldRet;
}

MyMatrix<double>
VerticalInterpolation_P2_R(ARVDtyp const &ARVD, MyMatrix<double> const &h,
                           MyMatrix<double> const &zeta,
                           MyMatrix<uint8_t> const &MSK, double const &dep,
                           Eigen::Tensor<double, 3> const &VertField_R,
                           int const &Choice) {
  std::cerr << "We are in VerticalInterpolation_P2_R\n";
  if (!IsEqualSizeMatrices(h, zeta)) {
    std::cerr << "   |h|=" << h.rows() << " / " << h.cols() << "\n";
    std::cerr << "|zeta|=" << zeta.rows() << " / " << zeta.cols() << "\n";
    std::cerr << "Most likely the grid array does not match the history file\n";
    throw TerminalException{1};
  }
  PairMSKfield ePair =
      VerticalInterpolation_P1_R(ARVD, h, zeta, MSK, dep, VertField_R, Choice);
  int eta = h.rows();
  int xi = h.cols();
  MyMatrix<double> FieldRet(eta, xi);
  for (int i = 0; i < eta; i++)
    for (int j = 0; j < xi; j++) {
      double eVal;
      if (ePair.MSK(i, j) == 0) {
        eVal = 0;
      } else {
        eVal = ePair.field(i, j);
      }
      FieldRet(i, j) = eVal;
    }
  return FieldRet;
}

double GetUnitInMeter(double const &eVal, std::string const &unit) {
  if (unit == "m")
    return eVal;
  if (unit == "dm")
    return eVal * double(0.1);
  if (unit == "cm")
    return eVal * double(0.01);
  if (unit == "mm")
    return eVal * double(0.001);
  if (unit == "km")
    return eVal * double(1000);
  std::cerr << "We should never reach that stage in GetUnitInMeter\n";
  throw TerminalException{1};
}

TotalArrGetData RetrieveTotalArr(TripleModelDesc const &eTriple) {
  GridArray GrdArr = RETRIEVE_GRID_ARRAY(eTriple);
  ArrayHistory eArr = ReadArrayHistory(eTriple);
  if (eTriple.ModelName == "ROMS" || eTriple.ModelName == "ROMS_IVICA") {
    int iFile = 0;
    std::string eFile = ARR_GetHisFileName(eArr, "irrelevant", iFile);
    GrdArr.ARVD = ReadROMSverticalStratification(eFile);
  }
  return {std::move(GrdArr), std::move(eArr)};
}

#endif
