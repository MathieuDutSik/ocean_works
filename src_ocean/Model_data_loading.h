// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_OCEAN_MODEL_DATA_LOADING_H_
#define SRC_OCEAN_MODEL_DATA_LOADING_H_

#include "Model_grids.h"
#include "ROMSfunctionality.h"
#include "SmoothingBathyBasic.h"
#include <algorithm>
#include <string>
#include <utility>
#include <vector>

struct PairMinMax {
  double TheMin;
  double TheMax;
};

PairMinMax ComputeMinMax(GridArray const &GrdArr, MyMatrix<double> const &F) {
  //  std::cerr << "Begin of ComputeMinMax\n";
  bool IsFirst = true;
  int eta_rho = F.rows();
  int xi_rho = F.cols();
  int eta_rho_msk = GrdArr.GrdArrRho.MSK.rows();
  int xi_rho_msk = GrdArr.GrdArrRho.MSK.cols();
  if (eta_rho != eta_rho_msk || xi_rho != xi_rho_msk) {
    std::cerr << "ComputeMinMax error : Inconsistency in dimension\n";
    std::cerr << "  F: eta_rho=" << eta_rho << " xi_rho=" << xi_rho << "\n";
    std::cerr << "MSK: eta_rho=" << eta_rho_msk << " xi_rho=" << xi_rho_msk
              << "\n";
    std::cerr << "Most likely the grid does not match the chosen history\n";
    throw TerminalException{1};
  }
  double TheMin = 0;
  double TheMax = 0;
  for (int i = 0; i < eta_rho; i++)
    for (int j = 0; j < xi_rho; j++) {
      int eMsk;
      if (GrdArr.IsFE == 1) {
        eMsk = 1;
      } else {
        eMsk = GrdArr.GrdArrRho.MSK(i, j);
      }
      if (eMsk == 1) {
        double eVal = F(i, j);
        if (IsFirst) {
          TheMin = eVal;
          TheMax = eVal;
        } else {
          if (eVal < TheMin)
            TheMin = eVal;
          if (eVal > TheMax)
            TheMax = eVal;
        }
        IsFirst = false;
      }
    }
  //  std::cerr << "  End of ComputeMinMax\n";
  return {TheMin, TheMax};
}

PairMinMax ComputeMinMaxMask(MyMatrix<uint8_t> const &MSK,
                             MyMatrix<double> const &F) {
  bool IsFirst = true;
  int eta_rho = F.rows();
  int xi_rho = F.cols();
  int eta_rho_msk = MSK.rows();
  int xi_rho_msk = MSK.cols();
  if (eta_rho != eta_rho_msk || xi_rho != xi_rho_msk) {
    std::cerr << "ComputeMinMax error : Inconsistency in dimension\n";
    std::cerr << "  F: eta_rho=" << eta_rho << " xi_rho=" << xi_rho << "\n";
    std::cerr << "MSK: eta_rho=" << eta_rho_msk << " xi_rho=" << xi_rho_msk
              << "\n";
    std::cerr << "Most likely the grid does not match the chosen history\n";
    throw TerminalException{1};
  }
  double TheMin = 0;
  double TheMax = 0;
  for (int i = 0; i < eta_rho; i++)
    for (int j = 0; j < xi_rho; j++) {
      int eMsk = MSK(i, j);
      if (eMsk == 1) {
        double eVal = F(i, j);
        if (IsFirst) {
          TheMin = eVal;
          TheMax = eVal;
        } else {
          if (eVal < TheMin)
            TheMin = eVal;
          if (eVal > TheMax)
            TheMax = eVal;
        }
        IsFirst = false;
      }
    }
  return {TheMin, TheMax};
}

Eigen::Tensor<double, 3>
ComputeNormPairOfTensor(Eigen::Tensor<double, 3> const &Uthree,
                        Eigen::Tensor<double, 3> const &Vthree) {
  auto LDim = Uthree.dimensions();
  int dim0 = LDim[0];
  int dim1 = LDim[1];
  int dim2 = LDim[2];
  Eigen::Tensor<double, 3> Tens3(dim0, dim1, dim2);
  for (int i0 = 0; i0 < dim0; i0++)
    for (int i1 = 0; i1 < dim1; i1++)
      for (int i2 = 0; i2 < dim2; i2++) {
        double eU = Uthree(i0, i1, i2);
        double eV = Vthree(i0, i1, i2);
        double eNorm = sqrt(eU * eU + eV * eV);
        Tens3(i0, i1, i2) = eNorm;
      }
  return Tens3;
}

PairMinMax ComputeMinMax_3D(GridArray const &GrdArr,
                            Eigen::Tensor<double, 3> const &F) {
  //  std::cerr << "Begin of ComputeMinMax\n";
  auto LDim = F.dimensions();
  int nVert = LDim[0];
  int eta_rho = LDim[1];
  int xi_rho = LDim[2];
  int eta_rho_msk = GrdArr.GrdArrRho.MSK.rows();
  int xi_rho_msk = GrdArr.GrdArrRho.MSK.cols();
  if (eta_rho != eta_rho_msk || xi_rho != xi_rho_msk) {
    std::cerr << "ComputeMinMax_3D error : Inconsistency in dimension\n";
    std::cerr << "  F: eta_rho=" << eta_rho << " xi_rho=" << xi_rho << "\n";
    std::cerr << "MSK: eta_rho=" << eta_rho_msk << " xi_rho=" << xi_rho_msk
              << "\n";
    std::cerr << "Most likely the grid does not match the chosen history\n";
    throw TerminalException{1};
  }
  bool IsFirst = true;
  double TheMin = 0;
  double TheMax = 0;
  for (int i = 0; i < eta_rho; i++)
    for (int j = 0; j < xi_rho; j++) {
      int eMsk;
      if (GrdArr.IsFE == 1)
        eMsk = 1;
      else
        eMsk = GrdArr.GrdArrRho.MSK(i, j);
      if (eMsk == 1) {
        for (int iVert = 0; iVert < nVert; iVert++) {
          double eVal = F(iVert, i, j);
          if (IsFirst) {
            TheMin = eVal;
            TheMax = eVal;
          } else {
            if (eVal < TheMin)
              TheMin = eVal;
            if (eVal > TheMax)
              TheMax = eVal;
          }
          IsFirst = false;
        }
      }
    }
  //  std::cerr << "  End of ComputeMinMax\n";
  return {TheMin, TheMax};
}

MyMatrix<double> COMPUTE_NORM(MyMatrix<double> const &U,
                              MyMatrix<double> const &V) {
  int eta = U.rows();
  int xi = U.cols();
  MyMatrix<double> WINDMAG(eta, xi);
  for (int i = 0; i < eta; i++)
    for (int j = 0; j < xi; j++) {
      double eU = U(i, j);
      double eV = V(i, j);
      WINDMAG(i, j) = sqrt(eU * eU + eV * eV);
    }
  return WINDMAG;
}

MyMatrix<double> FreqPeriodChange(MyMatrix<double> const &F) {
  int eta = F.rows();
  int xi = F.cols();
  MyMatrix<double> Fret(eta, xi);
  double pi = 3.1415926535;
  for (int i = 0; i < eta; i++)
    for (int j = 0; j < xi; j++) {
      double eVal = F(i, j);
      double NewVal = 2 * pi / eVal;
      Fret(i, j) = NewVal;
    }
  return Fret;
}

bool TOTALARR_IsVar(TotalArrGetData const &TotalArr, std::string const &eVar) {
  if (TotalArr.eArr.KindArchive == "NETCDF") {
    std::string HisFile = ARR_GetFirstFileName(TotalArr.eArr);
    return NC_IsVar(HisFile, eVar);
  }
  if (TotalArr.eArr.KindArchive == "GRIB") {
    for (auto &eVarName : TotalArr.eArr.RawVarNames)
      if (eVarName == eVar)
        return true;
    return false;
  }
  std::cerr << "Error in TOTALARR_IsVar\n";
  std::cerr << "The KindArchive does not allow to find the nature\n";
  std::cerr << "KindArchive=" << TotalArr.eArr.KindArchive << "\n";
  throw TerminalException{1};
}

std::string TOTALARR_GetUnits(TotalArrGetData const &TotalArr,
                              std::string const &eVar) {
  if (TotalArr.eArr.KindArchive == "NETCDF") {
    std::string HisFile = ARR_GetFirstFileName(TotalArr.eArr);
    netCDF::NcFile dataFile(HisFile, netCDF::NcFile::read);
    netCDF::NcVar data = dataFile.getVar(eVar);
    netCDF::NcVarAtt eAtt = data.getAtt("units");
    char eString[1024] = "";
    eAtt.getValues(eString);
    return std::string(eString);
  }
  if (TotalArr.eArr.KindArchive == "GRIB") {
    for (auto &eRec : TotalArr.eArr.ListListMessages[0])
      if (eRec.shortName == eVar)
        return eRec.units;
    std::cerr << "Failed to find the entry\n";
    throw TerminalException{1};
  }
  std::cerr << "Error in TOTALARR_GetUnits\n";
  std::cerr << "The KindArchive does not allow to find the nature\n";
  std::cerr << "KindArchive=" << TotalArr.eArr.KindArchive << "\n";
  throw TerminalException{1};
}

MyMatrix<double> Get2DvariableSpecTime(TotalArrGetData const &TotalArr,
                                       std::string const &VarName,
                                       double const &eTimeDay) {
  if (TotalArr.eArr.KindArchive == "NETCDF") {
    //    std::cerr << "Before call to NETCDF_Get2DvariableSpecTime\n";
    return NETCDF_Get2DvariableSpecTime(TotalArr, VarName, eTimeDay);
  }
  if (TotalArr.eArr.KindArchive == "GRIB") {
    //    std::cerr << "Before call to GRIB_Get2DvariableSpecTime\n";
    return GRIB_Get2DvariableSpecTime(TotalArr, "shortName", VarName, eTimeDay);
  }
  std::cerr << "The KindArchive does not allow to find the nature\n";
  std::cerr << "KindArchive=" << TotalArr.eArr.KindArchive << "\n";
  throw TerminalException{1};
}

std::string DecltypeString(std::string const &FullVarName) {
  std::vector<std::string> ListStr = STRING_Split(FullVarName, ":");
  return ListStr[0];
}

// Copied from rho_eos.F in ROMS
Eigen::Tensor<double, 3>
ComputeDensityAnomaly(Eigen::Tensor<double, 3> const &TsArr,
                      Eigen::Tensor<double, 3> const &TtArr,
                      Eigen::Tensor<double, 3> const &z_rArr)
{
  //  std::cerr << "|zeta|=" << zeta.rows() << " / " << zeta.cols() << "\n";
  //  std::cerr << "|TsArr|=" << TsArr.dimension(0)  << " / " <<
  //  TsArr.dimension(1)  << " / " << TsArr.dimension(2)  << "\n"; std::cerr <<
  //  "|TtArr|=" << TtArr.dimension(0)  << " / " << TtArr.dimension(1)  << " / "
  //  << TtArr.dimension(2)  << "\n";
  double A00 = +1.909256e+04;
  double A01 = +2.098925e+02;
  double A02 = -3.041638e+00;
  double A03 = -1.852732e-03;
  double A04 = -1.361629e-05;
  double B00 = +1.044077e+02;
  double B01 = -6.500517e+00;
  double B02 = +1.553190e-01;
  double B03 = +2.326469e-04;
  double D00 = -5.587545e+00;
  double D01 = +7.390729e-01;
  double D02 = -1.909078e-02;
  double E00 = +4.721788e-01;
  double E01 = +1.028859e-02;
  double E02 = -2.512549e-04;
  double E03 = -5.939910e-07;
  double F00 = -1.571896e-02;
  double F01 = -2.598241e-04;
  double F02 = +7.267926e-06;
  double G00 = +2.042967e-03;
  double G01 = +1.045941e-05;
  double G02 = -5.782165e-10;
  double G03 = +1.296821e-07;
  double H00 = -2.595994e-07;
  double H01 = -1.248266e-09;
  double H02 = -3.508914e-09;
  double Q00 = +9.99842594e+02;
  double Q01 = +6.793952e-02;
  double Q02 = -9.095290e-03;
  double Q03 = +1.001685e-04;
  double Q04 = -1.120083e-06;
  double Q05 = +6.536332e-09;
  double U00 = +8.24493e-01;
  double U01 = -4.08990e-03;
  double U02 = +7.64380e-05;
  double U03 = -8.24670e-07;
  double U04 = +5.38750e-09;
  double V00 = -5.72466e-03;
  double V01 = +1.02270e-04;
  double V02 = -1.65460e-06;
  double W00 = +4.8314e-04;
  //
  auto LDim = TsArr.dimensions();
  int nVert = LDim[0];
  int eta_rho = LDim[1];
  int xi_rho = LDim[2];
  //  std::cerr << "Before ROMS_ComputeVerticalGlobalCoordinate\n";
  MyVector<double> C(10);
  Eigen::Tensor<double, 3> retTens(nVert, eta_rho, xi_rho);
  for (int k = 0; k < nVert; k++)
    for (int i = 0; i < eta_rho; i++)
      for (int j = 0; j < xi_rho; j++) {
        double Tt = TtArr(k, i, j);
        double Ts = TsArr(k, i, j);
        double sqrtTs = sqrt(Ts);
        double Tp = z_rArr(k, i, j);
        double Tpr10 = Tp * 0.1;
        C(0) =
            Q00 + Tt * (Q01 + Tt * (Q02 + Tt * (Q03 + Tt * (Q04 + Tt * Q05))));
        C(1) = U00 + Tt * (U01 + Tt * (U02 + Tt * (U03 + Tt * U04)));
        C(2) = V00 + Tt * (V01 + Tt * V02);
        //
        double den1 = C(0) + Ts * (C(1) + sqrtTs * C(2) + Ts * W00);
        C(3) = A00 + Tt * (A01 + Tt * (A02 + Tt * (A03 + Tt * A04)));
        C(4) = B00 + Tt * (B01 + Tt * (B02 + Tt * B03));
        C(5) = D00 + Tt * (D01 + Tt * D02);
        C(6) = E00 + Tt * (E01 + Tt * (E02 + Tt * E03));
        C(7) = F00 + Tt * (F01 + Tt * F02);
        C(8) = G01 + Tt * (G02 + Tt * G03);
        C(9) = H00 + Tt * (H01 + Tt * H02);
        //
        double bulk0 = C(3) + Ts * (C(4) + sqrtTs * C(5));
        double bulk1 = C(6) + Ts * (C(7) + sqrtTs * G00);
        double bulk2 = C(8) + Ts * C(9);
        double bulk = bulk0 - Tp * (bulk1 - Tp * bulk2);
        //
        double cff = 1.0 / (bulk + Tpr10);
        double den = den1 * bulk * cff - 1000;
        //
        retTens(k, i, j) = den;
      }
  return retTens;
}

Eigen::Tensor<double, 3>
ROMS_ComputeDensityAnomaly(Eigen::Tensor<double, 3> const &TsArr,
                      Eigen::Tensor<double, 3> const &TtArr,
                      GridArray const &GrdArr, MyMatrix<double> const &zeta)
{
  Eigen::Tensor<double, 3> z_rArr =
      ROMS_ComputeVerticalGlobalCoordinate_r(GrdArr, zeta);
  return ComputeDensityAnomaly(TsArr, TtArr, z_rArr);
}

// According to Nastjenka Supic
// Dynamic depths (D) should be I think computed for both „disturbed“ (with LNG flow) and „undisturbed“ case.
// D=2,5*100000 / verticaly average of 0-25 m density
MyMatrix<double> GetDynamicDepth(MyMatrix<double> const &Avg25_Density) {
  int n_rows = Avg25_Density.rows();
  int n_cols = Avg25_Density.cols();
  MyMatrix<double> Mret(n_rows, n_cols);
  for (int i = 0; i < n_rows; i++) {
    for (int j = 0; j < n_cols; j++) {
      double eAvg = Avg25_Density(i,j);
      double eDynDep = 2.5 * 100000 / eAvg;
      Mret(i, j) = eDynDep;
    }
  }
  return Mret;
}

Eigen::Tensor<double, 3>
GetDensity(Eigen::Tensor<double, 3> const &DensAnomaly) {
  auto LDim = DensAnomaly.dimensions();
  int nVert = LDim[0];
  int eta_rho = LDim[1];
  int xi_rho = LDim[2];
  //  std::cerr << "Before ROMS_ComputeVerticalGlobalCoordinate\n";
  MyVector<double> C(10);
  Eigen::Tensor<double, 3> retDensity(nVert, eta_rho, xi_rho);
  for (int k = 0; k < nVert; k++) {
    for (int i = 0; i < eta_rho; i++) {
      for (int j = 0; j < xi_rho; j++) {
        retDensity(k, i, j) = 1000 + DensAnomaly(k, i, j);
      }
    }
  }
  return retDensity;
}




// From Ivica code.
MyMatrix<double> mixing_ratio2relative_humidity(MyMatrix<double> const &Q2,
                                                MyMatrix<double> const &PSFC,
                                                MyMatrix<double> const &T2K) {
  int eta = Q2.rows();
  int xi = Q2.cols();
  MyMatrix<double> retField(eta, xi);
  double t_sat = 373.16;
  double c1 = 7.90298;
  double c2 = 5.02808;
  double c3 = 1.3816e-7;
  double c4 = 11.344;
  double c5 = 8.1328e-3;
  double c6 = 3.49149;
  double water_air_mass_ratio = 0.621970585;
  for (int i = 0; i < eta; i++)
    for (int j = 0; j < xi; j++) {
      double mix_ratio = Q2(i, j) * 1000;
      double pres = PSFC(i, j) / 100;
      double temp = T2K(i, j);
      double t_ratio = t_sat / temp;
      double rt_ratio = 1 / t_ratio;
      double sl_pressure = 1013.246;
      double svp1 = -1.0 * c1 * (t_ratio - 1.0) + c2 * log10(t_ratio) -
                    c3 * (pow(10.0, c4 * (1.0 - rt_ratio)) - 1.0) +
                    c5 * (pow(10.0, -1.0 * c6 * (t_ratio - 1.0)) - 1.0) +
                    log10(sl_pressure);
      double svp = pow(10.0, svp1);
      double vapor_pressure_ratio = svp / (pres - svp);
      double smr = 1000.0 * water_air_mass_ratio * vapor_pressure_ratio;
      double rh = (mix_ratio / smr) * 100;
      retField(i, j) = rh;
    }
  return retField;
}

// A number of algorithms for computing relative humidity
MyMatrix<double> Algorithms_RelativeHumidity(TotalArrGetData const &TotalArr,
                                             std::string const &VarName2d,
                                             std::string const &VarName2r,
                                             std::string const &VarName2t,
                                             std::string const &VarNameMSL,
                                             std::string const &VarNameQ,
                                             double const &eTimeDay) {
  // directly available relative humidity. Returning it.
  if (TOTALARR_IsVar(TotalArr, VarName2r)) {
    return Get2DvariableSpecTime(TotalArr, VarName2r, eTimeDay);
  }
  //
  // Second try: Using Specific humidity to convert to relative humidity
  //
  /* This algo seems not finished and/or unreliable
    double airDens=1.225; // kg/m3 : Density of air, see
    https://en.wikipedia.org/wiki/Density_of_air double waterDens=0.804; // g/L
    = kg/m3 : Density of water vapor https://en.wikipedia.org/wiki/Water_vapor
    // We use formula from
    http://www.engineeringtoolbox.com/specific-relative-humidity-air-d_688.html
    // formula is then phi = 100 x /(0.622 * rho_ws/(rho - rho_ws) )
    // It seems not to work correctly.
    double quot=waterDens / (airDens - waterDens);
    MyMatrix<double> Fspecific=Get2DvariableSpecTime(TotalArr, VarNameQ,
    eTimeDay); double TheMult=100 / ( 0.622 * quot); F = Fspecific * TheMult; */
  if (TOTALARR_IsVar(TotalArr, VarNameQ) &&
      TOTALARR_IsVar(TotalArr, VarName2t)) {
    // We use formula from
    // http://earthscience.stackexchange.com/questions/2360/how-do-i-convert-specific-humidity-to-relative-humidity
    MyMatrix<double> F_q = Get2DvariableSpecTime(TotalArr, VarNameQ, eTimeDay);
    MyMatrix<double> F_p;
    if (TOTALARR_IsVar(TotalArr, VarNameMSL)) {
      F_p = Get2DvariableSpecTime(TotalArr, VarNameMSL, eTimeDay);
    } else {
      int eta = F_q.rows();
      int xi = F_q.cols();
      F_p.setConstant(eta, xi, 103000);
    }
    MyMatrix<double> F_TK =
        Get2DvariableSpecTime(TotalArr, VarName2t, eTimeDay);
    int eta = F_q.rows();
    int xi = F_q.cols();
    MyMatrix<double> F(eta, xi);
    for (int i = 0; i < eta; i++)
      for (int j = 0; j < xi; j++) {
        double eT = F_TK(i, j);
        double eQ = F_q(i, j);
        double eP = F_p(i, j);
        double eT0 = static_cast<double>(273.15);
        double TheQuot = static_cast<double>(17.67) * (eT - eT0) /
                         (eT - static_cast<double>(29.65));
        double eRH = 0.263 * eP * eQ / (exp(TheQuot));
        F(i, j) = std::min(eRH, static_cast<double>(100));
      }
    return F;
  }
  //
  // Using 2m_dew temperature + 2m temperature
  //
  // Following
  // https://github.com/dcherian/tools/blob/master/ROMS/arango/forcing/d_ecmwf2roms.m
  if (TOTALARR_IsVar(TotalArr, VarName2d) &&
      TOTALARR_IsVar(TotalArr, VarName2t)) {
    MyMatrix<double> F_tsur =
        Get2DvariableSpecTime(TotalArr, VarName2t, eTimeDay);
    MyMatrix<double> F_tdew =
        Get2DvariableSpecTime(TotalArr, VarName2d, eTimeDay);
    int eta = F_tsur.rows();
    int xi = F_tsur.cols();
    MyMatrix<double> F(eta, xi);
    for (int i = 0; i < eta; i++)
      for (int j = 0; j < xi; j++) {
        double tsur = F_tsur(i, j) - 273.15;
        double tdew = F_tdew(i, j) - 273.15;
        double E = 6.11 * std::pow(10.0, 7.5 * tdew / (237.7 + tdew));
        double Es = 6.11 * std::pow(10.0, 7.5 * tsur / (237.7 + tsur));
        double field = 100.0 * (E / Es);
        F(i, j) = field;
      }
    return F;
  }
  //
  // No algorithm found
  //
  std::cerr << "We could not find a relevant algorithm for relative humidity\n";
  std::cerr << "Following were tried\n";
  std::cerr << "--- Direct: Variable from the model\n";
  std::cerr << "--- Conversion from specific humidity to relative humidity\n";
  std::cerr << "--- Using 2m temperature and 2m dew temperature\n";
  std::cerr << "Please provide needed variables (or another algorithm)\n";
  throw TerminalException{1};
}

template <typename T> void RemoveNegativeValues(MyMatrix<T> &M) {
  size_t n_rows = M.rows();
  size_t n_cols = M.cols();
  for (size_t i = 0; i < n_rows; i++)
    for (size_t j = 0; j < n_cols; j++)
      if (M(i, j) < 0)
        M(i, j) = 0;
}

template <typename T>
void SetNegativeDeepNight(MyMatrix<T> &M, double const &eMJD) {
  std::vector<int> eVect = DATE_ConvertMjd2six(eMJD);
  int hour = eVect[3];
  if (hour <= 4 || hour >= 22) {
    size_t n_rows = M.rows();
    size_t n_cols = M.cols();
    for (size_t i = 0; i < n_rows; i++)
      for (size_t j = 0; j < n_cols; j++)
        M(i, j) = 0;
  }
}

template <typename T> void SetSmallToZero(MyMatrix<T> &M, T const &smallVal) {
  size_t n_rows = M.rows();
  size_t n_cols = M.cols();
  std::cerr << "SetSmallToZero, begin max=" << M.maxCoeff() << "\n";
  for (size_t i = 0; i < n_rows; i++)
    for (size_t j = 0; j < n_cols; j++)
      if (T_abs(M(i, j)) < smallVal)
        M(i, j) = 0;
  std::cerr << "SetSmallToZero, end max=" << M.maxCoeff() << "\n";
}

std::string GetBasicModelName(std::string const &eModelName) {
  if (eModelName == "ROMS_IVICA")
    return "ROMS";
  if (eModelName == "WWM_DAILY")
    return "WWM";
  return eModelName;
}

std::vector<std::string> Get_BFM_vars() {
  std::vector<std::string> LVar_BFM{"oxygen",
                                    "PO4",
                                    "NO3",
                                    "NH4",
                                    "NitrogenSink",
                                    "SiOH4",
                                    "ReductionEquivalent",
                                    "bacteriaC",
                                    "bacteriaN",
                                    "bacteriaP",
                                    "diatomsC",
                                    "diatomsN",
                                    "diatomsP",
                                    "diatomsL",
                                    "diatomsS",
                                    "flagellatesC",
                                    "flagellatesN",
                                    "flagellatesP",
                                    "flagellatesL",
                                    "picophytoplanktonC",
                                    "picophytoplanktonN",
                                    "picophytoplanktonP",
                                    "picophytoplanktonL",
                                    "largephytoplanktonC",
                                    "largephytoplanktonN",
                                    "largephytoplanktonP",
                                    "largephytoplanktonL",
                                    "CarnPlanktonC",
                                    "CarnPlanktonN",
                                    "CarnPlanktonP",
                                    "OmniPlanktonC",
                                    "OmniPlanktonN",
                                    "OmniPlanktonP",
                                    "MicroPlanktonC",
                                    "MicroPlanktonN",
                                    "MicroPlanktonP",
                                    "HeteroNanoflagelattesC",
                                    "HeteroNanoflagelattesN",
                                    "HeteroNanoflagelattesP",
                                    "LabileDOM1c",
                                    "LabileDOM1n",
                                    "LabileDOM1p",
                                    "LabileDOM2c",
                                    "RefractoryDOMc",
                                    "ParticleOMc",
                                    "DissolvedOMc",
                                    "ParticleOMn",
                                    "ParticleOMp",
                                    "ParticleOMs",
                                    "DissolvedICc",
                                    "DissolvedICh",
                                    "Irradiance",
                                    "DIC",
                                    "chlorophyll",
                                    "NetProductionP1",
                                    "NetProductionP2",
                                    "NetProductionP3",
                                    "NetProductionP4",
                                    "RegFactorP1",
                                    "RegFactorP2",
                                    "RegFactorP3",
                                    "RegFactorP4",
                                    "GrossPP",
                                    "SecondPP",
                                    "ExtinctionCoeff"};
  return LVar_BFM;
}

std::vector<std::string> GetAllPossibleVariables() {
  std::vector<std::string> ListVarOut{"IOBP",
                                      "MAPSTA",
                                      "FieldOut1",
                                      "CFL1",
                                      "CFL2",
                                      "CFL3",
                                      "pH",
                                      "alkalinity",
                                      "ThreeDfield1",
                                      "NbIterSolv",
                                      "WIND10",
                                      "Uwind",
                                      "Vwind",
                                      "WINDMAG",
                                      "Vorticity",
                                      "ChlorophyllConcOCI",
                                      "ChlorophyllConcOCX",
                                      "CalciteConc",
                                      "ParticleOrganicCarbon",
                                      "InstPhotosyntheticallyAvailableRad",
                                      "PhotosyntheticallyAvailableRad",
                                      "SurfStress",
                                      "SurfCurr",
                                      "SurfCurrMag",
                                      "Curr",
                                      "CurrMag",
                                      "HorizCurr",
                                      "ROMSbarotropicdefect",
                                      "CurrBaro",
                                      "CurrBaroMag",
                                      "ChlorophylA",
                                      "Temp",
                                      "Salt",
                                      "TempSurf",
                                      "DensAnomaly",
                                      "Density",
                                      "DynamicDepth",
                                      "AIRT2",
                                      "AIRT2K",
                                      "Rh2",
                                      "Rh2frac",
                                      "AIRD",
                                      "SurfPres",
                                      "ZetaOcean",
                                      "ZetaOceanDerivative",
                                      "DynBathy",
                                      "Bathymetry",
                                      "RoughnessFactor",
                                      "ZetaSetup",
                                      "phytoplanktonbiomass",
                                      "Dye1",
                                      "CdWave",
                                      "AlphaWave",
                                      "AirZ0",
                                      "AirFricVel",
                                      "CGwave",
                                      "shflux",
                                      "ssflux",
                                      "evaporation",
                                      "CloudFraction",
                                      "Hwave",
                                      "BreakingFraction",
                                      "rain",
                                      "swrad",
                                      "lwrad",
                                      "latent",
                                      "sensible",
                                      "SourceGainDustAerosolSML",
                                      "SourceGainDustAerosolSmall",
                                      "SourceGainDustAerosolMedium",
                                      "SourceGainDustAerosolLarge",
                                      "aermssdus",
                                      "aermssdum",
                                      "aermssdul",
                                      "aermssduSML",
                                      "aermssomhphil",
                                      "aermssomhphob",
                                      "MeanWaveFreq",
                                      "PeakWaveFreq",
                                      "TM02",
                                      "MeanWavePer",
                                      "PeakWavePer",
                                      "MeanWaveDirSpread",
                                      "PeakWaveDirSpread",
                                      "MeanWaveDir",
                                      "PeakWaveDir",
                                      "MeanWaveDirVect",
                                      "PeakWaveDirVect",
                                      "DiscPeakWaveDir",
                                      "MeanWaveLength",
                                      "PeakWaveLength",
                                      "MeanWaveNumber",
                                      "PeakWaveNumber",
                                      "TotSurfStr",
                                      "WaveSurfStr",
                                      "SurfStrHF"};
  for (auto &eVar : Get_BFM_vars()) {
    ListVarOut.push_back(eVar);
  }
  return ListVarOut;
}

Eigen::Tensor<double, 3>
RetrieveStandardVerticalCoordinate(TotalArrGetData const &TotalArr) {
  std::string eModelName = GetBasicModelName(TotalArr.GrdArr.ModelName);
  if (eModelName == "SCHISM_NETCDF_OUT") {
    double eTimeDay = MinimumTimeHistoryArray(TotalArr.eArr);
    Eigen::Tensor<double, 3> znl =
        NETCDF_Get3DvariableSpecTime(TotalArr, "znl", eTimeDay);
    return znl;
  }
  if (eModelName == "ROMS") {
    int eta_rho = TotalArr.GrdArr.GrdArrRho.LON.rows();
    int xi_rho = TotalArr.GrdArr.GrdArrRho.LON.cols();
    MyMatrix<double> zeta = ZeroMatrix<double>(eta_rho, xi_rho);
    return ROMS_ComputeVerticalGlobalCoordinate_r(TotalArr.GrdArr, zeta);
  }
  std::cerr << "Error in RetrieveStandardVerticalCoordinate\n";
  std::cerr << "eModelName=" << eModelName << "\n";
  std::cerr << "Programmed models : SCHISM_NETCDF_OUT and ROMS\n";
  std::cerr << "We cannot find a matching model\n";
  throw TerminalException{1};
}

struct VerticalLevelInfo {
  int Choice;
  double dep1;
  double dep2;
  std::string strNewVarName;
  std::string strDepth;
  std::string type;
};

// Use is
// HorizTemp:VR-2m for the temperature at 2m depth
// HorizTemp:VA-2m for the temperature at 2m depth
// VR stands for Relative depth (according to the zeta)
// VA stands for absolute depth (according to the geoide)
VerticalLevelInfo RetrieveVerticalInformation(std::string const &FullVarName,
                                              std::string const &eModelName) {
  int Choice = -1;
  std::string strDepth = "unset", type = "unset";
  double dep1 = -1000000;
  double dep2 = -1000000;
  std::string separator1 = ":";
  std::vector<std::string> ListStr = STRING_Split(FullVarName, separator1);
  std::string eVarName = ListStr[0];
  if (ListStr.size() != 2) {
    std::cerr << "Error in the variable name. ListStr should be of length 2\n";
    std::cerr << "FullVarName=" << FullVarName << "\n";
    std::cerr << "eVarName=" << eVarName << "\n";
    std::cerr << "separator1=" << separator1 << "\n";
    std::cerr << "eModelName=" << eModelName << "\n";
    std::cerr << "|ListStr|=" << ListStr.size() << "\n";
    for (size_t i_str = 0; i_str < ListStr.size(); i_str++)
      std::cerr << "i_str=" << i_str << " e_str=" << ListStr[i_str] << "\n";
    throw TerminalException{1};
  }
  // str should be "VA-2m" or "VR-2m"
  std::string str = ListStr[1];
  // ListStrB should be "VA", "-2", "m"
  std::vector<std::string> ListStrB = STRING_SplitCharNb(str);
  if (ListStrB.size() != 1 && ListStrB.size() != 3 && ListStrB.size() != 5) {
    std::cerr << "Error in the variable name should have 1, 3 or 5 blocks\n";
    std::cerr << "|ListStrB|=" << ListStrB.size() << "\n";
    std::cerr << "str=" << str << "\n";
    int siz = ListStrB.size();
    for (int j = 0; j < siz; j++)
      std::cerr << "j=" << j << " ListStrB[j]=" << ListStrB[j] << "\n";
    throw TerminalException{1};
  }
  if (ListStrB[0] == "VA")
    Choice = 1;
  if (ListStrB[0] == "VR")
    Choice = 2;
  if (ListStrB[0] == "TOTAVE")
    Choice = 3;
  if (ListStrB[0] == "BOT")
    Choice = 4;
  if (ListStrB[0] == "SURF")
    Choice = 5;
  if (ListStrB[0] == "AVE")
    Choice = 6;
  if (Choice == -1) {
    std::cerr << "We should have VA or VR as possible choice\n";
    throw TerminalException{1};
  }
  std::string strNewVarName;
  if (Choice == 1 || Choice == 2) {
    if (ListStrB.size() != 3) {
      std::cerr << "We should have |ListStrB|=3 for Choice=1 or 2\n";
      std::cerr << "|ListStrB|=" << ListStrB.size() << "\n";
      for (size_t i=0; i<ListStrB.size(); i++)
        std::cerr << "i=" << i << " str=" << ListStrB[i] << "\n";
      throw TerminalException{1};
    }
    double eVal;
    std::istringstream(ListStrB[1]) >> eVal;
    dep1 = GetUnitInMeter(eVal, ListStrB[2]);
    strDepth = " at " + ListStr[1];
    strNewVarName = eVarName + strDepth;
  }
  if (Choice == 3) {
    strDepth = " vertical average";
    strNewVarName = "Vertical average " + eVarName;
  }
  if (Choice == 4) {
    strDepth = " bottom";
    strNewVarName = "bottom " + eVarName;
  }
  if (Choice == 5) {
    strDepth = " surface";
    strNewVarName = "surface " + eVarName;
  }
  if (Choice == 6) {
    if (ListStrB.size() != 5) {
      std::cerr << "We should have |ListStrB|=5 for Choice=6\n";
      throw TerminalException{1};
    }
    double eVal1;
    std::istringstream(ListStrB[1]) >> eVal1;
    dep1 = GetUnitInMeter(eVal1, ListStrB[2]);
    double eVal2;
    std::istringstream(ListStrB[3]) >> eVal2;
    dep2 = GetUnitInMeter(eVal2, ListStrB[4]);
    std::string strDep1 = ListStrB[1] + ListStrB[2];
    std::string strDep2 = ListStrB[3] + ListStrB[4];
    strDepth = " from " + strDep1 + " to " + strDep2;
    strNewVarName = eVarName + strDepth;
  }
  //
  // Unclear as to what functionality is supported here
  if (eModelName == "WWM") {
    std::string str = ListStr[1];
    std::vector<std::string> ListStrB = STRING_SplitCharNb(str);
    std::string type = ListStrB[0];
    strDepth = ListStrB[1];
  }
  //
  return {Choice, dep1, dep2, strNewVarName, strDepth, type};
}

/* From a three dimensional array, compute the two-dimensional one.
   It can be several things:
   --- vertical level at absolute level
   --- vertical level at relative level
   --- vertical average of the level
*/
MyMatrix<double> ThreeDimensional_to_TwoDimensional(
    Eigen::Tensor<double, 3> const &F3, MyMatrix<double> const &zeta,
    TotalArrGetData const &TotalArr, VerticalLevelInfo const &VertInfo,
    double const &eTimeDay) {
  std::string eModelName = GetBasicModelName(TotalArr.GrdArr.ModelName);
  if (eModelName == "ROMS") {
    if (VertInfo.Choice == 1 || VertInfo.Choice == 2)
      return VerticalInterpolation_P2_R(
          TotalArr.GrdArr.ARVD, GetDEP(TotalArr.GrdArr.GrdArrRho), zeta,
          TotalArr.GrdArr.GrdArrRho.MSK, VertInfo.dep1, F3, VertInfo.Choice);
    if (VertInfo.Choice == 3)
      return ConvertBaroclinic_to_Barotropic(F3, zeta, TotalArr.GrdArr);
    if (VertInfo.Choice == 4)
      return DimensionExtraction(F3, 0, 0);
    if (VertInfo.Choice == 5) {
      int s_rho = F3.dimension(0);
      return DimensionExtraction(F3, 0, s_rho - 1);
    }
    if (VertInfo.Choice == 6) {
      return VerticalInterpolationAverage_P2_R(
          TotalArr.GrdArr.ARVD, GetDEP(TotalArr.GrdArr.GrdArrRho), zeta,
          TotalArr.GrdArr.GrdArrRho.MSK, VertInfo.dep1, VertInfo.dep2, F3);
    }
    std::cerr << "Failing to find matching entry for Choice\n";
    std::cerr << "Choice=" << VertInfo.Choice << "\n";
    throw TerminalException{1};
  }
  if (eModelName == "SCHISM_NETCDF_OUT") {
    Eigen::Tensor<double, 3> znl =
        NETCDF_Get3DvariableSpecTime(TotalArr, "znl", eTimeDay);
    return VerticalInterpolation_SCHISM_ZNL(znl, zeta, VertInfo.dep1, F3,
                                            VertInfo.Choice);
  }
  if (eModelName == "NEMO") {
    //    if (VertInfo.Choice == 3)
    //      return ConvertBaroclinic_to_Barotropic(F3, zeta, TotalArr.GrdArr);
    if (VertInfo.Choice == 4)
      return NEMO_GetDeepestField(F3);
    if (VertInfo.Choice == 5) {
      int s_rho = F3.dimension(0);
      return DimensionExtraction(F3, 0, s_rho - 1);
    }
    std::cerr << "Choice=" << VertInfo.Choice << " not supported apparently\n";
    throw TerminalException{1};
  }
  std::cerr << "Error in ThreeDimensional_to_TwoDimensional\n";
  std::cerr << "eModelName=" << eModelName << " not supported\n";
  throw TerminalException{1};
}

std::string TransformVarName(std::string const &str) {
  std::string RetStr;
  int len = str.size();
  for (int i = 0; i < len; i++) {
    std::string eChar = str.substr(i, 1);
    if (eChar == ":")
      eChar = "_";
    RetStr += eChar;
  }
  return RetStr;
}

MyMatrix<double> GetNormMatrix(MyMatrix<double> const &U,
                               MyMatrix<double> const &V) {
  int nbRow = U.rows();
  int nbCol = U.cols();
  MyMatrix<double> Fwr(nbRow, nbCol);
  for (int iRow = 0; iRow < nbRow; iRow++)
    for (int iCol = 0; iCol < nbCol; iCol++) {
      double eU = U(iRow, iCol);
      double eV = V(iRow, iCol);
      double eNorm = sqrt(eU * eU + eV * eV);
      Fwr(iRow, iCol) = eNorm;
    }
  return Fwr;
}

RecVar Average_RecVar(std::vector<RecVar> const &ListRecVar) {
  int nEnt = ListRecVar.size();
  RecVar RetRecVar = ListRecVar[0];
  // U
  if (RetRecVar.U.size() > 0) {
    for (int iEnt = 1; iEnt < nEnt; iEnt++)
      RetRecVar.U += ListRecVar[iEnt].U;
    RetRecVar.U /= nEnt;
  }
  // V
  if (RetRecVar.V.size() > 0) {
    for (int iEnt = 1; iEnt < nEnt; iEnt++)
      RetRecVar.V += ListRecVar[iEnt].V;
    RetRecVar.V /= nEnt;
  }
  // F
  if (RetRecVar.F.size() > 0) {
    for (int iEnt = 1; iEnt < nEnt; iEnt++)
      RetRecVar.F += ListRecVar[iEnt].F;
    RetRecVar.F /= nEnt;
  }
  // Uthree
  if (RetRecVar.Uthree.size() > 0) {
    for (int iEnt = 1; iEnt < nEnt; iEnt++)
      RetRecVar.Uthree += ListRecVar[iEnt].Uthree;
    for (int i = 0; i < RetRecVar.Uthree.size(); i++)
      RetRecVar.Uthree(i) /= nEnt;
  }
  // Vthree
  if (RetRecVar.Vthree.size() > 0) {
    for (int iEnt = 1; iEnt < nEnt; iEnt++)
      RetRecVar.Vthree += ListRecVar[iEnt].Vthree;
    for (int i = 0; i < RetRecVar.Vthree.size(); i++)
      RetRecVar.Vthree(i) /= nEnt;
  }
  // Tens3
  if (RetRecVar.Tens3.size() > 0) {
    for (int iEnt = 1; iEnt < nEnt; iEnt++)
      RetRecVar.Tens3 += ListRecVar[iEnt].Tens3;
    for (int i = 0; i < RetRecVar.Tens3.size(); i++)
      RetRecVar.Tens3(i) /= nEnt;
  }
  return RetRecVar;
}

RecSymbolic GetRecSymbolic(double const &eTimeDay,
                           std::string const &FullVarName) {
  std::string strPres = DATE_ConvertMjd2mystringPres(eTimeDay);
  std::string strFile = DATE_ConvertMjd2mystringFile(eTimeDay);
  RecSymbolic RecS;
  RecS.eTimeDay = eTimeDay;
  RecS.strPres = "at " + strPres;
  RecS.strFile = strFile;
  RecS.VarNature = "rho";
  RecS.FullVarName = FullVarName;
  RecS.VarName1 = TransformVarName(FullVarName);
  RecS.VarName2 = "unset";
  RecS.strTime_ROMS = "unset";
  RecS.varName_GRIB = "unset";
  return RecS;
}

RecVar ModelSpecificVarSpecificTime_Kernel(TotalArrGetData const &TotalArr,
                                           std::string const &FullVarName,
                                           double const &eTimeDay) {
  std::string eModelName = GetBasicModelName(TotalArr.GrdArr.ModelName);
  std::vector<std::string> ListStr = STRING_Split(FullVarName, ":");
  std::string eVarName = ListStr[0];
  //  std::cerr << "   ModelSpecificVarSpecificTime_Kernel, FullVarName=" <<
  //  FullVarName << "\n"; std::cerr << "   ModelSpecificVarSpecificTime_Kernel,
  //  eVarName=" << eVarName << " eTimeDay=" << eTimeDay << " eModelName=" <<
  //  eModelName << "\n";
  int eta_rho = TotalArr.GrdArr.GrdArrRho.LON.rows();
  int xi_rho = TotalArr.GrdArr.GrdArrRho.LON.cols();
  RecSymbolic RecS = GetRecSymbolic(eTimeDay, FullVarName);
  MyMatrix<double> F;
  MyMatrix<double> U;
  MyMatrix<double> V;
  Eigen::Tensor<double, 3> Tens3;
  Eigen::Tensor<double, 3> Uthree;
  Eigen::Tensor<double, 3> Vthree;
  //
  // Generic code for handling the interpolation from 3D to 2D
  //
  if (ListStr.size() != 1) {
    VerticalLevelInfo VertInfo =
        RetrieveVerticalInformation(FullVarName, eModelName);
    RecVar Rec3D =
        ModelSpecificVarSpecificTime_Kernel(TotalArr, eVarName, eTimeDay);
    if (Rec3D.RecS.VarNature != "3Drho" && Rec3D.RecS.VarNature != "3Duv") {
      std::cerr << "In order to apply the VR, VA, Surf transformations, the "
                   "variable should be 3D\n";
      throw TerminalException{1};
    }
    RecS = Rec3D.RecS;
    RecS.FullVarName = FullVarName;
    RecS.VarName1 = TransformVarName(FullVarName);
    RecS.VarName2 = VertInfo.strNewVarName;
    //
    //    std::cerr << "VarName1=" << RecS.VarName1 << "\n";
    //    std::cerr << "VarName2=" << RecS.VarName2 << "\n";
    auto get_reczeta = [&]() -> RecVar {
      try {
        return ModelSpecificVarSpecificTime_Kernel(TotalArr, "ZetaOcean",
                                                   eTimeDay);
      } catch (TerminalException const &e) {
        std::cerr
            << "We failed to obtain ZetaOcean, so instead, we set it to zero\n";
        RecVar RecZeta;
        RecZeta.RecS = RecS;
        RecZeta.F = ZeroMatrix<double>(eta_rho, xi_rho);
        return RecZeta;
      }
    };
    RecVar RecZeta = get_reczeta();
    if (Rec3D.RecS.VarNature == "3Drho") {
      F = ThreeDimensional_to_TwoDimensional(Rec3D.Tens3, RecZeta.F, TotalArr,
                                             VertInfo, eTimeDay);
      RecS.VarNature = "rho";
    }
    if (Rec3D.RecS.VarNature == "3Duv") {
      U = ThreeDimensional_to_TwoDimensional(Rec3D.Uthree, RecZeta.F, TotalArr,
                                             VertInfo, eTimeDay);
      V = ThreeDimensional_to_TwoDimensional(Rec3D.Vthree, RecZeta.F, TotalArr,
                                             VertInfo, eTimeDay);
      RecS.VarNature = "uv";
    }
  }
  //
  // Generic model kind of variables
  //
  if (FullVarName == "NbIterSolv") {
    if (eModelName == "WWM")
      F = Get2DvariableSpecTime(TotalArr, "NB_ITER_SOLV", eTimeDay);
    RecS.VarName2 = "nb Iteration Solver";
    RecS.minval = 0;
    RecS.maxval = 50;
    RecS.mindiff = -5;
    RecS.maxdiff = 5;
    RecS.Unit = "nondim.";
  }
  if (FullVarName == "Dye1") {
    if (eModelName == "ROMS") {
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "dye_01", eTimeDay);
    }
    RecS.VarName2 = "Dye concentration";
    RecS.minval = 0;
    RecS.maxval = 1;
    RecS.mindiff = -0.00001;
    RecS.maxdiff = 0.00001;
    RecS.VarNature = "3Drho";
    RecS.Unit = "nondim.e";
  }
  if (FullVarName == "CFL1") {
    if (eModelName == "WWM")
      F = Get2DvariableSpecTime(TotalArr, "CFL1", eTimeDay);
    RecS.VarName2 = "CFL1";
    RecS.minval = 0;
    RecS.maxval = 5;
    RecS.mindiff = -1;
    RecS.maxdiff = 1;
    RecS.Unit = "nondim.";
  }
  if (FullVarName == "CFL2") {
    if (eModelName == "WWM")
      F = Get2DvariableSpecTime(TotalArr, "CFL2", eTimeDay);
    RecS.VarName2 = "CFL2";
    RecS.minval = 0;
    RecS.maxval = 5;
    RecS.mindiff = -1;
    RecS.maxdiff = 1;
    RecS.Unit = "nondim.";
  }
  if (FullVarName == "CFL3") {
    if (eModelName == "WWM")
      F = Get2DvariableSpecTime(TotalArr, "CFL3", eTimeDay);
    RecS.VarName2 = "CFL3";
    RecS.minval = 0;
    RecS.maxval = 5;
    RecS.mindiff = -1;
    RecS.maxdiff = 1;
    RecS.Unit = "nondim.";
  }
  if (FullVarName == "ThreeDfield1") {
    if (eModelName == "WWM" || eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "ThreeDfield1", eTimeDay);
    RecS.VarName2 = "Generic three dim. field 1";
    RecS.minval = 0;
    RecS.maxval = 1;
    RecS.mindiff = -1;
    RecS.maxdiff = 1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "unspecified";
  }
  if (FullVarName == "CGwave") {
    if (eModelName == "WWM")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "CG", eTimeDay);
    RecS.VarName2 = "group velocity";
    RecS.minval = 0;
    RecS.maxval = 30;
    RecS.mindiff = -1;
    RecS.maxdiff = 1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "m/s";
  }
  if (FullVarName == "FieldOut1") {
    if (eModelName == "WWM")
      F = Get2DvariableSpecTime(TotalArr, "FieldOut1", eTimeDay);
    RecS.VarName2 = "Generic Field Out 1";
    RecS.minval = 0;
    RecS.maxval = 1;
    RecS.mindiff = -1;
    RecS.maxdiff = 1;
    RecS.Unit = "unspecified";
  }
  if (FullVarName == "ChlorophyllConcOCI") {
    RecS.VarName2 = "Chlorophyll Concentration, OCI Algorithm";
    RecS.CFshortName = "chlor_a";
    RecS.minval = 0;
    RecS.maxval = 100;
    RecS.mindiff = -50;
    RecS.maxdiff = 50;
    RecS.Unit = "mg m-3";
  }
  if (FullVarName == "ChlorophyllConcOCX") {
    RecS.VarName2 = "Chlorophyll Concentration OCX Algorithm";
    RecS.CFshortName = "chl_ocx";
    RecS.minval = 0;
    RecS.maxval = 100;
    RecS.mindiff = -50;
    RecS.maxdiff = 50;
    RecS.Unit = "mg m-3";
  }
  if (FullVarName == "CalciteConc") {
    RecS.VarName2 = "Calcite Concentration";
    RecS.CFshortName = "pic";
    RecS.minval = 0;
    RecS.maxval = 0.125;
    RecS.mindiff = -0.03;
    RecS.maxdiff = 0.03;
    RecS.Unit = "mol m-3";
  }
  if (FullVarName == "ParticleOrganicCarbon") {
    RecS.VarName2 = "Particle Organic Carbon";
    RecS.CFshortName = "poc";
    RecS.minval = 0;
    RecS.maxval = 1000;
    RecS.mindiff = -300;
    RecS.maxdiff = 300;
    RecS.Unit = "mg m-3";
  }
  if (FullVarName == "InstPhotosyntheticallyAvailableRad") {
    RecS.VarName2 = "Instantaneous Photosynthetically Available Radiation";
    RecS.CFshortName = "ipar";
    RecS.minval = 0;
    RecS.maxval = 0.0032;
    RecS.mindiff = -0.001;
    RecS.maxdiff = 0.001;
    RecS.Unit = "einstein m-2 s-1";
  }
  if (FullVarName == "PhotosyntheticallyAvailableRad") {
    RecS.VarName2 = "Photosynthetically Available Radiation";
    RecS.CFshortName = "par";
    RecS.minval = 0;
    RecS.maxval = 130;
    RecS.mindiff = -50;
    RecS.maxdiff = 50;
    RecS.Unit = "einstein m-2 day-1";
  }
  if (FullVarName == "ChlorophylA") {
    if (eModelName == "CFCONVENTION")
      F = Get2DvariableSpecTime(TotalArr, "FieldOut1", eTimeDay);
    RecS.VarName2 = "Chlorophyl A concentration";
    RecS.minval = 0;
    RecS.maxval = 300;
    RecS.mindiff = 0;
    RecS.maxdiff = 300;
    RecS.Unit = "mg m-3";
  }
  if (FullVarName == "IOBP") {
    if (eModelName == "WWM") {
      // we should have a NETCDF_WW3
      F = Get2DvariableSpecTime(TotalArr, "IOBP_WW3", eTimeDay);
    }
    if (eModelName == "UNRUNOFF") {
      int mnp = TotalArr.GrdArr.IOBP.size();
      if (mnp == 0) {
        std::cerr << "The IOBP has not been assigned originally. So cannot "
                     "plot IOBP\n";
        std::cerr << "Most likely the BoundFile is set to \"unset\" and that "
                     "is wrong\n";
        throw TerminalException{1};
      }
      MyVector<double> Fret(mnp);
      for (int i = 0; i < mnp; i++)
        Fret(i) = static_cast<double>(TotalArr.GrdArr.IOBP(i));
      F = Fret;
    }
    RecS.VarName2 = "IOBP of the unstructured model";
    RecS.minval = 0;
    RecS.maxval = 1;
    RecS.mindiff = -1;
    RecS.maxdiff = 1;
    RecS.Unit = "nondim.";
  }
  if (FullVarName == "RoughnessFactor") {
    if (eModelName == "WWM") {
      F = GetRoughnessFactor(GetDEP(TotalArr.GrdArr.GrdArrRho),
                             TotalArr.GrdArr.INE);
    }
    RecS.VarName2 = "Roughness factor";
    RecS.minval = 0;
    RecS.maxval = 0.2;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.Unit = "nondim.";
  }
  if (FullVarName == "MAPSTA") {
    if (eModelName == "WWM")
      F = Get2DvariableSpecTime(TotalArr, "MAPSTA", eTimeDay);
    RecS.VarName2 = "MAPSTA of wavewatchIII";
    RecS.minval = -2;
    RecS.maxval = 2;
    RecS.mindiff = -1;
    RecS.maxdiff = 1;
    RecS.Unit = "nondim.";
  }
  //
  // Atmospheric variables
  //
  if (FullVarName == "Uwind") {
    RecVar RecVarWork =
        ModelSpecificVarSpecificTime_Kernel(TotalArr, "WIND10", eTimeDay);
    F = RecVarWork.U;
    RecS.VarName2 = "Eastward wind";
    RecS.minval = -10;
    RecS.maxval = 10;
    RecS.mindiff = -2;
    RecS.maxdiff = 2;
    RecS.Unit = "m/s";
  }
  if (FullVarName == "Vwind") {
    RecVar RecVarWork =
        ModelSpecificVarSpecificTime_Kernel(TotalArr, "WIND10", eTimeDay);
    F = RecVarWork.V;
    RecS.VarName2 = "Northward wind";
    RecS.minval = -10;
    RecS.maxval = 10;
    RecS.mindiff = -2;
    RecS.maxdiff = 2;
    RecS.Unit = "m/s";
  }
  if (FullVarName == "WIND10") {
    if (eModelName == "UNRUNOFF") {
      U = Get2DvariableSpecTime(TotalArr, "Uwind", eTimeDay);
      V = Get2DvariableSpecTime(TotalArr, "Vwind", eTimeDay);
    }
    if (eModelName == "SCHISM_SFLUX") {
      U = Get2DvariableSpecTime(TotalArr, "uwind", eTimeDay);
      V = Get2DvariableSpecTime(TotalArr, "vwind", eTimeDay);
    }
    if (eModelName == "SCHISM_NETCDF_OUT") {
      U = Get2DvariableSpecTime(TotalArr, "Uwind", eTimeDay);
      V = Get2DvariableSpecTime(TotalArr, "Vwind", eTimeDay);
    }
    if (eModelName == "ROMS" || eModelName == "WWM") {
      U = Get2DvariableSpecTime(TotalArr, "Uwind", eTimeDay);
      V = Get2DvariableSpecTime(TotalArr, "Vwind", eTimeDay);
    }
    if (eModelName == "COSMO" || eModelName == "WAM") {
      U = Get2DvariableSpecTime(TotalArr, "U_10", eTimeDay);
      V = Get2DvariableSpecTime(TotalArr, "V_10", eTimeDay);
    }
    if (eModelName == "WRF") {
      U = Get2DvariableSpecTime(TotalArr, "U10", eTimeDay);
      V = Get2DvariableSpecTime(TotalArr, "V10", eTimeDay);
    }
    if (eModelName == "WW3") {
      U = Get2DvariableSpecTime(TotalArr, "uwnd", eTimeDay);
      V = Get2DvariableSpecTime(TotalArr, "vwnd", eTimeDay);
      /*
        That part seems old code, but you never know with WW3.
      Eigen::Tensor<double, 3> Utens =
          NETCDF_Get3DvariableSpecTime(TotalArr, "u10m", eTimeDay);
      Eigen::Tensor<double, 3> Vtens =
          NETCDF_Get3DvariableSpecTime(TotalArr, "v10m", eTimeDay);
      U = DimensionExtraction(Utens, 0, 0);
      V = DimensionExtraction(Vtens, 0, 0);
      */
    }
    if (eModelName == "GRIB_DWD" || eModelName == "GRIB_GFS" ||
        eModelName == "GRIB_ECMWF" || eModelName == "GRIB_COSMO" ||
        eModelName == "GRIB_ALADIN" || eModelName == "GRIB_IFS") {
      U = Get2DvariableSpecTime(TotalArr, "10u", eTimeDay);
      V = Get2DvariableSpecTime(TotalArr, "10v", eTimeDay);
    }
    AngleRhoRot(U, V, TotalArr.GrdArr.GrdArrRho.ANG);
    RecS.VarName2 = "10m wind";
    RecS.minval = 0;
    RecS.maxval = 13;
    RecS.mindiff = -2;
    RecS.maxdiff = 2;
    RecS.Unit = "m/s";
    RecS.VarNature = "uv";
    RecS.nameU = "Uwind";
    RecS.nameV = "Vwind";
    RecS.strTime_ROMS = "wind_time";
    RecS.varName_ROMS_U = "Uwind";
    RecS.varName_ROMS_V = "Vwind";
  }
  if (FullVarName == "aermssomhphil") {
    if (eModelName == "GRIB_ECMWF")
      F = Get2DvariableSpecTime(TotalArr, "aermssomhphil", eTimeDay);
    RecS.VarName2 = "Vert. int. hydrophilic organic matter aerosol";
    RecS.minval = 0;
    RecS.maxval = 0.0001;
    RecS.mindiff = -0.00001;
    RecS.maxdiff = 0.00001;
    RecS.Unit = "kg m^2";
  }
  if (FullVarName == "aermssomhphob") {
    if (eModelName == "GRIB_ECMWF")
      F = Get2DvariableSpecTime(TotalArr, "aermssomhphob", eTimeDay);
    RecS.VarName2 = "Vert. int. hydrophobic organic matter aerosol";
    RecS.minval = 0;
    RecS.maxval = 0.0001;
    RecS.mindiff = -0.00001;
    RecS.maxdiff = 0.00001;
    RecS.Unit = "kg m^2";
  }
  if (FullVarName == "aermssdus") {
    if (eModelName == "GRIB_ECMWF")
      F = Get2DvariableSpecTime(TotalArr, "aermssdus", eTimeDay);
    RecS.VarName2 = "Vert. int. aerosol dust (0.03 - 0.55)";
    RecS.minval = 0;
    RecS.maxval = 0.0001;
    RecS.mindiff = -0.00001;
    RecS.maxdiff = 0.00001;
    RecS.Unit = "kg m^2";
  }
  if (FullVarName == "aermssdum") {
    if (eModelName == "GRIB_ECMWF")
      F = Get2DvariableSpecTime(TotalArr, "aermssdum", eTimeDay);
    RecS.VarName2 = "Vert. int. aerosol dust (0.55 - 9)";
    RecS.minval = 0;
    RecS.maxval = 0.0001;
    RecS.mindiff = -0.00001;
    RecS.maxdiff = 0.00001;
    RecS.Unit = "kg m^2";
  }
  if (FullVarName == "aermssdul") {
    if (eModelName == "GRIB_ECMWF")
      F = Get2DvariableSpecTime(TotalArr, "aermssdul", eTimeDay);
    RecS.VarName2 = "Vert. int. aerosol dust (9 - 20)";
    RecS.minval = 0;
    RecS.maxval = 0.0001;
    RecS.mindiff = -0.00001;
    RecS.maxdiff = 0.00001;
    RecS.Unit = "kg m^2";
  }
  if (FullVarName == "aermssduSML") {
    if (eModelName == "GRIB_ECMWF") {
      F = Get2DvariableSpecTime(TotalArr, "aermssdus", eTimeDay);
      F += Get2DvariableSpecTime(TotalArr, "aermssdum", eTimeDay);
      F += Get2DvariableSpecTime(TotalArr, "aermssdul", eTimeDay);
    }
    RecS.VarName2 = "Vert. int. aerosol dust";
    RecS.minval = 0;
    RecS.maxval = 0.0001;
    RecS.mindiff = -0.00001;
    RecS.maxdiff = 0.00001;
    RecS.Unit = "kg m^2";
  }
  if (FullVarName == "SourceGainDustAerosolSmall") {
    if (eModelName == "GRIB_ECMWF")
      F = Get2DvariableSpecTime(TotalArr, "aersrcdus", eTimeDay);
    RecS.VarName2 = "Source/gain dust aerosol (0.03 - 0.55)";
    RecS.minval = 0;
    RecS.maxval = 13;
    RecS.mindiff = -2;
    RecS.maxdiff = 2;
    RecS.Unit = "kg m^2 s^-1";
  }
  if (FullVarName == "SourceGainDustAerosolMedium") {
    if (eModelName == "GRIB_ECMWF")
      F = Get2DvariableSpecTime(TotalArr, "aersrcdum", eTimeDay);
    RecS.VarName2 = "Source/gain dust aerosol (0.55 - 9)";
    RecS.minval = 0;
    RecS.maxval = 13;
    RecS.mindiff = -2;
    RecS.maxdiff = 2;
    RecS.Unit = "kg m^2 s^-1";
  }
  if (FullVarName == "SourceGainDustAerosolLarge") {
    if (eModelName == "GRIB_ECMWF")
      F = Get2DvariableSpecTime(TotalArr, "aersrcdul", eTimeDay);
    RecS.VarName2 = "Source/gain dust aerosol (9 - 20)";
    RecS.minval = 0;
    RecS.maxval = 13;
    RecS.mindiff = -2;
    RecS.maxdiff = 2;
    RecS.Unit = "kg m^2 s^-1";
  }
  if (FullVarName == "SourceGainDustAerosolSML") {
    if (eModelName == "GRIB_ECMWF") {
      F = Get2DvariableSpecTime(TotalArr, "aersrcdus", eTimeDay);
      F += Get2DvariableSpecTime(TotalArr, "aersrcdum", eTimeDay);
      F += Get2DvariableSpecTime(TotalArr, "aersrcdul", eTimeDay);
    }
    if (eModelName == "GEOS") {
      F = Get2DvariableSpecTime(TotalArr, "DUDP001", eTimeDay);
      F += Get2DvariableSpecTime(TotalArr, "DUDP002", eTimeDay);
      F += Get2DvariableSpecTime(TotalArr, "DUDP003", eTimeDay);
      F += Get2DvariableSpecTime(TotalArr, "DUDP004", eTimeDay);
      F += Get2DvariableSpecTime(TotalArr, "DUDP005", eTimeDay);
      F *= 86400;
    }
    RecS.VarName2 = "Total dust dry deposition";
    RecS.minval = 0;
    RecS.maxval = 13;
    RecS.mindiff = -2;
    RecS.maxdiff = 2;
    RecS.Unit = "kg m^2 day^-1";
  }
  if (FullVarName == "WINDMAG") {
    if (eModelName == "ROMS" || eModelName == "WWM") {
      if (TOTALARR_IsVar(TotalArr, "Uwind") &&
          TOTALARR_IsVar(TotalArr, "Vwind")) {
        MyMatrix<double> Us =
            Get2DvariableSpecTime(TotalArr, "Uwind", eTimeDay);
        MyMatrix<double> Vs =
            Get2DvariableSpecTime(TotalArr, "Vwind", eTimeDay);
        F = COMPUTE_NORM(Us, Vs);
      } else {
        if (eModelName == "WWM")
          F = Get2DvariableSpecTime(TotalArr, "WINDMAG", eTimeDay);
        else
          F = Get2DvariableSpecTime(TotalArr, "WNDMAG", eTimeDay);
      }
    }
    if (eModelName == "WW3") {
      MyMatrix<double> Us = Get2DvariableSpecTime(TotalArr, "uwnd", eTimeDay);
      MyMatrix<double> Vs = Get2DvariableSpecTime(TotalArr, "vwnd", eTimeDay);
      F = COMPUTE_NORM(Us, Vs);
    }
    if (eModelName == "UNRUNOFF") {
      MyMatrix<double> Us = Get2DvariableSpecTime(TotalArr, "Uwind", eTimeDay);
      MyMatrix<double> Vs = Get2DvariableSpecTime(TotalArr, "Vwind", eTimeDay);
      F = COMPUTE_NORM(Us, Vs);
    }
    if (eModelName == "SCHISM_SFLUX") {
      MyMatrix<double> Us = Get2DvariableSpecTime(TotalArr, "uwind", eTimeDay);
      MyMatrix<double> Vs = Get2DvariableSpecTime(TotalArr, "vwind", eTimeDay);
      F = COMPUTE_NORM(Us, Vs);
    }
    if (eModelName == "SCHISM_NETCDF_OUT") {
      MyMatrix<double> Us = Get2DvariableSpecTime(TotalArr, "Uwind", eTimeDay);
      MyMatrix<double> Vs = Get2DvariableSpecTime(TotalArr, "Vwind", eTimeDay);
      F = COMPUTE_NORM(Us, Vs);
    }
    if (eModelName == "COSMO" || eModelName == "WAM") {
      MyMatrix<double> Us = Get2DvariableSpecTime(TotalArr, "U_10", eTimeDay);
      MyMatrix<double> Vs = Get2DvariableSpecTime(TotalArr, "V_10", eTimeDay);
      F = COMPUTE_NORM(Us, Vs);
    }
    if (eModelName == "GRIB_DWD" || eModelName == "GRIB_GFS" ||
        eModelName == "GRIB_ECMWF" || eModelName == "GRIB_COSMO" ||
        eModelName == "GRIB_ALADIN" || eModelName == "GRIB_IFS") {
      MyMatrix<double> Us = Get2DvariableSpecTime(TotalArr, "10u", eTimeDay);
      MyMatrix<double> Vs = Get2DvariableSpecTime(TotalArr, "10v", eTimeDay);
      F = COMPUTE_NORM(Us, Vs);
    }
    if (eModelName == "GRIB_WAM_FORT30")
      F = Get2DvariableSpecTime(TotalArr, "wind", eTimeDay);
    RecS.VarName2 = "10m wind speed";
    RecS.minval = 0;
    RecS.maxval = 13;
    RecS.mindiff = -2;
    RecS.maxdiff = 2;
    RecS.Unit = "m/s";
  }
  if (FullVarName == "AIRD") {
    if (eModelName == "COSMO" || eModelName == "WAM")
      F = Get2DvariableSpecTime(TotalArr, "AIRD", eTimeDay);
    RecS.VarName2 = "air density";
    RecS.minval = 1.12;
    RecS.maxval = 1.20;
    RecS.mindiff = -0.02;
    RecS.maxdiff = 0.02;
    RecS.Unit = "kg/m3";
  }
  if (FullVarName == "rain") {
    if (eModelName == "UNRUNOFF")
      F = Get2DvariableSpecTime(TotalArr, "rain", eTimeDay);
    if (eModelName == "SCHISM_SFLUX")
      F = Get2DvariableSpecTime(TotalArr, "prate", eTimeDay);
    if (eModelName == "ROMS")
      F = Get2DvariableSpecTime(TotalArr, "rain", eTimeDay);
    //    if (eModelName == "GRIB_DWD")
    //      F=Get2DvariableSpecTime(TotalArr, "tp", eTimeDay);
    if (eModelName == "WRF") {
      F = GRID_Get2DVariableTimeDifferentiate(TotalArr, "RAINNC", eTimeDay);
    }
    std::vector<std::string> ListGRIBmodel{"GRIB_COSMO", "GRIB_DWD",
                                           "GRIB_ECMWF", "GRIB_ALADIN"};
    if (PositionVect(ListGRIBmodel, eModelName) != -1) {
      std::cerr << "Calling GRID_Get2DVariableTimeDifferentiate for tp\n";
      // Conversion from m/s to kg/m^2/s
      double fact = 0.001;
      F = fact * GRID_Get2DVariableTimeDifferentiate(TotalArr, "tp", eTimeDay);
    }
    RemoveNegativeValues(F);
    RecS.VarName2 = "rainfall rate";
    RecS.minval = 0;
    RecS.maxval = 0.001;
    RecS.mindiff = -0.001;
    RecS.maxdiff = 0.001;
    RecS.Unit = "kg/m^2/s";
    RecS.strTime_ROMS = "rain_time";
    RecS.varName_ROMS = "rain";
  }
  if (FullVarName == "swrad") {
    if (eModelName == "SCHISM_SFLUX")
      F = Get2DvariableSpecTime(TotalArr, "dswrf", eTimeDay);
    if (eModelName == "WRF")
      F = Get2DvariableSpecTime(TotalArr, "SWDOWN", eTimeDay);
    if (eModelName == "ROMS")
      F = Get2DvariableSpecTime(TotalArr, "swrad", eTimeDay);
    if (eModelName == "SCHISM_NETCDF_OUT")
      F = Get2DvariableSpecTime(TotalArr, "srad", eTimeDay);
    if (eModelName == "GRIB_COSMO")
      F = Get2DvariableSpecTime(TotalArr, "sobs_rad", eTimeDay);
    if (eModelName == "GRIB_ALADIN") {
      F = GRID_Get2DVariableTimeDifferentiate(TotalArr, "nswrs", eTimeDay);
      RemoveNegativeValues(F);
      SetNegativeDeepNight(F, eTimeDay);
    }
    if (eModelName == "GRIB_ECMWF") {
      std::cerr << "Passing here for swrad\n";
      F = GRID_Get2DVariableTimeDifferentiate(TotalArr, "ssrd", eTimeDay);
      RemoveNegativeValues(F);
      //      SetNegativeDeepNight(F, eTimeDay);
      double smallVal = 0.1;
      SetSmallToZero(F, smallVal);
    }
    //      F=Get2DvariableSpecTime(TotalArr, "ssrd", eTimeDay);
    RecS.VarName2 = "Shortwave flux";
    RecS.minval = 100;
    RecS.maxval = 1000;
    RecS.mindiff = -100;
    RecS.maxdiff = 100;
    RecS.Unit = "W/m2";
    RecS.strTime_ROMS = "srf_time";
    RecS.varName_ROMS = "swrad";
  }
  if (FullVarName == "lwrad") {
    if (eModelName == "SCHISM_SFLUX")
      F = Get2DvariableSpecTime(TotalArr, "dlwrf", eTimeDay);
    if (eModelName == "SCHISM_NETCDF_OUT") {
      MyMatrix<double> F1 = Get2DvariableSpecTime(TotalArr, "hradu", eTimeDay);
      MyMatrix<double> F2 = Get2DvariableSpecTime(TotalArr, "hradd", eTimeDay);
      F = F1 - F2;
    }
    if (eModelName == "WRF")
      F = Get2DvariableSpecTime(TotalArr, "GLW", eTimeDay);
    if (eModelName == "ROMS")
      F = Get2DvariableSpecTime(TotalArr, "lwrad", eTimeDay);
    if (eModelName == "GRIB_COSMO")
      F = Get2DvariableSpecTime(TotalArr, "thbs_rad", eTimeDay);
    if (eModelName == "GRIB_ECMWF")
      F = GRID_Get2DVariableTimeDifferentiate(TotalArr, "strd", eTimeDay);
    //      F=Get2DvariableSpecTime(TotalArr, "thbs_rad", eTimeDay);
    RecS.VarName2 = "Longwave flux";
    RecS.minval = 200;
    RecS.maxval = 500;
    RecS.mindiff = -50;
    RecS.maxdiff = 50;
    RecS.Unit = "W/m2";
    RecS.strTime_ROMS = "lrf_time";
    RecS.varName_ROMS = "lwrad_down";
  }
  if (FullVarName == "latent") {
    if (eModelName == "UNRUNOFF")
      F = Get2DvariableSpecTime(TotalArr, "lat_flux", eTimeDay);
    if (eModelName == "ROMS")
      F = Get2DvariableSpecTime(TotalArr, "latent", eTimeDay);
    if (eModelName == "SCHISM_NETCDF_OUT")
      F = Get2DvariableSpecTime(TotalArr, "fluxlu", eTimeDay);
    RecS.VarName2 = "Latent flux";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.Unit = "W/m2";
  }
  if (FullVarName == "SurfPres") {
    //    std::cerr << "eModelName=" << eModelName << "\n";
    if (eModelName == "ROMS") {
      MyMatrix<double> Fin = Get2DvariableSpecTime(TotalArr, "Pair", eTimeDay);
      F = 100 * Fin;
    }
    if (eModelName == "UNRUNOFF")
      F = Get2DvariableSpecTime(TotalArr, "Pair", eTimeDay);
    if (eModelName == "WRF")
      F = Get2DvariableSpecTime(TotalArr, "PSFC", eTimeDay);
    if (eModelName == "SCHISM_SFLUX")
      F = Get2DvariableSpecTime(TotalArr, "prmsl", eTimeDay);
    if (eModelName == "SCHISM_NETCDF_OUT")
      F = Get2DvariableSpecTime(TotalArr, "pr", eTimeDay);
    if (eModelName == "GRIB_DWD" || eModelName == "GRIB_GFS")
      F = Get2DvariableSpecTime(TotalArr, "prmsl", eTimeDay);
    if (eModelName == "GRIB_ECMWF" || eModelName == "GRIB_ALADIN") {
      //      std::cerr << "Retrieving the msl from ECMWF\n";
      F = Get2DvariableSpecTime(TotalArr, "msl", eTimeDay);
    }
    if (eModelName == "GRIB_COSMO") {
      if (TOTALARR_IsVar(TotalArr, "pmsl"))
        F = Get2DvariableSpecTime(TotalArr, "pmsl", eTimeDay);
      if (TOTALARR_IsVar(TotalArr, "msl")) {
        std::cerr << "Retrieving MSL\n";
        F = Get2DvariableSpecTime(TotalArr, "msl", eTimeDay);
      }
    }
    RecS.VarName2 = "mean sea level pressure";
    RecS.minval = 100000;
    RecS.maxval = 103000;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.Unit = "Pa";
    RecS.strTime_ROMS = "pair_time";
    RecS.varName_ROMS = "Pair";
  }
  if (FullVarName == "sensible") {
    if (eModelName == "UNRUNOFF")
      F = Get2DvariableSpecTime(TotalArr, "sen_flux", eTimeDay);
    if (eModelName == "ROMS")
      F = Get2DvariableSpecTime(TotalArr, "sensible", eTimeDay);
    if (eModelName == "SCHISM_NETCDF_OUT")
      F = Get2DvariableSpecTime(TotalArr, "fluxsu", eTimeDay);
    RecS.VarName2 = "Sensible heat flux";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.Unit = "W/m2";
  }
  if (FullVarName == "shflux") {
    if (eModelName == "ROMS")
      F = Get2DvariableSpecTime(TotalArr, "shflux", eTimeDay);
    if (eModelName == "GRIB_ECMWF") {
      // follows
      // https://github.com/dcherian/tools/blob/master/ROMS/arango/forcing/d_ecmwf2roms.m
      MyMatrix<double> sensbl =
          Get2DvariableSpecTime(TotalArr, "sshf", eTimeDay);
      MyMatrix<double> latent =
          Get2DvariableSpecTime(TotalArr, "slhf", eTimeDay);
      MyMatrix<double> nlwrad =
          Get2DvariableSpecTime(TotalArr, "str", eTimeDay);
      MyMatrix<double> nswrad =
          Get2DvariableSpecTime(TotalArr, "ssr", eTimeDay);
      double scale = 1 / (3 * static_cast<double>(3600));
      F = (sensbl + latent + nlwrad + nswrad) * scale;
    }
    RecS.VarName2 = "Surface heat flux";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.Unit = "W/m2";
  }
  if (FullVarName == "ssflux") {
    if (eModelName == "ROMS")
      F = Get2DvariableSpecTime(TotalArr, "ssflux", eTimeDay);
    RecS.VarName2 = "Surface salinity flux";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.Unit = "PSU/m2s";
  }
  if (FullVarName == "evaporation") {
    if (eModelName == "ROMS")
      F = Get2DvariableSpecTime(TotalArr, "evaporation", eTimeDay);
    RecS.VarName2 = "Evaporation rate";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.Unit = "kg/m2s";
  }
  if (FullVarName == "Vorticity") {
    if (eModelName == "ROMS") {
      Eigen::Tensor<double, 3> Utot =
          NETCDF_Get3DvariableSpecTime(TotalArr, "u", eTimeDay);
      Eigen::Tensor<double, 3> Vtot =
          NETCDF_Get3DvariableSpecTime(TotalArr, "v", eTimeDay);
      MyMatrix<double> Usurf = DimensionExtraction(Utot, 0, 0);
      MyMatrix<double> Vsurf = DimensionExtraction(Vtot, 0, 0);
      F = GRID_VorticityRho(TotalArr.GrdArr, Usurf, Vsurf);
    }
    RecS.VarName2 = "Vorticity";
    RecS.minval = 0;
    RecS.maxval = 100;
    RecS.mindiff = -20;
    RecS.maxdiff = 20;
    RecS.Unit = "1/s";
  }
  if (FullVarName == "CloudFraction") {
    if (eModelName == "ROMS")
      F = Get2DvariableSpecTime(TotalArr, "cloud", eTimeDay);
    if (eModelName == "GRIB_ALADIN" || eModelName == "GRIB_ECMWF")
      F = Get2DvariableSpecTime(TotalArr, "tcc", eTimeDay);
    RecS.VarName2 = "Cloud fraction";
    RecS.minval = 0;
    RecS.maxval = 100;
    RecS.mindiff = -20;
    RecS.maxdiff = 20;
    RecS.Unit = "nondimensional";
    RecS.strTime_ROMS = "cloud_time";
    RecS.varName_ROMS = "cloud";
  }
  if (FullVarName == "AIRT2K") {
    RecVar RecVarWork =
        ModelSpecificVarSpecificTime_Kernel(TotalArr, "AIRT2", eTimeDay);
    F = RecVarWork.F;
    int siz = F.size();
    for (int i = 0; i < siz; i++)
      F(i) += static_cast<double>(273.15);
    RecS.VarName2 = "2m air temperature (K)";
    RecS.minval = 273.15 + 10;
    RecS.maxval = 273.15 + 20;
    RecS.mindiff = -2;
    RecS.maxdiff = 2;
    RecS.Unit = "deg K";
  }
  if (FullVarName == "AIRT2") {
    if (eModelName == "SCHISM_SFLUX") {
      F = Get2DvariableSpecTime(TotalArr, "stmp", eTimeDay);
      int siz = F.size();
      for (int i = 0; i < siz; i++)
        F(i) -= static_cast<double>(273.15);
    }
    if (eModelName == "SCHISM_NETCDF_OUT") {
      F = Get2DvariableSpecTime(TotalArr, "airt1", eTimeDay);
      int siz = F.size();
      for (int i = 0; i < siz; i++)
        F(i) -= static_cast<double>(273.15);
    }
    if (eModelName == "ROMS") {
      F = Get2DvariableSpecTime(TotalArr, "Tair", eTimeDay);
    }
    if (eModelName == "COSMO") {
      F = Get2DvariableSpecTime(TotalArr, "t_2m", eTimeDay);
      int siz = F.size();
      for (int i = 0; i < siz; i++)
        F(i) -= static_cast<double>(273.15);
    }
    if (eModelName == "WRF") {
      F = Get2DvariableSpecTime(TotalArr, "T2", eTimeDay);
      int siz = F.size();
      for (int i = 0; i < siz; i++)
        F(i) -= static_cast<double>(273.15);
    }
    std::vector<std::string> ListModel{"GRIB_DWD", "GRIB_ECMWF", "GRIB_GFS",
                                       "GRIB_COSMO", "GRIB_ALADIN"};
    if (PositionVect(ListModel, eModelName) != -1) {
      F = Get2DvariableSpecTime(TotalArr, "2t", eTimeDay);
      int siz = F.size();
      for (int i = 0; i < siz; i++)
        F(i) -= static_cast<double>(273.15);
    }
    RecS.VarName2 = "2m air temperature";
    RecS.minval = 10;
    RecS.maxval = 20;
    RecS.mindiff = -2;
    RecS.maxdiff = 2;
    RecS.Unit = "deg C";
    RecS.strTime_ROMS = "tair_time";
    RecS.varName_ROMS = "Tair";
  }
  if (FullVarName == "Rh2frac") {
    RecVar RecVarWork =
        ModelSpecificVarSpecificTime_Kernel(TotalArr, "Rh2", eTimeDay);
    F = RecVarWork.F / static_cast<double>(100);
    RecS.VarName2 = "2m relative humidity";
    RecS.minval = 0;
    RecS.maxval = 1;
    RecS.mindiff = -0.2;
    RecS.maxdiff = 0.2;
    RecS.Unit = "nondim.";
  }
  if (FullVarName == "Rh2") {
    if (eModelName == "COSMO")
      F = Get2DvariableSpecTime(TotalArr, "rh_2m", eTimeDay);
    if (eModelName == "GRIB_ALADIN")
      F = Get2DvariableSpecTime(TotalArr, "r", eTimeDay);
    if (eModelName == "GRIB_DWD")
      F = Get2DvariableSpecTime(TotalArr, "RELHUM_2M", eTimeDay);
    if (eModelName == "SCHISM_NETCDF_OUT")
      F = Get2DvariableSpecTime(TotalArr, "shum1", eTimeDay);
    if (eModelName == "SCHISM_SFLUX") {
      MyMatrix<double> Fin = Get2DvariableSpecTime(TotalArr, "spfh", eTimeDay);
      F = 100 * Fin;
    }
    if (eModelName == "ROMS")
      F = Get2DvariableSpecTime(TotalArr, "Qair", eTimeDay);
    if (eModelName == "GRIB_ECMWF")
      F = Algorithms_RelativeHumidity(TotalArr, "2d", "2r", "2t", "msl", "q",
                                      eTimeDay);
    if (eModelName == "GRIB_COSMO")
      F = Algorithms_RelativeHumidity(TotalArr, "2d", "2r", "2t", "msl", "QV_S",
                                      eTimeDay);
    if (eModelName == "WRF") {
      MyMatrix<double> Q2 = Get2DvariableSpecTime(TotalArr, "Q2", eTimeDay);
      MyMatrix<double> PSFC = Get2DvariableSpecTime(TotalArr, "PSFC", eTimeDay);
      MyMatrix<double> T2K = Get2DvariableSpecTime(TotalArr, "T2", eTimeDay);
      F = mixing_ratio2relative_humidity(Q2, PSFC, T2K);
    }
    RecS.VarName2 = "2m relative humidity";
    RecS.minval = 0;
    RecS.maxval = 100;
    RecS.mindiff = -20;
    RecS.maxdiff = 20;
    RecS.Unit = "nondim.";
    RecS.strTime_ROMS = "qair_time";
    RecS.varName_ROMS = "Qair";
  }
  //
  // Oceanic variables
  //
  if (FullVarName == "HorizCurr") {
    VerticalLevelInfo VertInfo;
    if (eModelName != "TRIVIAL")
      VertInfo = RetrieveVerticalInformation(FullVarName, eModelName);
    if (eModelName == "WWM") {
      std::string strCallU =
          "UhorizCurr_" + VertInfo.type + VertInfo.strDepth + "m";
      std::string strCallV =
          "VhorizCurr_" + VertInfo.type + VertInfo.strDepth + "m";
      U = Get2DvariableSpecTime(TotalArr, strCallU, eTimeDay);
      V = Get2DvariableSpecTime(TotalArr, strCallV, eTimeDay);
    }
    RecS.VarName2 = "horizontal current" + VertInfo.strDepth;
    RecS.minval = 0;
    RecS.maxval = 0.5;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.Unit = "m/s";
    RecS.VarNature = "uv";
    RecS.nameU = "UhorizCurr";
    RecS.nameV = "VhorizCurr";
  }
  if (FullVarName == "Curr") {
    if (eModelName == "ROMS") {
      Eigen::Tensor<double, 3> Utot =
          NETCDF_Get3DvariableSpecTime(TotalArr, "u", eTimeDay);
      Eigen::Tensor<double, 3> Vtot =
          NETCDF_Get3DvariableSpecTime(TotalArr, "v", eTimeDay);
      Uthree = My_u2rho_3D(Utot, TotalArr.GrdArr.GrdArrU.MSK);
      Vthree = My_v2rho_3D(Vtot, TotalArr.GrdArr.GrdArrV.MSK);
    }
    if (eModelName == "AREG") {
      Eigen::Tensor<double, 3> Utot =
          NETCDF_Get3DvariableSpecTime(TotalArr, "U", eTimeDay);
      Eigen::Tensor<double, 3> Vtot =
          NETCDF_Get3DvariableSpecTime(TotalArr, "V", eTimeDay);
      Uthree = My_u2rho_3D(Utot, TotalArr.GrdArr.GrdArrU.MSK);
      Vthree = My_v2rho_3D(Vtot, TotalArr.GrdArr.GrdArrV.MSK);
    }
    if (eModelName == "WWM") {
      Uthree = NETCDF_Get3DvariableSpecTime(TotalArr, "Ucurr", eTimeDay);
      Vthree = NETCDF_Get3DvariableSpecTime(TotalArr, "Vcurr", eTimeDay);
    }
    if (eModelName == "NEMO") {
      Uthree = NETCDF_Get3DvariableSpecTime(TotalArr, "uo", eTimeDay);
      Vthree = NETCDF_Get3DvariableSpecTime(TotalArr, "uo", eTimeDay);
    }
    if (eModelName == "HYCOM") {
      Uthree = NETCDF_Get3DvariableSpecTime(TotalArr, "water_u", eTimeDay);
      Vthree = NETCDF_Get3DvariableSpecTime(TotalArr, "water_v", eTimeDay);
    }
    AngleRhoRot_3D(Uthree, Vthree, TotalArr.GrdArr.GrdArrRho.ANG);
    RecS.VarName2 = "baroclinic current";
    RecS.minval = -0.2; // We need a negative value for plotting u or v.
    RecS.maxval = 0.2;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Duv";
    RecS.Unit = "m/s";
  }
  if (FullVarName == "Temp") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "temp", eTimeDay);
    if (eModelName == "AREG")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "T", eTimeDay);
    if (eModelName == "NEMO")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "thetao", eTimeDay);
    if (eModelName == "HYCOM")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "water_temp", eTimeDay);
    RecS.VarName2 = "temperature";
    RecS.minval = 17;
    RecS.maxval = 25;
    RecS.mindiff = -2;
    RecS.maxdiff = 2;
    RecS.VarNature = "3Drho";
    RecS.Unit = "deg C";
  }
  if (FullVarName == "Salt") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "salt", eTimeDay);
    if (eModelName == "AREG")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "S", eTimeDay);
    if (eModelName == "NEMO")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "so", eTimeDay);
    if (eModelName == "HYCOM")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "salinity", eTimeDay);
    RecS.VarName2 = "salinity";
    RecS.minval = 35;
    RecS.maxval = 38;
    RecS.mindiff = -1;
    RecS.maxdiff = 1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "PSU";
  }
  if (FullVarName == "CurrMag") {
    RecVar RecVarWork =
        ModelSpecificVarSpecificTime_Kernel(TotalArr, "Curr", eTimeDay);
    Tens3 = ComputeNormPairOfTensor(RecVarWork.Uthree, RecVarWork.Vthree);
    RecS.VarName2 = "baroclinic current magnitude";
    RecS.minval = 0;
    RecS.maxval = 0.2;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "m/s";
  }
  if (FullVarName == "CurrBaro") {
    if (eModelName == "UNRUNOFF") {
      MyMatrix<double> Ucurr =
          Get2DvariableSpecTime(TotalArr, "CURTX", eTimeDay);
      MyMatrix<double> Vcurr =
          Get2DvariableSpecTime(TotalArr, "CURTY", eTimeDay);
      MyMatrix<double> Helev = Get2DvariableSpecTime(TotalArr, "H", eTimeDay);
      U = Helev.cwiseProduct(Ucurr);
      V = Helev.cwiseProduct(Vcurr);
    }
    if (eModelName == "ROMS") {
      MyMatrix<double> UBAR = Get2DvariableSpecTime(TotalArr, "ubar", eTimeDay);
      MyMatrix<double> VBAR = Get2DvariableSpecTime(TotalArr, "vbar", eTimeDay);
      U = My_u2rho(UBAR, TotalArr.GrdArr.GrdArrU.MSK);
      V = My_v2rho(VBAR, TotalArr.GrdArr.GrdArrV.MSK);
    }
    AngleRhoRot(U, V, TotalArr.GrdArr.GrdArrRho.ANG);
    RecS.VarName2 = "barotropic current";
    RecS.minval = 0;
    RecS.maxval = 0.5;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.Unit = "m^2/s";
    RecS.VarNature = "uv";
    RecS.nameU = "UsurfCurr";
    RecS.nameV = "VsurfCurr";
  }
  if (FullVarName == "ROMSbarotropicdefect") {
    if (eModelName == "ROMS") {
      MyMatrix<double> UBAR_mod =
          Get2DvariableSpecTime(TotalArr, "ubar", eTimeDay);
      MyMatrix<double> VBAR_mod =
          Get2DvariableSpecTime(TotalArr, "vbar", eTimeDay);
      Eigen::Tensor<double, 3> Utot =
          NETCDF_Get3DvariableSpecTime(TotalArr, "u", eTimeDay);
      Eigen::Tensor<double, 3> Vtot =
          NETCDF_Get3DvariableSpecTime(TotalArr, "v", eTimeDay);
      MyMatrix<double> zeta_rho =
          Get2DvariableSpecTime(TotalArr, "zeta", eTimeDay);
      MyMatrix<double> zeta_u = My_rho2u_2D(TotalArr.GrdArr, zeta_rho);
      MyMatrix<double> zeta_v = My_rho2v_2D(TotalArr.GrdArr, zeta_rho);
      MyMatrix<double> UBAR_int = ConvertBaroclinic_to_Barotropic_ARVD_Coord(
          Utot, zeta_u, TotalArr.GrdArr.ARVD, TotalArr.GrdArr.GrdArrU);
      MyMatrix<double> VBAR_int = ConvertBaroclinic_to_Barotropic_ARVD_Coord(
          Vtot, zeta_v, TotalArr.GrdArr.ARVD, TotalArr.GrdArr.GrdArrV);
      //
      MyMatrix<double> UBAR_pred = UBAR_mod - UBAR_int;
      MyMatrix<double> VBAR_pred = VBAR_mod - VBAR_int;
      MyMatrix<double> UBAR_diff = UBAR_pred.cwiseAbs();
      MyMatrix<double> VBAR_diff = VBAR_pred.cwiseAbs();
      //
      F = My_u2rho(UBAR_diff, TotalArr.GrdArr.GrdArrU.MSK) +
          My_v2rho(VBAR_diff, TotalArr.GrdArr.GrdArrV.MSK);
    }
    RecS.VarName2 = "ROMS barotropic defect";
    RecS.minval = 0;
    RecS.maxval = 0.5;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.Unit = "m^2/s";
  }
  if (FullVarName == "CurrBaroMag") {
    RecVar RecVarWork =
        ModelSpecificVarSpecificTime_Kernel(TotalArr, "CurrBaro", eTimeDay);
    F = COMPUTE_NORM(RecVarWork.U, RecVarWork.V);
    RecS.VarName2 = "barotropic current magnitude";
    RecS.minval = 0;
    RecS.maxval = 0.5;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.Unit = "m^2/s";
  }
  if (FullVarName == "SurfCurr") {
    if (eModelName == "ROMS") {
      Eigen::Tensor<double, 3> Utot =
          NETCDF_Get3DvariableSpecTime(TotalArr, "u", eTimeDay);
      Eigen::Tensor<double, 3> Vtot =
          NETCDF_Get3DvariableSpecTime(TotalArr, "v", eTimeDay);
      int s_rho = Utot.dimension(0);
      MyMatrix<double> Usurf = DimensionExtraction(Utot, 0, s_rho - 1);
      MyMatrix<double> Vsurf = DimensionExtraction(Vtot, 0, s_rho - 1);
      U = My_u2rho(Usurf, TotalArr.GrdArr.GrdArrU.MSK);
      V = My_v2rho(Vsurf, TotalArr.GrdArr.GrdArrV.MSK);
    }
    if (eModelName == "UNRUNOFF") {
      U = Get2DvariableSpecTime(TotalArr, "CURTX", eTimeDay);
      V = Get2DvariableSpecTime(TotalArr, "CURTY", eTimeDay);
    }
    if (eModelName == "WWM") {
      if (TOTALARR_IsVar(TotalArr, "CURTX") &&
          TOTALARR_IsVar(TotalArr, "CURTX")) {
        U = Get2DvariableSpecTime(TotalArr, "CURTX", eTimeDay);
        V = Get2DvariableSpecTime(TotalArr, "CURTY", eTimeDay);
      } else {
        if (TOTALARR_IsVar(TotalArr, "UsurfCurr") &&
            TOTALARR_IsVar(TotalArr, "VsurfCurr")) {
          U = Get2DvariableSpecTime(TotalArr, "UsurfCurr", eTimeDay);
          V = Get2DvariableSpecTime(TotalArr, "VsurfCurr", eTimeDay);
        } else {
          Eigen::Tensor<double, 3> Utot =
              NETCDF_Get3DvariableSpecTime(TotalArr, "Ucurr", eTimeDay);
          Eigen::Tensor<double, 3> Vtot =
              NETCDF_Get3DvariableSpecTime(TotalArr, "Vcurr", eTimeDay);
          int s_rho = Utot.dimension(0);
          U = DimensionExtraction(Utot, 0, s_rho - 1);
          V = DimensionExtraction(Vtot, 0, s_rho - 1);
        }
      }
    }
    if (eModelName == "COSMO" || eModelName == "WAM") {
      U = Get2DvariableSpecTime(TotalArr, "ucurr", eTimeDay);
      V = Get2DvariableSpecTime(TotalArr, "vcurr", eTimeDay);
    }
    if (eModelName == "WW3") {
      U = Get2DvariableSpecTime(TotalArr, "U", eTimeDay);
      V = Get2DvariableSpecTime(TotalArr, "V", eTimeDay);
    }
    if (eModelName == "HYCOM") {
      Uthree = NETCDF_Get3DvariableSpecTime(TotalArr, "water_u", eTimeDay);
      Vthree = NETCDF_Get3DvariableSpecTime(TotalArr, "water_v", eTimeDay);
      // surface is the first entry for HYCOM
      U = DimensionExtraction(Uthree, 0, 0);
      V = DimensionExtraction(Vthree, 0, 0);
    }
    AngleRhoRot(U, V, TotalArr.GrdArr.GrdArrRho.ANG);
    RecS.VarName2 = "surface current";
    RecS.minval = 0;
    RecS.maxval = 0.5;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.Unit = "m/s";
    RecS.VarNature = "uv";
    RecS.nameU = "UsurfCurr";
    RecS.nameV = "VsurfCurr";
  }
  if (FullVarName == "SurfStress") {
    if (eModelName == "UNRUNOFF") {
      U = Get2DvariableSpecTime(TotalArr, "sustr", eTimeDay);
      V = Get2DvariableSpecTime(TotalArr, "svstr", eTimeDay);
    }
    if (eModelName == "ROMS") {
      U = Get2DvariableSpecTime(TotalArr, "sustr", eTimeDay);
      V = Get2DvariableSpecTime(TotalArr, "svstr", eTimeDay);
    }
    AngleRhoRot(U, V, TotalArr.GrdArr.GrdArrRho.ANG);
    RecS.VarName2 = "surface stress";
    RecS.minval = 0;
    RecS.maxval = 0.5;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.Unit = "m/s";
    RecS.VarNature = "uv";
    RecS.nameU = "UsurfStress";
    RecS.nameV = "VsurfStress";
  }
  if (FullVarName == "SurfCurrMag") {
    RecVar RecVarWork =
        ModelSpecificVarSpecificTime_Kernel(TotalArr, "SurfCurr", eTimeDay);
    F = COMPUTE_NORM(RecVarWork.U, RecVarWork.V);
    RecS.VarName2 = "surface current magnitude";
    RecS.minval = 0;
    RecS.maxval = 0.5;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.Unit = "m/s";
  }
  if (FullVarName == "TempSurf") {
    if (eModelName == "COSMO") {
      F = Get2DvariableSpecTime(TotalArr, "t_s", eTimeDay);
      int siz = F.size();
      for (int i = 0; i < siz; i++)
        F(i) -= static_cast<double>(273.15);
    }
    if (eModelName == "SCHISM_NETCDF_OUT") {
      Eigen::Tensor<double, 3> TheTemp =
          NETCDF_Get3DvariableSpecTime(TotalArr, "tr_nd1", eTimeDay);
      int s_rho = TheTemp.dimension(0);
      F = DimensionExtraction(TheTemp, 0, s_rho - 1);
    }
    if (eModelName == "ROMS") {
      return ModelSpecificVarSpecificTime_Kernel(TotalArr, "Temp:SURF", eTimeDay);
    }
    RecS.VarName2 = "sea surface temperature";
    RecS.minval = 10;
    RecS.maxval = 20;
    RecS.mindiff = -2;
    RecS.maxdiff = 2;
    RecS.Unit = "deg C";
  }
  if (FullVarName == "DensAnomaly") {
    if (eModelName == "ROMS") {
      bool test = NETCDF_TestVariableAccessSpecTime(TotalArr, "rho", eTimeDay);
      if (test) {
        Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "rho", eTimeDay);
      } else {
        Eigen::Tensor<double, 3> TtArr =
          NETCDF_Get3DvariableSpecTime(TotalArr, "temp", eTimeDay);
        Eigen::Tensor<double, 3> TsArr =
          NETCDF_Get3DvariableSpecTime(TotalArr, "salt", eTimeDay);
        MyMatrix<double> zeta =
          NETCDF_Get2DvariableSpecTime(TotalArr, "zeta", eTimeDay);
        Tens3 = ROMS_ComputeDensityAnomaly(TsArr, TtArr, TotalArr.GrdArr, zeta);
      }
    }
    RecS.VarName2 = "density anomaly";
    RecS.minval = 30;
    RecS.maxval = 40;
    RecS.mindiff = -2;
    RecS.maxdiff = 2;
    RecS.VarNature = "3Drho";
    RecS.Unit = "kg m-3";
  }
  if (FullVarName == "Density") {
    RecVar RecVarWork =
      ModelSpecificVarSpecificTime_Kernel(TotalArr, "DensAnomaly", eTimeDay);
    Tens3 = GetDensity(RecVarWork.Tens3);
    RecS.VarName2 = "density";
    RecS.minval = 1030;
    RecS.maxval = 1040;
    RecS.mindiff = -2;
    RecS.maxdiff = 2;
    RecS.VarNature = "3Drho";
    RecS.Unit = "kg m-3";
  }
  if (FullVarName == "DynamicDepth") {
    if (eModelName != "TRIVIAL") {
      RecVar RecVarWork =
        ModelSpecificVarSpecificTime_Kernel(TotalArr, "Density:AVE-25m0m", eTimeDay);
      F = GetDynamicDepth(RecVarWork.F);
    }
    RecS.VarName2 = "dynamic depth";
    RecS.minval = 230;
    RecS.maxval = 240;
    RecS.mindiff = -10;
    RecS.maxdiff = 10;
    RecS.Unit = "m2 s-2";
  }
  if (FullVarName == "ZetaOcean") {
    if (eModelName == "COSMO")
      F = Get2DvariableSpecTime(TotalArr, "ZetaOcean", eTimeDay);
    if (eModelName == "ROMS")
      F = Get2DvariableSpecTime(TotalArr, "zeta", eTimeDay);
    if (eModelName == "AREG")
      F = Get2DvariableSpecTime(TotalArr, "SSH", eTimeDay);
    if (eModelName == "WWM")
      F = Get2DvariableSpecTime(TotalArr, "WATLEV", eTimeDay);
    if (eModelName == "WW3") {
      //      F = Get2DvariableSpecTime(TotalArr, "XE", eTimeDay);
      F = Get2DvariableSpecTime(TotalArr, "wlv", eTimeDay);
    }
    if (eModelName == "UNRUNOFF") {
      MyMatrix<double> TotalElev =
          Get2DvariableSpecTime(TotalArr, "H", eTimeDay);
      const MyMatrix<double> &DEP = GetDEP(TotalArr.GrdArr.GrdArrRho);
      if (!IsEqualSizeMatrices(TotalElev, DEP)) {
        std::cerr << "The matrices TotalElev and DEP have different sizes\n";
        std::cerr << "Most likely the grid does not match the history used\n";
        throw TerminalException{1};
      }
      F = TotalElev + DEP;
      //      std::cerr << "TotalElev min=" << TotalElev.minCoeff() << " / " <<
      //      TotalElev.maxCoeff() << "\n"; std::cerr << "      DEP min=" <<
      //      DEP.minCoeff() << " / " << DEP.maxCoeff() << "\n"; std::cerr << "
      //      F min=" << F.minCoeff() << " / " << F.maxCoeff() << "\n";
    }
    if (eModelName == "SCHISM_NETCDF_OUT")
      F = Get2DvariableSpecTime(TotalArr, "WATLEV", eTimeDay);
    if (eModelName == "NEMO")
      F = NETCDF_Get2DvariableSpecTime(TotalArr, "zos", eTimeDay);
    if (eModelName == "HYCOM")
      F = NETCDF_Get2DvariableSpecTime(TotalArr, "surf_el", eTimeDay);
    RecS.VarName2 = "free surface elevation";
    RecS.minval = -0.2;
    RecS.maxval = 0.2;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.Unit = "m";
  }
  if (FullVarName == "ZetaOceanDerivative") {
    if (eModelName == "WWM")
      F = Get2DvariableSpecTime(TotalArr, "DEPDT", eTimeDay);



    RecS.VarName2 = "free surface elevation derivative";
    RecS.minval = -0.01;
    RecS.maxval = 0.01;
    RecS.mindiff = -0.001;
    RecS.maxdiff = 0.001;
    RecS.Unit = "m/s";
  }
  //
  // Wave variables
  //
  if (FullVarName == "MeanWaveLength") {
    if (eModelName == "WWM")
      F = Get2DvariableSpecTime(TotalArr, "WLM", eTimeDay);
    RecS.VarName2 = "mean wave length";
    RecS.minval = 2;
    RecS.maxval = 30;
    RecS.mindiff = -5;
    RecS.maxdiff = 5;
    RecS.Unit = "m";
  }
  if (FullVarName == "PeakWaveLength") {
    if (eModelName == "WWM")
      F = Get2DvariableSpecTime(TotalArr, "LPP", eTimeDay);
    RecS.VarName2 = "peak wave length";
    RecS.minval = 2;
    RecS.maxval = 30;
    RecS.mindiff = -5;
    RecS.maxdiff = 5;
    RecS.Unit = "m";
  }
  if (FullVarName == "MeanWaveNumber") {
    if (eModelName == "WWM")
      F = Get2DvariableSpecTime(TotalArr, "KLM", eTimeDay);
    RecS.VarName2 = "mean wave number";
    RecS.minval = 0;
    RecS.maxval = 1;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.Unit = "m-1";
  }
  if (FullVarName == "PeakWaveNumber") {
    if (eModelName == "WWM")
      F = Get2DvariableSpecTime(TotalArr, "KPP", eTimeDay);
    RecS.VarName2 = "peak wave number";
    RecS.minval = 0;
    RecS.maxval = 1;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.Unit = "m-1";
  }
  if (FullVarName == "MeanWaveDir") {
    if (eModelName == "WWM")
      F = Get2DvariableSpecTime(TotalArr, "DM", eTimeDay);
    RecS.VarName2 = "mean wave direction";
    RecS.minval = 0;
    RecS.maxval = 360;
    RecS.mindiff = -30;
    RecS.maxdiff = 30;
    RecS.Unit = "deg";
  }
  if (FullVarName == "PeakWaveDir") {
    if (eModelName == "WWM")
      F = Get2DvariableSpecTime(TotalArr, "PEAKD", eTimeDay);
    RecS.VarName2 = "peak wave direction";
    RecS.minval = 0;
    RecS.maxval = 360;
    RecS.mindiff = -30;
    RecS.maxdiff = 30;
    RecS.Unit = "deg";
  }
  if (FullVarName == "MeanWaveDirVect") {
    if (eModelName == "WWM") {
      F = Get2DvariableSpecTime(TotalArr, "DM", eTimeDay);
      int nbRow = F.rows();
      int nbCol = F.cols();
      double deg2rad = 3.1415926535 / static_cast<double>(180);
      U = MyMatrix<double>(nbRow, nbCol);
      V = MyMatrix<double>(nbRow, nbCol);
      for (int iRow = 0; iRow < nbRow; iRow++)
        for (int iCol = 0; iCol < nbCol; iCol++) {
          double eAngRad = deg2rad * F(iRow, iCol);
          U(iRow, iCol) = cos(eAngRad);
          V(iRow, iCol) = sin(eAngRad);
        }
    }
    RecS.VarName2 = "mean wave direction";
    RecS.minval = 0;
    RecS.maxval = 360;
    RecS.mindiff = -30;
    RecS.maxdiff = 30;
    RecS.VarNature = "uv";
    RecS.Unit = "deg";
  }
  if (FullVarName == "PeakWaveDirVect") {
    if (eModelName == "WWM") {
      F = Get2DvariableSpecTime(TotalArr, "PEAKD", eTimeDay);
      int nbRow = F.rows();
      int nbCol = F.cols();
      double deg2rad = 3.1415926535 / static_cast<double>(180);
      U = MyMatrix<double>(nbRow, nbCol);
      V = MyMatrix<double>(nbRow, nbCol);
      for (int iRow = 0; iRow < nbRow; iRow++)
        for (int iCol = 0; iCol < nbCol; iCol++) {
          double eAngRad = deg2rad * F(iRow, iCol);
          U(iRow, iCol) = cos(eAngRad);
          V(iRow, iCol) = sin(eAngRad);
        }
    }
    RecS.VarName2 = "peak wave direction";
    RecS.minval = 0;
    RecS.maxval = 360;
    RecS.mindiff = -30;
    RecS.maxdiff = 30;
    RecS.VarNature = "uv";
    RecS.Unit = "deg";
  }
  if (FullVarName == "DiscPeakWaveDir") {
    if (eModelName == "WWM")
      F = Get2DvariableSpecTime(TotalArr, "DPEAK", eTimeDay);
    RecS.VarName2 = "discrete peak wave direction";
    RecS.minval = 0;
    RecS.maxval = 360;
    RecS.mindiff = -30;
    RecS.maxdiff = 30;
    RecS.Unit = "deg";
  }
  if (FullVarName == "ZetaSetup") {
    if (eModelName == "WWM")
      F = Get2DvariableSpecTime(TotalArr, "ZETA_SETUP", eTimeDay);
    if (eModelName == "WW3") {
      F = Get2DvariableSpecTime(TotalArr, "dpt", eTimeDay);
      F -= GetDEP(TotalArr.GrdArr.GrdArrRho);
    }
    RecS.VarName2 = "free surface setup";
    RecS.minval = 0;
    RecS.maxval = 0.76;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.Unit = "m";
  }
  if (FullVarName == "BreakingFraction") {
    MyMatrix<double> Fhs, Fzeta;
    auto f_set=[&]() -> void {
      const MyMatrix<double> &DEP = GetDEP(TotalArr.GrdArr.GrdArrRho);
      if (!IsEqualSizeMatrices(Fhs, DEP)) {
        std::cerr << "The matrices Fhs and DEP have different sizes\n";
        std::cerr << "Most likely the grid does not match the history used\n";
        throw TerminalException{1};
      }
      F = MyMatrix<double>(eta_rho, xi_rho);
      for (int i = 0; i < eta_rho; i++)
        for (int j = 0; j < xi_rho; j++)
          F(i, j) = Fhs(i, j) / (Fzeta(i, j) + DEP(i, j));
    };
    if (eModelName == "WWM") {
      Fhs = Get2DvariableSpecTime(TotalArr, "HS", eTimeDay);
      Fzeta = Get2DvariableSpecTime(TotalArr, "WATLEV", eTimeDay);
      f_set();
    }
    if (eModelName == "WW3") {
      Fhs = Get2DvariableSpecTime(TotalArr, "hs", eTimeDay);
      Fzeta = Get2DvariableSpecTime(TotalArr, "wlv", eTimeDay);
      f_set();
    }
    RecS.VarName2 = "Breaking fraction";
    RecS.minval = 0;
    RecS.maxval = 0.76;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.Unit = "nondim.";
  }
  if (FullVarName == "Hwave") {
    if (eModelName == "WWM")
      F = Get2DvariableSpecTime(TotalArr, "HS", eTimeDay);
    if (eModelName == "WW3")
      F = Get2DvariableSpecTime(TotalArr, "hs", eTimeDay);
    if (eModelName == "COSMO" || eModelName == "WAM")
      F = Get2DvariableSpecTime(TotalArr, "Hwave", eTimeDay);
    if (eModelName == "GRIB_WAM_FORT30" || eModelName == "GRIB_IFS")
      F = Get2DvariableSpecTime(TotalArr, "swh", eTimeDay);
    RecS.VarName2 = "Significant wave height";
    RecS.varName_GRIB = "swh";
    RecS.minval = 0;
    RecS.maxval = 4.5;
    RecS.mindiff = -0.5;
    RecS.maxdiff = 0.5;
    RecS.Unit = "m";
  }
  if (FullVarName == "MeanWaveFreq") {
    if (eModelName == "COSMO" || eModelName == "WAM")
      F = Get2DvariableSpecTime(TotalArr, "MwaveFreq", eTimeDay);
    if (eModelName == "WWM") {
      MyMatrix<double> Fin = Get2DvariableSpecTime(TotalArr, "TM01", eTimeDay);
      F = FreqPeriodChange(Fin);
    }
    RecS.VarName2 = "mean wave frequency";
    RecS.minval = 0;
    RecS.maxval = 0.9;
    RecS.mindiff = -0.2;
    RecS.maxdiff = 0.2;
    RecS.Unit = "Hz";
  }
  if (FullVarName == "PeakWaveFreq") {
    if (eModelName == "COSMO" || eModelName == "WAM")
      F = Get2DvariableSpecTime(TotalArr, "PwaveFreq", eTimeDay);
    if (eModelName == "WWM") {
      MyMatrix<double> Fin = Get2DvariableSpecTime(TotalArr, "TPP", eTimeDay);
      F = FreqPeriodChange(Fin);
    }
    RecS.VarName2 = "peak wave frequency";
    RecS.minval = 0;
    RecS.maxval = 0.9;
    RecS.mindiff = -0.2;
    RecS.maxdiff = 0.2;
    RecS.Unit = "Hz";
  }
  if (FullVarName == "MeanWavePer") {
    if (eModelName == "COSMO" || eModelName == "WAM") {
      MyMatrix<double> Fin =
          Get2DvariableSpecTime(TotalArr, "MwaveFreq", eTimeDay);
      F = FreqPeriodChange(Fin);
    }
    if (eModelName == "WWM")
      F = Get2DvariableSpecTime(TotalArr, "TM01", eTimeDay);
    if (eModelName == "WW3")
      F = Get2DvariableSpecTime(TotalArr, "t01", eTimeDay);
    RecS.VarName2 = "mean wave period";
    RecS.minval = 2;
    RecS.maxval = 10;
    RecS.mindiff = -1;
    RecS.maxdiff = 1;
    RecS.Unit = "s";
  }
  if (FullVarName == "PeakWavePer") {
    if (eModelName == "COSMO" || eModelName == "WAM") {
      MyMatrix<double> Fin =
          Get2DvariableSpecTime(TotalArr, "PwaveFreq", eTimeDay);
      F = FreqPeriodChange(Fin);
    }
    if (eModelName == "WWM")
      F = Get2DvariableSpecTime(TotalArr, "TPP", eTimeDay);
    RecS.VarName2 = "peak wave period";
    RecS.minval = 2;
    RecS.maxval = 10;
    RecS.mindiff = -1;
    RecS.maxdiff = 1;
    RecS.Unit = "s";
  }
  if (FullVarName == "TM02") {
    if (eModelName == "WWM")
      F = Get2DvariableSpecTime(TotalArr, "TM02", eTimeDay);
    if (eModelName == "WW3")
      F = Get2DvariableSpecTime(TotalArr, "t02", eTimeDay);
    RecS.VarName2 = "zero crossing wave period";
    RecS.minval = 2;
    RecS.maxval = 10;
    RecS.mindiff = -1;
    RecS.maxdiff = 1;
    RecS.Unit = "s";
  }
  if (FullVarName == "DynBathy") {
    if (eModelName == "UNRUNOFF")
      F = Get2DvariableSpecTime(TotalArr, "H", eTimeDay);
    if (eModelName == "WWM")
      F = Get2DvariableSpecTime(TotalArr, "DW", eTimeDay);
    if (eModelName == "WW3") {
      F = Get2DvariableSpecTime(TotalArr, "wlv", eTimeDay);
      F += GetDEP(TotalArr.GrdArr.GrdArrRho);
    }
    RecS.VarName2 = "dynamic bathymetry";
    RecS.minval = 0;
    RecS.maxval = 30;
    RecS.mindiff = -5;
    RecS.maxdiff = 5;
    RecS.Unit = "m";
  }
  if (FullVarName == "Bathymetry") {
    RecS.VarName2 = "bathymetry";
    RecS.minval = 0;
    RecS.maxval = 1;
    RecS.mindiff = -50;
    RecS.maxdiff = 50;
    RecS.Unit = "m";
    //
    std::vector<std::string> ListModel = {"ROMS", "UNRUNOFF", "WWM", "NEMO"};
    if (PositionVect(ListModel, eModelName) != -1) {
      F = GetDEP(TotalArr.GrdArr.GrdArrRho);
      PairMinMax ePair = ComputeMinMax(TotalArr.GrdArr, F);
      RecS.maxval = ePair.TheMax;
    }
  }
  if (FullVarName == "MeanWaveDirSpread") {
    if (eModelName == "WWM")
      F = Get2DvariableSpecTime(TotalArr, "DSPR", eTimeDay);
    RecS.VarName2 = "directional spreading";
    RecS.minval = 0;
    RecS.maxval = 30;
    RecS.mindiff = -5;
    RecS.maxdiff = 5;
    RecS.Unit = "deg";
  }
  if (FullVarName == "PeakWaveDirSpread") {
    if (eModelName == "WWM")
      F = Get2DvariableSpecTime(TotalArr, "PEAKDSPR", eTimeDay);
    RecS.VarName2 = "peak directional spreading";
    RecS.minval = 0;
    RecS.maxval = 30;
    RecS.mindiff = -5;
    RecS.maxdiff = 5;
    RecS.Unit = "deg";
  }
  if (FullVarName == "AirZ0") {
    if (eModelName == "WWM")
      F = Get2DvariableSpecTime(TotalArr, "Z0", eTimeDay);
    RecS.VarName2 = "air roughness length";
    RecS.minval = 0;
    RecS.maxval = 0.0002;
    RecS.mindiff = -0.00005;
    RecS.maxdiff = 0.00005;
    RecS.Unit = "m";
  }
  if (FullVarName == "AirFricVel") {
    if (eModelName == "WWM")
      F = Get2DvariableSpecTime(TotalArr, "UFRIC", eTimeDay);
    RecS.VarName2 = "air roughness length";
    RecS.minval = 0;
    RecS.maxval = 0.3;
    RecS.mindiff = -0.05;
    RecS.maxdiff = 0.05;
    RecS.Unit = "m";
  }
  if (FullVarName == "CdWave") {
    if (eModelName == "COSMO" || eModelName == "WAM")
      F = Get2DvariableSpecTime(TotalArr, "CdWave", eTimeDay);
    if (eModelName == "WWM")
      F = Get2DvariableSpecTime(TotalArr, "CD", eTimeDay);
    RecS.VarName2 = "drag coefficient from the wave model";
    RecS.minval = 0;
    RecS.maxval = 0.20;
    RecS.mindiff = -0.05;
    RecS.maxdiff = 0.05;
    RecS.Unit = "nondim.";
  }
  if (FullVarName == "AlphaWave") {
    if (eModelName == "COSMO" || eModelName == "WAM")
      F = Get2DvariableSpecTime(TotalArr, "AlphaWave", eTimeDay);
    if (eModelName == "WWM")
      F = Get2DvariableSpecTime(TotalArr, "ALPHA_CH", eTimeDay);
    RecS.VarName2 = "Charnock coefficient from the wave model";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.Unit = "nondim.";
  }
  if (FullVarName == "TotSurfStr") {
    if (eModelName == "WWM")
      F = Get2DvariableSpecTime(TotalArr, "TAUTOT", eTimeDay);
    RecS.VarName2 = "Total Surface stress";
    RecS.minval = 0;
    RecS.maxval = 0.06;
    RecS.mindiff = -0.01;
    RecS.maxdiff = 0.01;
    RecS.Unit = "unknown";
  }
  if (FullVarName == "WaveSurfStr") {
    if (eModelName == "WWM")
      F = Get2DvariableSpecTime(TotalArr, "TAUW", eTimeDay);
    RecS.VarName2 = "wave supported Surface stress";
    RecS.minval = 0;
    RecS.maxval = 0.06;
    RecS.mindiff = -0.01;
    RecS.maxdiff = 0.01;
    RecS.Unit = "unknown";
  }
  if (FullVarName == "SurfStrHF") {
    if (eModelName == "WWM")
      F = Get2DvariableSpecTime(TotalArr, "TAUHF", eTimeDay);
    RecS.VarName2 = "high frequency Surface stress";
    RecS.minval = 0;
    RecS.maxval = 0.06;
    RecS.mindiff = -0.01;
    RecS.maxdiff = 0.01;
    RecS.Unit = "unknown";
  }
  //
  // BFM model variables.
  //
  if (FullVarName == "oxygen") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "oxygen", eTimeDay);
    if (eModelName == "NEMO")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "o2", eTimeDay);
    RecS.VarName2 = "dissolved oxygen concentration";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mm/m3";
  }
  if (FullVarName == "PO4") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "PO4", eTimeDay);
    if (eModelName == "NEMO")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "po4", eTimeDay);
    RecS.VarName2 = "phosphate concentration";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mm/m3";
    RecS.varName_ROMS = "PO4";
  }
  if (FullVarName == "pH") {
    if (eModelName == "NEMO")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "ph", eTimeDay);
    RecS.VarName2 = "pH";
    RecS.minval = 0;
    RecS.maxval = 14;
    RecS.mindiff = -1;
    RecS.maxdiff = 1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "nondim";
    RecS.varName_ROMS = "pH";
  }
  if (FullVarName == "alkalinity") {
    if (eModelName == "NEMO")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "talk", eTimeDay);
    RecS.VarName2 = "alkalinity";
    RecS.minval = 2;
    RecS.maxval = 3;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mol/m3";
    RecS.varName_ROMS = "talk";
  }
  if (FullVarName == "NO3") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "NO3", eTimeDay);
    if (eModelName == "NEMO")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "no3", eTimeDay);
    RecS.VarName2 = "nitrate concentration";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mm/m3";
    RecS.varName_ROMS = "NO3";
  }
  if (FullVarName == "NH4") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "NH4", eTimeDay);
    if (eModelName == "NEMO")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "nh4", eTimeDay);
    RecS.VarName2 = "ammonium concentration";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mm/m3";
    RecS.varName_ROMS = "NH4";
  }
  if (FullVarName == "NitrogenSink") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "NitrogenSink", eTimeDay);
    RecS.VarName2 = "aerobic and anaerobic bacteria";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mm/m3";
    RecS.varName_ROMS = "NitrogenSink";
  }
  if (FullVarName == "SiOH4") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "SiOH4", eTimeDay);
    if (eModelName == "NEMO")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "si", eTimeDay);
    RecS.VarName2 = "silicate concentration";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mm/m3";
    RecS.varName_ROMS = "SiOH4";
  }
  if (FullVarName == "ReductionEquivalent") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "ReductionEquivalent",
                                           eTimeDay);
    RecS.VarName2 = "reduction equivalent";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mm/m3";
    RecS.varName_ROMS = "ReductionEquivalent";
  }
  if (FullVarName == "bacteriaC") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "bacteria_c", eTimeDay);
    if (eModelName == "NEMO") {
      double fact = (22.0 / 0.03) * 0.03 / (0.03 + 0.01 + 0.01 + 0.07);
      Tens3 = fact * NETCDF_Get3DvariableSpecTime(TotalArr, "chl", eTimeDay);
    }
    RecS.VarName2 = "aerobic and anaerobic bacteria(C)";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mg/m3";
    RecS.varName_ROMS = "bacteria_c";
  }
  if (FullVarName == "bacteriaN") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "bacteria_n", eTimeDay);
    if (eModelName == "NEMO") {
      double fact = (0.017 / 0.03) * 0.03 / (0.03 + 0.01 + 0.01 + 0.07);
      Tens3 = fact * NETCDF_Get3DvariableSpecTime(TotalArr, "chl", eTimeDay);
    }
    RecS.VarName2 = "aerobic and anaerobic bacteria(N)";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mg/m3";
    RecS.varName_ROMS = "bacteria_n";
  }
  if (FullVarName == "bacteriaP") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "bacteria_p", eTimeDay);
    if (eModelName == "NEMO") {
      double fact = (0.0019 / 0.03) * 0.03 / (0.03 + 0.01 + 0.01 + 0.07);
      Tens3 = fact * NETCDF_Get3DvariableSpecTime(TotalArr, "chl", eTimeDay);
    }
    RecS.VarName2 = "aerobic and anaerobic bacteria(P)";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mg/m3";
    RecS.varName_ROMS = "bacteria_p";
  }
  if (FullVarName == "diatomsC") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "diatoms_c", eTimeDay);
    if (eModelName == "NEMO") {
      double fact = (5.2 / 0.03) * 0.03 / (0.03 + 0.01 + 0.01 + 0.07);
      Tens3 = fact * NETCDF_Get3DvariableSpecTime(TotalArr, "chl", eTimeDay);
    }
    RecS.VarName2 = "diatoms carbon";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mg/m3";
    RecS.varName_ROMS = "diatoms_c";
  }
  if (FullVarName == "diatomsL") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "diatoms_l", eTimeDay);
    if (eModelName == "NEMO") {
      double fact = 0.03 / (0.03 + 0.01 + 0.01 + 0.07);
      Tens3 = fact * NETCDF_Get3DvariableSpecTime(TotalArr, "chl", eTimeDay);
    }
    RecS.VarName2 = "diatoms chlorophyl";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mg/m3";
    RecS.varName_ROMS = "diatoms_l";
  }
  if (FullVarName == "diatomsN") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "diatoms_n", eTimeDay);
    if (eModelName == "NEMO") {
      double fact = (0.07 / 0.03) * 0.03 / (0.03 + 0.01 + 0.01 + 0.07);
      Tens3 = fact * NETCDF_Get3DvariableSpecTime(TotalArr, "chl", eTimeDay);
    }
    RecS.VarName2 = "diatoms nitrogen";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mm/m3";
    RecS.varName_ROMS = "diatoms_n";
  }
  if (FullVarName == "diatomsP") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "diatoms_p", eTimeDay);
    if (eModelName == "NEMO") {
      double fact = (0.006 / 0.03) * 0.03 / (0.03 + 0.01 + 0.01 + 0.07);
      Tens3 = fact * NETCDF_Get3DvariableSpecTime(TotalArr, "chl", eTimeDay);
    }
    RecS.VarName2 = "diatoms phosphate";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mm/m3";
    RecS.varName_ROMS = "diatoms_p";
  }
  if (FullVarName == "diatomsS") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "diatoms_s", eTimeDay);
    if (eModelName == "NEMO") {
      double fact = (0.046 / 0.03) * 0.03 / (0.03 + 0.01 + 0.01 + 0.07);
      Tens3 = fact * NETCDF_Get3DvariableSpecTime(TotalArr, "chl", eTimeDay);
    }
    RecS.VarName2 = "diatoms silicate";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mm/m3";
    RecS.varName_ROMS = "diatoms_s";
  }
  if (FullVarName == "flagellatesC") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "flagellates_c", eTimeDay);
    if (eModelName == "NEMO") {
      double fact = (4.23 / 0.01) * 0.01 / (0.03 + 0.01 + 0.01 + 0.07);
      Tens3 = fact * NETCDF_Get3DvariableSpecTime(TotalArr, "chl", eTimeDay);
    }
    RecS.VarName2 = "flagellates carbon";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mg/m3";
    RecS.varName_ROMS = "flagellates_c";
  }
  if (FullVarName == "flagellatesN") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "flagellates_n", eTimeDay);
    if (eModelName == "NEMO") {
      double fact = (0.03 / 0.01) * 0.01 / (0.03 + 0.01 + 0.01 + 0.07);
      Tens3 = fact * NETCDF_Get3DvariableSpecTime(TotalArr, "chl", eTimeDay);
    }
    RecS.VarName2 = "flagellates nitrogen";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mm/m3";
    RecS.varName_ROMS = "flagellates_n";
  }
  if (FullVarName == "flagellatesP") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "flagellates_p", eTimeDay);
    if (eModelName == "NEMO") {
      double fact = (0.004 / 0.01) * 0.01 / (0.03 + 0.01 + 0.01 + 0.07);
      Tens3 = fact * NETCDF_Get3DvariableSpecTime(TotalArr, "chl", eTimeDay);
    }
    RecS.VarName2 = "flagellates phosphate";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mm/m3";
    RecS.varName_ROMS = "flagellates_p";
  }
  if (FullVarName == "flagellatesL") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "flagellates_l", eTimeDay);
    if (eModelName == "NEMO") {
      double fact = 0.01 / (0.03 + 0.01 + 0.01 + 0.07);
      Tens3 = fact * NETCDF_Get3DvariableSpecTime(TotalArr, "chl", eTimeDay);
    }
    RecS.VarName2 = "flagellates chlorophyl";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mg/m3";
    RecS.varName_ROMS = "flagellates_l";
  }
  if (FullVarName == "phytoplanktonbiomass") {
    if (eModelName == "NEMO") {
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "phyc", eTimeDay);
    }
    RecS.VarName2 = "phytoplankton biomass";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mmol/m3";
  }
  if (FullVarName == "picophytoplanktonC") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "picophytoplankton_c",
                                           eTimeDay);
    if (eModelName == "NEMO") {
      double fact = 15.6 / (0.03 + 0.01 + 0.01 + 0.07);
      Tens3 = fact * NETCDF_Get3DvariableSpecTime(TotalArr, "chl", eTimeDay);
    }
    RecS.VarName2 = "picophytoplankton carbon";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mg/m3";
    RecS.varName_ROMS = "picophytoplankton_c";
  }
  if (FullVarName == "picophytoplanktonN") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "picophytoplankton_n",
                                           eTimeDay);
    if (eModelName == "NEMO") {
      double fact = 0.1 / (0.03 + 0.01 + 0.01 + 0.07);
      Tens3 = fact * NETCDF_Get3DvariableSpecTime(TotalArr, "chl", eTimeDay);
    }
    RecS.VarName2 = "picophytoplankton nitrogen";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mm/m3";
    RecS.varName_ROMS = "picophytoplankton_n";
  }
  if (FullVarName == "picophytoplanktonP") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "picophytoplankton_p",
                                           eTimeDay);
    if (eModelName == "NEMO") {
      double fact = 0.015 / (0.03 + 0.01 + 0.01 + 0.07);
      Tens3 = fact * NETCDF_Get3DvariableSpecTime(TotalArr, "chl", eTimeDay);
    }
    RecS.VarName2 = "picophytoplankton phosphate";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mm/m3";
    RecS.varName_ROMS = "picophytoplankton_p";
  }
  if (FullVarName == "picophytoplanktonL") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "picophytoplankton_l",
                                           eTimeDay);
    if (eModelName == "NEMO") {
      double fact = 0.01 / (0.03 + 0.01 + 0.01 + 0.07);
      Tens3 = fact * NETCDF_Get3DvariableSpecTime(TotalArr, "chl", eTimeDay);
    }
    RecS.VarName2 = "picophytoplankton chlorophyl";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mg/m3";
    RecS.varName_ROMS = "picophytoplankton_l";
  }
  if (FullVarName == "largephytoplanktonC") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "largephytoplankton_c",
                                           eTimeDay);
    if (eModelName == "NEMO") {
      double fact = 15.0 / (0.03 + 0.01 + 0.01 + 0.07);
      Tens3 = fact * NETCDF_Get3DvariableSpecTime(TotalArr, "chl", eTimeDay);
    }
    RecS.VarName2 = "largephytoplankton carbon";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mg/m3";
    RecS.varName_ROMS = "largephytoplankton_c";
  }
  if (FullVarName == "largephytoplanktonN") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "largephytoplankton_n",
                                           eTimeDay);
    if (eModelName == "NEMO") {
      double fact = 0.11 / (0.03 + 0.01 + 0.01 + 0.07);
      Tens3 = fact * NETCDF_Get3DvariableSpecTime(TotalArr, "chl", eTimeDay);
    }
    RecS.VarName2 = "largephytoplankton nitrogen";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mm/m3";
    RecS.varName_ROMS = "largephytoplankton_n";
  }
  if (FullVarName == "largephytoplanktonP") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "largephytoplankton_p",
                                           eTimeDay);
    if (eModelName == "NEMO") {
      double fact = 0.025 / (0.03 + 0.01 + 0.01 + 0.07);
      Tens3 = fact * NETCDF_Get3DvariableSpecTime(TotalArr, "chl", eTimeDay);
    }
    RecS.VarName2 = "largephytoplankton phosphate";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mm/m3";
    RecS.varName_ROMS = "largephytoplankton_p";
  }
  if (FullVarName == "largephytoplanktonL") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "largephytoplankton_l",
                                           eTimeDay);
    if (eModelName == "NEMO") {
      double fact = 0.07 / (0.03 + 0.01 + 0.01 + 0.07);
      Tens3 = fact * NETCDF_Get3DvariableSpecTime(TotalArr, "chl", eTimeDay);
    }
    RecS.VarName2 = "largephytoplankton chlorophyl";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mg/m3";
    RecS.varName_ROMS = "largephytoplankton_l";
  }
  if (FullVarName == "CarnPlanktonC") {
    if (eModelName == "ROMS")
      Tens3 =
          NETCDF_Get3DvariableSpecTime(TotalArr, "CarnPlankton_c", eTimeDay);
    if (eModelName == "NEMO") {
      double fact = 4.2 / (0.03 + 0.01 + 0.01 + 0.07);
      Tens3 = fact * NETCDF_Get3DvariableSpecTime(TotalArr, "chl", eTimeDay);
    }
    RecS.VarName2 = "Carnivorous Mesozooplankton(C)";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mg/m3";
    RecS.varName_ROMS = "CarnPlankton_c";
  }
  if (FullVarName == "CarnPlanktonN") {
    if (eModelName == "ROMS")
      Tens3 =
          NETCDF_Get3DvariableSpecTime(TotalArr, "CarnPlankton_n", eTimeDay);
    if (eModelName == "NEMO") {
      double fact = 0.014 / (0.03 + 0.01 + 0.01 + 0.07);
      Tens3 = fact * NETCDF_Get3DvariableSpecTime(TotalArr, "chl", eTimeDay);
    }
    RecS.VarName2 = "Carnivorous Mesozooplankton(N)";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mg/m3";
    RecS.varName_ROMS = "CarnPlankton_n";
  }
  if (FullVarName == "CarnPlanktonP") {
    if (eModelName == "ROMS")
      Tens3 =
          NETCDF_Get3DvariableSpecTime(TotalArr, "CarnPlankton_p", eTimeDay);
    if (eModelName == "NEMO") {
      double fact = 0.0016 / (0.03 + 0.01 + 0.01 + 0.07);
      Tens3 = fact * NETCDF_Get3DvariableSpecTime(TotalArr, "chl", eTimeDay);
    }
    RecS.VarName2 = "Carnivorous Mesozooplankton(P)";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mg/m3";
    RecS.varName_ROMS = "CarnPlankton_p";
  }
  if (FullVarName == "OmniPlanktonC") {
    if (eModelName == "ROMS")
      Tens3 =
          NETCDF_Get3DvariableSpecTime(TotalArr, "OmniPlankton_c", eTimeDay);
    if (eModelName == "NEMO") {
      double fact = 6.29 / (0.03 + 0.01 + 0.01 + 0.07);
      Tens3 = fact * NETCDF_Get3DvariableSpecTime(TotalArr, "chl", eTimeDay);
    }
    RecS.VarName2 = "Omnivorous Mesozooplankton(C)";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mg/m3";
    RecS.varName_ROMS = "OmniPlankton_c";
  }
  if (FullVarName == "OmniPlanktonN") {
    if (eModelName == "ROMS")
      Tens3 =
          NETCDF_Get3DvariableSpecTime(TotalArr, "OmniPlankton_n", eTimeDay);
    if (eModelName == "NEMO") {
      double fact = 0.015 / (0.03 + 0.01 + 0.01 + 0.07);
      Tens3 = fact * NETCDF_Get3DvariableSpecTime(TotalArr, "chl", eTimeDay);
    }
    RecS.VarName2 = "Omnivorous Mesozooplankton(N)";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mg/m3";
    RecS.varName_ROMS = "OmniPlankton_n";
  }
  if (FullVarName == "OmniPlanktonP") {
    if (eModelName == "ROMS")
      Tens3 =
          NETCDF_Get3DvariableSpecTime(TotalArr, "OmniPlankton_p", eTimeDay);
    if (eModelName == "NEMO") {
      double fact = 0.0016 / (0.03 + 0.01 + 0.01 + 0.07);
      Tens3 = fact * NETCDF_Get3DvariableSpecTime(TotalArr, "chl", eTimeDay);
    }
    RecS.VarName2 = "Omnivorous Mesozooplankton(P)";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mg/m3";
    RecS.varName_ROMS = "OmniPlankton_p";
  }
  if (FullVarName == "MicroPlanktonC") {
    if (eModelName == "ROMS")
      Tens3 =
          NETCDF_Get3DvariableSpecTime(TotalArr, "MicroPlankton_c", eTimeDay);
    if (eModelName == "NEMO") {
      double fact = 17.5 / (0.03 + 0.01 + 0.01 + 0.07);
      Tens3 = fact * NETCDF_Get3DvariableSpecTime(TotalArr, "chl", eTimeDay);
    }
    RecS.VarName2 = "Microzooplankton(C)";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mg/m3";
    RecS.varName_ROMS = "MicroPlankton_c";
  }
  if (FullVarName == "MicroPlanktonN") {
    if (eModelName == "ROMS")
      Tens3 =
          NETCDF_Get3DvariableSpecTime(TotalArr, "MicroPlankton_n", eTimeDay);
    if (eModelName == "NEMO") {
      double fact = 0.016 / (0.03 + 0.01 + 0.01 + 0.07);
      Tens3 = fact * NETCDF_Get3DvariableSpecTime(TotalArr, "chl", eTimeDay);
    }
    RecS.VarName2 = "Microzooplankton(N)";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mg/m3";
    RecS.varName_ROMS = "MicroPlankton_n";
  }
  if (FullVarName == "MicroPlanktonP") {
    if (eModelName == "ROMS")
      Tens3 =
          NETCDF_Get3DvariableSpecTime(TotalArr, "MicroPlankton_p", eTimeDay);
    if (eModelName == "NEMO") {
      double fact = 0.0018 / (0.03 + 0.01 + 0.01 + 0.07);
      Tens3 = fact * NETCDF_Get3DvariableSpecTime(TotalArr, "chl", eTimeDay);
    }
    RecS.VarName2 = "Microzooplankton(P)";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mg/m3";
    RecS.varName_ROMS = "MicroPlankton_p";
  }
  if (FullVarName == "HeteroNanoflagelattesC") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "HeteroNanoflagelattes_c",
                                           eTimeDay);
    if (eModelName == "NEMO") {
      double fact = 10.9 / (0.03 + 0.01 + 0.01 + 0.07);
      Tens3 = fact * NETCDF_Get3DvariableSpecTime(TotalArr, "chl", eTimeDay);
    }
    RecS.VarName2 = "Heterotrophic Nanoflagellates (HNAN) C";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mg/m3";
    RecS.varName_ROMS = "HeteroNanoflagelattes_c";
  }
  if (FullVarName == "HeteroNanoflagelattesN") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "HeteroNanoflagelattes_n",
                                           eTimeDay);
    if (eModelName == "NEMO") {
      double fact = 0.016 / (0.03 + 0.01 + 0.01 + 0.07);
      Tens3 = fact * NETCDF_Get3DvariableSpecTime(TotalArr, "chl", eTimeDay);
    }
    RecS.VarName2 = "Heterotrophic Nanoflagellates (HNAN) N";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mg/m3";
    RecS.varName_ROMS = "HeteroNanoflagelattes_n";
  }
  if (FullVarName == "HeteroNanoflagelattesP") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "HeteroNanoflagelattes_p",
                                           eTimeDay);
    if (eModelName == "NEMO") {
      double fact = 0.0018 / (0.03 + 0.01 + 0.01 + 0.07);
      Tens3 = fact * NETCDF_Get3DvariableSpecTime(TotalArr, "chl", eTimeDay);
    }
    RecS.VarName2 = "Heterotrophic Nanoflagellates (HNAN) P";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mg/m3";
    RecS.varName_ROMS = "HeteroNanoflagelattes_p";
  }
  if (FullVarName == "LabileDOM1c") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "LabileDOM1_c", eTimeDay);
    RecS.VarName2 = "Labile Dissolved Organic Matter (Carbon)";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mg/m3";
    RecS.varName_ROMS = "LabileDOM1_c";
  }
  if (FullVarName == "LabileDOM1n") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "LabileDOM1_n", eTimeDay);
    RecS.VarName2 = "Labile Dissolved Organic Matter (Nitrogen)";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mg/m3";
    RecS.varName_ROMS = "LabileDOM1_n";
  }
  if (FullVarName == "LabileDOM1p") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "LabileDOM1_p", eTimeDay);
    RecS.VarName2 = "Labile Dissolved Organic Matter (Phosphate)";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mg/m3";
    RecS.varName_ROMS = "LabileDOM1_p";
  }
  if (FullVarName == "LabileDOM2c") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "LabileDOM2_c", eTimeDay);
    RecS.VarName2 = "Semi-labile Dissolved Organic Carbon";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mg/m3";
    RecS.varName_ROMS = "LabileDOM2_c";
  }
  if (FullVarName == "RefractoryDOMc") {
    if (eModelName == "ROMS")
      Tens3 =
          NETCDF_Get3DvariableSpecTime(TotalArr, "RefractoryDOM_c", eTimeDay);
    RecS.VarName2 = "Semi-refractory Dissolved Organic Carbon";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mg/m3";
    RecS.varName_ROMS = "RefractoryDOM_c";
  }
  if (FullVarName == "DissolvedOMc") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "DissolvedOM_c", eTimeDay);
    RecS.VarName2 = "Dissolved Organic Matter (Carbon)";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mg/m3";
    RecS.varName_ROMS = "DissolvedOM_c";
  }
  if (FullVarName == "ParticleOMc") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "ParticleOM_c", eTimeDay);
    RecS.VarName2 = "Particle Organic Matter (Carbon)";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mg/m3";
    RecS.varName_ROMS = "ParticleOM_c";
  }
  if (FullVarName == "ParticleOMn") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "ParticleOM_n", eTimeDay);
    RecS.VarName2 = "Particle Organic Matter (Nitrogen)";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mm/m3";
    RecS.varName_ROMS = "ParticleOM_n";
  }
  if (FullVarName == "ParticleOMp") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "ParticleOM_p", eTimeDay);
    RecS.VarName2 = "Particle Organic Matter (Phosphate)";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mm/m3";
    RecS.varName_ROMS = "ParticleOM_p";
  }
  if (FullVarName == "ParticleOMs") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "ParticleOM_s", eTimeDay);
    RecS.VarName2 = "Particle Organic Matter (Silicate)";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mm/m3";
    RecS.varName_ROMS = "ParticleOM_s";
  }
  if (FullVarName == "DissolvedICc") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "DissolvedIC_c", eTimeDay);
    RecS.VarName2 = "Dissolved Inorganic Carbon C";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mg/m3";
    RecS.varName_ROMS = "DissolvedIC_c";
  }
  if (FullVarName == "DissolvedICh") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "DissolvedIC_h", eTimeDay);
    RecS.VarName2 = "Dissolved Inorganic Carbon H";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mm/m3";
    RecS.varName_ROMS = "DissolvedIC_h";
  }
  if (FullVarName == "Irradiance") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "Irradiance", eTimeDay);
    RecS.VarName2 = "Irradiance";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "uE/m2s";
    RecS.varName_ROMS = "Irradiance";
  }
  if (FullVarName == "DIC") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "DIC", eTimeDay);
    if (eModelName == "NEMO")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "dissic", eTimeDay);
    RecS.VarName2 = "DIC concentration";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mm/m3";
    RecS.varName_ROMS = "DIC";
  }
  if (FullVarName == "chlorophyll") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "chlorophyll", eTimeDay);
    if (eModelName == "NEMO")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "chl", eTimeDay);
    RecS.VarName2 = "chlorophyll concentration";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mg/m3";
    RecS.varName_ROMS = "chlorophyll";
  }
  if (FullVarName == "NetProductionP1") {
    if (eModelName == "ROMS")
      Tens3 =
          NETCDF_Get3DvariableSpecTime(TotalArr, "NetProductionP1", eTimeDay);
    RecS.VarName2 = "Specific Net Production of P1(Phytoplankton)";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "1/d";
    RecS.varName_ROMS = "NetProductionP1";
  }
  if (FullVarName == "NetProductionP2") {
    if (eModelName == "ROMS")
      Tens3 =
          NETCDF_Get3DvariableSpecTime(TotalArr, "NetProductionP2", eTimeDay);
    RecS.VarName2 = "Specific Net Production of P2(Phytoplankton)";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "1/d";
    RecS.varName_ROMS = "NetProductionP2";
  }
  if (FullVarName == "NetProductionP3") {
    if (eModelName == "ROMS")
      Tens3 =
          NETCDF_Get3DvariableSpecTime(TotalArr, "NetProductionP3", eTimeDay);
    RecS.VarName2 = "Specific Net Production of P3(Phytoplankton)";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "1/d";
    RecS.varName_ROMS = "NetProductionP3";
  }
  if (FullVarName == "NetProductionP4") {
    if (eModelName == "ROMS")
      Tens3 =
          NETCDF_Get3DvariableSpecTime(TotalArr, "NetProductionP4", eTimeDay);
    RecS.VarName2 = "Specific Net Production of P4(Phytoplankton)";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "1/d";
    RecS.varName_ROMS = "NetProductionP4";
  }
  if (FullVarName == "RegFactorP1") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "RegFactorP1", eTimeDay);
    RecS.VarName2 = "Regular Factor for Light in P1(Phytoplankton)";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "nondim.";
    RecS.varName_ROMS = "RegFactorP1";
  }
  if (FullVarName == "RegFactorP2") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "RegFactorP2", eTimeDay);
    RecS.VarName2 = "Regular Factor for Light in P2(Phytoplankton)";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "nondim.";
    RecS.varName_ROMS = "RegFactorP2";
  }
  if (FullVarName == "RegFactorP3") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "RegFactorP3", eTimeDay);
    RecS.VarName2 = "Regular Factor for Light in P3(Phytoplankton)";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "nondim.";
    RecS.varName_ROMS = "RegFactorP3";
  }
  if (FullVarName == "RegFactorP4") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "RegFactorP4", eTimeDay);
    RecS.VarName2 = "Regular Factor for Light in P4(Phytoplankton)";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "nondim.";
    RecS.varName_ROMS = "RegFactorP4";
  }
  if (FullVarName == "GrossPP") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "GrossPP", eTimeDay);
    if (eModelName == "NEMO")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "nppv", eTimeDay);
    RecS.VarName2 = "Gross Primary Production";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mg/m3d";
    RecS.varName_ROMS = "GrossPP";
  }
  if (FullVarName == "SecondPP") {
    if (eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "SecondPP", eTimeDay);
    RecS.VarName2 = "Second Primary Production";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "mg/m3d";
    RecS.varName_ROMS = "SecondPP";
  }
  if (FullVarName == "ExtinctionCoeff") {
    if (eModelName == "ROMS")
      Tens3 =
          NETCDF_Get3DvariableSpecTime(TotalArr, "ExtinctionCoeff", eTimeDay);
    RecS.VarName2 = "Total Extinction Coefficient";
    RecS.minval = 0;
    RecS.maxval = 0.033;
    RecS.mindiff = -0.1;
    RecS.maxdiff = 0.1;
    RecS.VarNature = "3Drho";
    RecS.Unit = "1/m";
    RecS.varName_ROMS = "ExtinctionCoeff";
  }
  //
  // Now error parsing and assignations
  //
  RecVar eRecVar;
  eRecVar.RecS = RecS;
  if (RecS.VarName2 == "unset") {
    std::cerr << "We did not find the variable\n";
    std::cerr << "eVarName = " << eVarName << "\n";
    std::cerr << "in the list of allowed ones\n";
    std::cerr << "possibly missspelling or lack of relevant code\n";
    throw TerminalException{1};
  }
  if (eModelName != "TRIVIAL") {
    if (RecS.VarNature == "rho") {
      if (F.size() == 0) {
        std::cerr << "VarNature = " << RecS.VarNature << "\n";
        std::cerr << "Variable eVarName = " << eVarName << "\n";
        std::cerr << "is recognized by the program as authorized variable\n";
        std::cerr << "But it has not been assigned.\n";
        std::cerr << "Possibly because of missing facility for\n";
        std::cerr << "eModelName = " << eModelName << "\n";
        throw TerminalException{1};
      }
      eRecVar.F = F;
    }
    if (RecS.VarNature == "uv") {
      if (U.size() == 0 || V.size() == 0) {
        std::cerr << "VarNature = " << RecS.VarNature << "\n";
        std::cerr << "Variable eVarName = " << eVarName << "\n";
        std::cerr << "is recognized by the program\n";
        std::cerr << "But it has not been assigned.\n";
        std::cerr << "Possibly because of missing facility for\n";
        std::cerr << "eModelName = " << eModelName << "\n";
        throw TerminalException{1};
      }
      eRecVar.U = U;
      eRecVar.V = V;
      if (!IsEqualSizeMatrices(U, V)) {
        std::cerr << "Matrices U and V are not equal sized\n";
        std::cerr << "|U|=" << StringSizeMatrix(U) << "\n";
        std::cerr << "|V|=" << StringSizeMatrix(V) << "\n";
        throw TerminalException{1};
      }
      eRecVar.F = GetNormMatrix(U, V);
    }
    if (RecS.VarNature == "3Drho") {
      auto LDim = Tens3.dimensions();
      if (LDim[0] == 0) {
        std::cerr << "VarNature = " << RecS.VarNature << "\n";
        std::cerr << "Variable eVarName = " << eVarName << "\n";
        std::cerr << "is recognized by the program\n";
        std::cerr << "But it has not been assigned.\n";
        std::cerr << "Possibly because of missing facility for\n";
        std::cerr << "eModelName = " << eModelName << "\n";
        throw TerminalException{1};
      }
      eRecVar.Tens3 = Tens3;
    }
    if (RecS.VarNature == "3Duv") {
      auto LDim = Uthree.dimensions();
      if (LDim[0] == 0) {
        std::cerr << "VarNature = " << RecS.VarNature << "\n";
        std::cerr << "Variable eVarName = " << eVarName << "\n";
        std::cerr << "is recognized by the program\n";
        std::cerr << "But it has not been assigned.\n";
        std::cerr << "Possibly because of missing facility for\n";
        std::cerr << "eModelName = " << eModelName << "\n";
        throw TerminalException{1};
      }
      eRecVar.Tens3 = ComputeNormPairOfTensor(Uthree, Vthree);
      eRecVar.Uthree = std::move(Uthree);
      eRecVar.Vthree = std::move(Vthree);
    }
  }
  return eRecVar;
}

RecVar ModelSpecificVarSpecificTime(TotalArrGetData const &TotalArr,
                                    std::string const &eVarName,
                                    double const &eTimeDay) {
  std::string eSep = "_";
  std::vector<std::string> ListStr = STRING_Split(eVarName, eSep);
  int len = ListStr.size();
  //  std::cerr << "len=" << len << "\n";
  if (len == 1)
    return ModelSpecificVarSpecificTime_Kernel(TotalArr, eVarName, eTimeDay);
  if (len != 2) {
    std::cerr << "We should clearly rethink len. Expected value is 1 or 2\n";
    std::cerr << "But len=" << len << "\n";
    std::cerr << "eVarName=" << eVarName << "\n";
    for (int i = 0; i < len; i++) {
      std::cerr << "i=" << i << " LStr[i]=" << ListStr[i] << "\n";
    }
    throw TerminalException{1};
  }
  std::string eVar_rho = ListStr[0];
  std::string eVar_uv = ListStr[1];
  RecVar RecVar_rho =
      ModelSpecificVarSpecificTime_Kernel(TotalArr, eVar_rho, eTimeDay);
  RecVar RecVar_uv =
      ModelSpecificVarSpecificTime_Kernel(TotalArr, eVar_uv, eTimeDay);
  std::string VarNat_rho = RecVar_rho.RecS.VarNature;
  std::string VarNat_uv = RecVar_uv.RecS.VarNature;
  if (VarNat_rho != "rho" && VarNat_rho != "3Drho") {
    std::cerr << "The RecVar_rho is not a rho type variable. Error!\n";
    std::cerr << "Correct way to call is Var_rho _ Var_uv\n";
    std::cerr << "for Example Hwave_SurfCurr for Hwave as rho variable and "
                 "SurfCurr as uv variable\n";
    std::cerr << "The call was with eVarName=" << eVarName << "\n";
    throw TerminalException{1};
  }
  if (VarNat_uv != "uv" && VarNat_uv != "3Duv") {
    std::cerr << "The RecVar_uv is not a uv type variable. Error!\n";
    std::cerr << "Correct way to call is Var_rho _ Var_uv\n";
    std::cerr << "for Example Hwave_SurfCurr for Hwave as rho variable and "
                 "SurfCurr as uv variable\n";
    std::cerr << "The call was with eVarName=" << eVarName << "\n";
    throw TerminalException{1};
  }
  if ((VarNat_rho == "3Drho" && VarNat_uv == "uv") ||
      (VarNat_rho == "rho" && VarNat_uv == "3Duv")) {
    std::cerr << "Error. variables do not have the same dimensionality\n";
    std::cerr << "It should be both 3D or both 2D\n";
    std::cerr << "Right now, we have VarNat_rho=" << VarNat_rho << "\n";
    std::cerr << "Right now, we have VarNat_uv =" << VarNat_uv << "\n";
    throw TerminalException{1};
  }
  if (VarNat_rho == "rho") {
    RecVar_rho.U = RecVar_uv.U;
    RecVar_rho.V = RecVar_uv.V;
  }
  if (VarNat_rho == "3Drho") {
    RecVar_rho.Uthree = RecVar_uv.Uthree;
    RecVar_rho.Vthree = RecVar_uv.Vthree;
  }
  RecVar_rho.RecS.VarName1 += "_" + RecVar_uv.RecS.VarName1;
  RecVar_rho.RecS.VarName2 += " + " + RecVar_uv.RecS.VarName2;
  RecVar_rho.RecS.VarNature = "uv";
  return RecVar_rho;
}

RecVar RetrieveTrivialRecVar(std::string const &FullVarName) {
  TotalArrGetData TotalArrTrivial;
  TotalArrTrivial.GrdArr.ModelName = "TRIVIAL";
  std::string eVarName = DecltypeString(FullVarName);
  double eTimeDayTrivial = 0;
  return ModelSpecificVarSpecificTime(TotalArrTrivial, eVarName,
                                      eTimeDayTrivial);
}

std::vector<std::string> GetAllPossibleVariables_with_pairs() {
  std::vector<std::string> ListVar = GetAllPossibleVariables();
  std::vector<std::string> ListVar_rho;
  std::vector<std::string> ListVar_uv;
  std::vector<std::string> ListVar_3Drho;
  std::vector<std::string> ListVar_3Duv;
  std::vector<std::string> ListVar_Ret;
  for (auto &eVar : ListVar) {
    RecVar eRec = RetrieveTrivialRecVar(eVar);
    if (eRec.RecS.VarNature == "rho")
      ListVar_rho.push_back(eVar);
    if (eRec.RecS.VarNature == "uv")
      ListVar_uv.push_back(eVar);
    if (eRec.RecS.VarNature == "3Drho")
      ListVar_3Drho.push_back(eVar);
    if (eRec.RecS.VarNature == "3Duv")
      ListVar_3Duv.push_back(eVar);
    ListVar_Ret.push_back(eVar);
  }
  for (auto &eVar_uv : ListVar_uv) {
    for (auto &eVar_rho : ListVar_rho) {
      std::string eVarTot = eVar_rho + "_" + eVar_uv;
      ListVar_Ret.push_back(eVarTot);
    }
  }
  for (auto &eVar_uv : ListVar_3Duv) {
    for (auto &eVar_rho : ListVar_3Drho) {
      std::string eVarTot = eVar_rho + "_" + eVar_uv;
      ListVar_Ret.push_back(eVarTot);
    }
  }
  return ListVar_Ret;
}

double ApplyRounding(double const& the_input, std::string const& the_rounding) {
  auto f_approx=[](double const& val, int const& n_digit) -> double {
    int epow = 1;
    for (int i=0; i<n_digit; i++) {
      epow *= 10;
    }
    double epow_d = static_cast<double>(epow);
    double val2 = epow_d * val;
    int val3 = static_cast<int>(round(val2));
    double val4 = static_cast<double>(val3);
    return val4 / epow_d;
  };
  if (the_rounding == "exact") {
    return the_input;
  }
  if (the_rounding == "onedot") {
    return f_approx(the_input, 1);
  }
  if (the_rounding == "empirical") {
    if (T_abs(the_input) > 10) {
      return f_approx(the_input,0);
    }
    if (the_input > 0.4) {
      double min_err = 20.0;
      double choice = the_input;
      for (int i=1; i< 10; i++) {
        double poss = 0.5 * i;
        double err = T_abs(poss - the_input);
        if (err < min_err) {
          min_err = err;
          choice = poss;
        }
      }
      return choice;
    }
    return f_approx(the_input,1);
  }
  std::cerr << "Failed to find the matching entry in ApplyRounding\n";
  std::cerr << "the_input=" << the_input << " the_rounding=" << the_rounding << "\n";
  throw TerminalException{1};
}

void ApplyPlotBound(TotalArrGetData const &TotalArr, RecVar &eRecVar,
                    std::string const &eVarName, PlotBound const &ePlotBound) {
  //
  // Setting up bounds for the plots.
  //
  int nbSingle = ePlotBound.BoundSingle_var.size();
  int nbSingleMin = ePlotBound.BoundSingle_min.size();
  int nbSingleMax = ePlotBound.BoundSingle_max.size();
  //  std::cerr << "nbSingle=" << nbSingle << "\n";
  if (nbSingle != nbSingleMin || nbSingle != nbSingleMax) {
    std::cerr << "Number of entries in BoundSingle_var, BoundSingle_min, "
                 "BoundSingle_max\n";
    std::cerr << "Should all be the same. Now,\n";
    std::cerr << "nbSingle    = " << nbSingle << "\n";
    std::cerr << "nbSingleMin = " << nbSingleMin << "\n";
    std::cerr << "nbSingleMax = " << nbSingleMax << "\n";
  }
  std::string DeclVarName = DecltypeString(eVarName);
  //  std::cerr << "eVarName=" << eVarName << " DeclVarName=" << DeclVarName
  //            << "\n";
  for (int iS = 0; iS < nbSingle; iS++)
    if (ePlotBound.BoundSingle_var[iS] == eVarName ||
        ePlotBound.BoundSingle_var[iS] == DeclVarName) {
      eRecVar.RecS.minval = ePlotBound.BoundSingle_min[iS];
      eRecVar.RecS.maxval = ePlotBound.BoundSingle_max[iS];
    }
  int nbDiff = ePlotBound.BoundDiff_var.size();
  int nbDiffMin = ePlotBound.BoundDiff_min.size();
  int nbDiffMax = ePlotBound.BoundDiff_max.size();
  //  std::cerr << "nbDiff=" << nbDiff << "\n";
  if (nbDiff != nbDiffMin || nbDiff != nbDiffMax) {
    std::cerr
        << "Number of entries in BoundDiff_var, BoundDiff_min, BoundDiff_max\n";
    std::cerr << "Should all be the same. Now,\n";
    std::cerr << "nbDiff    = " << nbDiff << "\n";
    std::cerr << "nbDiffMin = " << nbDiffMin << "\n";
    std::cerr << "nbDiffMax = " << nbDiffMax << "\n";
  }
  for (int iD = 0; iD < nbDiff; iD++)
    if (ePlotBound.BoundDiff_var[iD] == eVarName) {
      eRecVar.RecS.mindiff = ePlotBound.BoundDiff_min[iD];
      eRecVar.RecS.maxdiff = ePlotBound.BoundDiff_max[iD];
    }
  int eSize = eRecVar.F.size();
  if (ePlotBound.VariableRange && eSize > 0) {
    PairMinMax ePair = ComputeMinMax(TotalArr.GrdArr, eRecVar.F);
    ePair.TheMin = ApplyRounding(ePair.TheMin, ePlotBound.VariableRangeRounding);
    ePair.TheMax = ApplyRounding(ePair.TheMax, ePlotBound.VariableRangeRounding);
    eRecVar.RecS.minval = ePair.TheMin;
    eRecVar.RecS.maxval = ePair.TheMax;
  }
}

void ApplyPlotBoundPair(TotalArrGetData const &TotalArr1,
                        TotalArrGetData const &TotalArr2, RecVar &eRecVar1,
                        RecVar &eRecVar2, std::string const &eVarName,
                        PlotBound const &ePlotBound) {

  int eSize = eRecVar1.F.size();
  //  std::cerr << "ApplyPlotBoundPair : eSize = " << eSize << "\n";
  if (ePlotBound.VariableRange && eSize > 0) {
    MyMatrix<double> eDiff12 = eRecVar1.F - eRecVar2.F;
    PairMinMax ePair = ComputeMinMax(TotalArr1.GrdArr, eDiff12);
    double MaxChange = std::max(ePair.TheMax, -ePair.TheMin);
    MaxChange = ApplyRounding(MaxChange, ePlotBound.VariableRangeRounding);
    //    std::cerr << "ApplyPlotBoundPair : min/max = " << ePair.TheMin << " /
    //    " << ePair.TheMax << " MaxChange=" << MaxChange << "\n";
    eRecVar1.RecS.mindiff = -MaxChange;
    eRecVar1.RecS.maxdiff = MaxChange;
    eRecVar2.RecS.mindiff = -MaxChange;
    eRecVar2.RecS.maxdiff = MaxChange;
  }
}

RecVar ModelSpecificVarSpecificTimeBound(TotalArrGetData const &TotalArr,
                                         std::string const &eVarName,
                                         double const &eTimeDay,
                                         PlotBound const &ePlotBound) {
  RecVar eRecVar = ModelSpecificVarSpecificTime(TotalArr, eVarName, eTimeDay);
  ApplyPlotBound(TotalArr, eRecVar, eVarName, ePlotBound);
  return eRecVar;
}

std::string GetStrAllOfPlot(VarQuery const &eQuery) {
  int iTime = eQuery.iTime;
  std::string strAll;
  std::string strFile;
  if (eQuery.typeQuery == "direct") {
    //    strFile=DATE_ConvertMjd2mystringFile(eQuery.eTimeDay);
    strFile = DATE_ConvertMjd2mystringFileMilisecond(eQuery.eTimeDay);
  } else {
    std::vector<int> eDate = DATE_ConvertMjd2six(eQuery.eTimeDay);
    int iYear = eDate[0];
    int iMon = eDate[1];
    if (eQuery.typeQuery == "monthly") {
      std::vector<std::string> ListMon{"Jan", "Feb", "Mar", "Apr",
                                       "May", "Jun", "Jul", "Aug",
                                       "Sep", "Oct", "Nov", "Dec"};
      strFile = ListMon[iMon - 1] + std::to_string(iYear);
    }
    if (eQuery.typeQuery == "seasonal") {
      std::vector<std::string> ListSeas{"Winter", "Spring", "Summer", "Autumn"};
      int iSeas = (iMon - 1) / 3;
      strFile = ListSeas[iSeas] + std::to_string(iYear);
    }
    if (eQuery.typeQuery == "yearly") {
      strFile = std::to_string(iYear);
    }
  }
  if (strFile.size() == 0) {
    std::cerr << "Failed to assign strFile\n";
    throw TerminalException{1};
  }
  if (iTime == -1) {
    strAll = strFile;
  } else {
    if (iTime > 10000) {
      std::cerr << "Error in the code\n";
      std::cerr << "iTime is too large\n";
      throw TerminalException{1};
    }
    strAll = StringNumber(iTime, 4) + "_" + strFile;
  }
  if (eQuery.NatureQuery != "instant")
    strAll += "_" + eQuery.NatureQuery;
  return strAll;
}

std::string GetStrPresOfPlot(VarQuery const &eQuery) {
  //  std::string strPres1=DATE_ConvertMjd2mystringPresReduced(eQuery.eTimeDay);
  std::string strPres1 =
      DATE_ConvertMjd2mystringPresReducedMilisecond(eQuery.eTimeDay);
  if (eQuery.NatureQuery == "instant")
    return "at " + strPres1;
  double TimeFrameDay = eQuery.TimeFrameDay;
  std::string strPres2 =
      DATE_ConvertMjd2mystringPresReduced(eQuery.eTimeDay + TimeFrameDay);
  if (eQuery.NatureQuery == "swathMax")
    return "max from " + strPres1 + " to " + strPres2;
  if (eQuery.NatureQuery == "swathMin")
    return "min from " + strPres1 + " to " + strPres2;
  if (eQuery.NatureQuery == "average")
    return "avg. from " + strPres1 + " to " + strPres2;
  std::cerr << "Failed to find NatureQuery in list of available options\n";
  std::cerr << "eQuery.NatureQuery=" << eQuery.NatureQuery << "\n";
  throw TerminalException{1};
}

RecVar ModelSpecificVarSpecificTimeGeneral(TotalArrGetData const &TotalArr,
                                           std::string const &eVarName,
                                           VarQuery const &eQuery,
                                           PlotBound const &ePlotBound) {
  //
  // Check correctness
  //
  std::vector<std::string> ListAllow{"instant", "average", "swathMax",
                                     "swathMin"};
  if (std::find(ListAllow.begin(), ListAllow.end(), eQuery.NatureQuery) ==
      ListAllow.end()) {
    std::cerr << "We failed to find NatureQuery=" << eQuery.NatureQuery << "\n";
    std::cerr << "List of allowed queries:\n";
    for (auto &eStr : ListAllow)
      std::cerr << "  eStr=" << eStr << "\n";
    std::cerr << "Error occurs in ModelSpecificVarSpecificTimeGeneral\n";
    throw TerminalException{1};
  }
  std::string strPres = DATE_ConvertMjd2mystringPres(eQuery.eTimeDay);
  std::cerr << "Query ModelSpecificVarSpecificTimeGeneral NatureQuery="
            << eQuery.NatureQuery << " date=" << strPres
            << " VarName=" << eVarName << "\n";
  //
  // Reading array and doing operations "instant", "average", etc.
  //
  RecVar eRecVar;
  if (eQuery.NatureQuery == "instant") {
    eRecVar = ModelSpecificVarSpecificTimeBound(TotalArr, eVarName,
                                                eQuery.eTimeDay, ePlotBound);
  } else {
    double eTimeDay = eQuery.eTimeDay;
    double TimeFrameDay = eQuery.TimeFrameDay;
    std::vector<int> ListRelITime =
        GetIntervalListITime(TotalArr.eArr, eTimeDay, TimeFrameDay);
    int nbTimeRel = ListRelITime.size();
    if (nbTimeRel == 0) {
      std::cerr << "  NatureQuery = " << eQuery.NatureQuery << "\n";
      std::cerr
          << "  We have nbTimeRel = 0 therefore we cannot process rhe query\n";
      std::cerr << "  eTimeDay=" << eTimeDay << "\n";
      std::cerr << "  TimeFrameDay=" << TimeFrameDay << "\n";
      throw TerminalException{1};
    }
    RecVar RecVarTrivial = RetrieveTrivialRecVar(eVarName);
    MyMatrix<double> F, U, V;
    Eigen::Tensor<double, 3> Tens3;
    Eigen::Tensor<double, 3> Uthree;
    Eigen::Tensor<double, 3> Vthree;
    std::cerr << "nbTimeRel = " << nbTimeRel << "\n";
    for (int iTimeRel = 0; iTimeRel < nbTimeRel; iTimeRel++) {
      int iTime = ListRelITime[iTimeRel];
      double eTimeDayB = ARR_GetTime(TotalArr.eArr, iTime);
      std::string strPres = DATE_ConvertMjd2mystringPres(eTimeDayB);
      eRecVar = ModelSpecificVarSpecificTimeBound(TotalArr, eVarName, eTimeDayB,
                                                  ePlotBound);
      std::cerr << "iTimeRel=" << iTimeRel << "/" << nbTimeRel
                << " date=" << strPres << "\n";
      //      std::cerr << "iTimeRel=" << iTimeRel << " iTime=" << iTime << "
      //      VarNature=" << eRecVar.RecS.VarNature << "\n";
      if (iTimeRel == 0) {
        if (RecVarTrivial.RecS.VarNature == "rho")
          F = eRecVar.F;
        if (RecVarTrivial.RecS.VarNature == "uv") {
          U = eRecVar.U;
          V = eRecVar.V;
          F = eRecVar.F;
        }
        if (RecVarTrivial.RecS.VarNature == "3Drho")
          Tens3 = eRecVar.Tens3;
        if (RecVarTrivial.RecS.VarNature == "3Duv") {
          Uthree = eRecVar.Uthree;
          Vthree = eRecVar.Vthree;
          Tens3 = eRecVar.Tens3;
        }
      } else {
        if (eQuery.NatureQuery == "average") {
          if (RecVarTrivial.RecS.VarNature == "rho")
            F += eRecVar.F;
          if (RecVarTrivial.RecS.VarNature == "uv") {
            U += eRecVar.U;
            V += eRecVar.V;
            F += eRecVar.F;
          }
          if (RecVarTrivial.RecS.VarNature == "3Drho")
            Tens3 += eRecVar.Tens3;
          if (RecVarTrivial.RecS.VarNature == "3Duv") {
            Uthree += eRecVar.Uthree;
            Vthree += eRecVar.Vthree;
            Tens3 += eRecVar.Tens3;
          }
        }
        if (eQuery.NatureQuery == "swathMax") {
          if (RecVarTrivial.RecS.VarNature == "rho")
            F = F.cwiseMax(eRecVar.F);
          if (RecVarTrivial.RecS.VarNature == "3Drho")
            Tens3 = Tens3.cwiseMax(eRecVar.Tens3);
          if (RecVarTrivial.RecS.VarNature == "uv" ||
              RecVarTrivial.RecS.VarNature == "3Duv") {
            std::cerr << "swathMax for uv does not have any sense\n";
            throw TerminalException{1};
          }
        }
        if (eQuery.NatureQuery == "swathMin") {
          if (RecVarTrivial.RecS.VarNature == "rho")
            F = F.cwiseMin(eRecVar.F);
          if (RecVarTrivial.RecS.VarNature == "3Drho")
            Tens3 = Tens3.cwiseMin(eRecVar.Tens3);
          if (RecVarTrivial.RecS.VarNature == "uv" ||
              RecVarTrivial.RecS.VarNature == "3Duv") {
            std::cerr << "swathMin for uv does not have any sense\n";
            throw TerminalException{1};
          }
        }
      }
    }
    if (eQuery.NatureQuery == "average") {
      double fact = 1 / static_cast<double>(nbTimeRel);
      if (RecVarTrivial.RecS.VarNature == "rho") {
        F *= fact;
      }
      if (RecVarTrivial.RecS.VarNature == "uv") {
        U *= fact;
        V *= fact;
        F *= fact;
      }
      if (RecVarTrivial.RecS.VarNature == "3Drho") {
        // does not compile
        // Tens3 /= double(nbTimeRel);
        auto LDim = Tens3.dimensions();
        int dim0 = LDim[0];
        int dim1 = LDim[1];
        int dim2 = LDim[2];
        for (int i0 = 0; i0 < dim0; i0++)
          for (int i1 = 0; i1 < dim1; i1++)
            for (int i2 = 0; i2 < dim2; i2++)
              Tens3(i0, i1, i2) *= fact;
      }
      if (RecVarTrivial.RecS.VarNature == "3Duv") {
        // does not compile
        // Tens3 /= double(nbTimeRel);
        auto LDim = Tens3.dimensions();
        int dim0 = LDim[0];
        int dim1 = LDim[1];
        int dim2 = LDim[2];
        for (int i0 = 0; i0 < dim0; i0++)
          for (int i1 = 0; i1 < dim1; i1++)
            for (int i2 = 0; i2 < dim2; i2++) {
              Uthree(i0, i1, i2) *= fact;
              Vthree(i0, i1, i2) *= fact;
              Tens3(i0, i1, i2) *= fact;
            }
      }
    }
    if (RecVarTrivial.RecS.VarNature == "rho") {
      eRecVar.F = F;
    }
    if (RecVarTrivial.RecS.VarNature == "3Drho") {
      eRecVar.Tens3 = Tens3;
    }
    if (RecVarTrivial.RecS.VarNature == "uv") {
      eRecVar.U = U;
      eRecVar.V = V;
      eRecVar.F = F;
    }
    if (RecVarTrivial.RecS.VarNature == "3Duv") {
      eRecVar.Uthree = Uthree;
      eRecVar.Vthree = Vthree;
      eRecVar.Tens3 = Tens3;
    }
  }
  //  std::cerr << "End of ModelSpecificVarSpecificTimeGeneral\n";
  std::string strAll = GetStrAllOfPlot(eQuery);
  std::string strPresPlot = GetStrPresOfPlot(eQuery);
  //  std::cerr << "strAll=" << strAll << "\n";
  eRecVar.RecS.strAll = strAll;
  eRecVar.RecS.strPres = strPresPlot;
  ApplyPlotBound(TotalArr, eRecVar, eVarName, ePlotBound);
  //  std::cerr << "After ApplyPlotBound\n";
  return eRecVar;
}

PairRecVar ModelPairSpecificVarSpecificTimeGeneral(
    TotalArrGetData const &TotalArr1, TotalArrGetData const &TotalArr2,
    std::string const &eVarName, VarQuery const &eQuery,
    PlotBound const &ePlotBound) {
  //  std::cerr << "NatureQuery=" << eQuery.NatureQuery << "\n";
  std::vector<std::string> ListAllow{"instant", "average", "swathMax",
                                     "swathMin"};
  RecVar eRecVar1, eRecVar2;
  //  std::cerr << "Begin ModelPairSpecificVarSpecificTimeGeneral\n";
  if (std::find(ListAllow.begin(), ListAllow.end(), eQuery.NatureQuery) !=
      ListAllow.end()) {
    eRecVar1 = ModelSpecificVarSpecificTimeGeneral(TotalArr1, eVarName, eQuery,
                                                   ePlotBound);
    eRecVar2 = ModelSpecificVarSpecificTimeGeneral(TotalArr2, eVarName, eQuery,
                                                   ePlotBound);
    ApplyPlotBoundPair(TotalArr1, TotalArr2, eRecVar1, eRecVar2, eVarName,
                       ePlotBound);
    //    std::cerr << "Begin ModelPairSpecificVarSpecificTimeGeneral, case
    //    1\n";
  } else {
    //    std::cerr << "Begin ModelPairSpecificVarSpecificTimeGeneral, case
    //    2\n";
    std::vector<std::string> ListAllowB{"MaxDiff", "MinDiff"};
    if (std::find(ListAllowB.begin(), ListAllowB.end(), eQuery.NatureQuery) ==
        ListAllowB.end()) {
      std::cerr
          << "Only two difference operator are allowed: MaxDiff and MinDiff\n";
      std::cerr << "NatureQuery=" << eQuery.NatureQuery << "\n";
    }
    RecVar RecVarTrivial = RetrieveTrivialRecVar(eVarName);
    //
    double eTimeDay = eQuery.eTimeDay;
    double TimeFrameDay = eQuery.TimeFrameDay;
    std::vector<int> ListRelITime1 =
        GetIntervalListITime(TotalArr1.eArr, eTimeDay, TimeFrameDay);
    std::vector<int> ListRelITime2 =
        GetIntervalListITime(TotalArr2.eArr, eTimeDay, TimeFrameDay);
    int nbTimeRel1 = ListRelITime1.size();
    int nbTimeRel2 = ListRelITime2.size();
    if (nbTimeRel1 == 0 || nbTimeRel2 == 0) {
      std::cerr << "We have nbTimeRel1 = " << nbTimeRel1 << "\n";
      std::cerr << "We have nbTimeRel2 = " << nbTimeRel2 << "\n";
      std::cerr << "We need both to be non-zero\n";
      throw TerminalException{1};
    }
    double DeltaTime1 =
        (ARR_GetTime(TotalArr1.eArr, ListRelITime1[nbTimeRel1 - 1]) -
         ARR_GetTime(TotalArr1.eArr, ListRelITime1[0])) /
        static_cast<double>(nbTimeRel1 - 1);
    double DeltaTime2 =
        (ARR_GetTime(TotalArr2.eArr, ListRelITime2[nbTimeRel2 - 1]) -
         ARR_GetTime(TotalArr2.eArr, ListRelITime2[0])) /
        static_cast<double>(nbTimeRel2 - 1);
    double DeltaTime = std::min(DeltaTime1, DeltaTime2);
    double FirstTime = eTimeDay;
    double LastTime = eTimeDay + TimeFrameDay;
    std::vector<double> ListRelTimeDay =
        GetIntervalFLD(FirstTime, LastTime, DeltaTime);
    int nbTimeRel = ListRelTimeDay.size();
    MyMatrix<double> diffF;
    if (RecVarTrivial.RecS.VarNature == "uv") {
      std::cerr << "Cannot do MaxDiff or MinDiff for variables of type uv\n";
      throw TerminalException{1};
    }
    for (int iTimeRel = 0; iTimeRel < nbTimeRel; iTimeRel++) {
      double eTimeDayB = ListRelTimeDay[iTimeRel];
      eRecVar1 = ModelSpecificVarSpecificTimeBound(TotalArr1, eVarName,
                                                   eTimeDayB, ePlotBound);
      eRecVar2 = ModelSpecificVarSpecificTimeBound(TotalArr2, eVarName,
                                                   eTimeDayB, ePlotBound);
      if (iTimeRel == 0) {
        if (RecVarTrivial.RecS.VarNature == "rho") {
          diffF = eRecVar1.F - eRecVar2.F;
        }
      } else {
        if (eQuery.NatureQuery == "MaxDiff") {
          if (RecVarTrivial.RecS.VarNature == "rho")
            diffF = diffF.cwiseMax(eRecVar1.F - eRecVar2.F);
        }
        if (eQuery.NatureQuery == "MinDiff") {
          if (RecVarTrivial.RecS.VarNature == "rho")
            diffF = diffF.cwiseMin(eRecVar1.F - eRecVar2.F);
        }
      }
    }
    if (RecVarTrivial.RecS.VarNature == "rho") {
      eRecVar1.F = diffF;
      eRecVar2.F.fill(static_cast<double>(0));
    }
    ApplyPlotBound(TotalArr1, eRecVar1, eVarName, ePlotBound);
    ApplyPlotBound(TotalArr2, eRecVar2, eVarName, ePlotBound);
    ApplyPlotBoundPair(TotalArr1, TotalArr2, eRecVar1, eRecVar2, eVarName,
                       ePlotBound);
  }
  std::string strAll = GetStrAllOfPlot(eQuery);
  eRecVar1.RecS.strAll = strAll;
  eRecVar2.RecS.strAll = strAll;
  return {eRecVar1, eRecVar2};
}

RecVar GetTrivialArrayPlot(GridArray const &GrdArr) {
  int eta_rho = GrdArr.GrdArrRho.LON.rows();
  int xi_rho = GrdArr.GrdArrRho.LON.cols();
  MyMatrix<double> F(eta_rho, xi_rho);
  double MaxLon = GrdArr.GrdArrRho.LON.maxCoeff();
  double MinLon = GrdArr.GrdArrRho.LON.minCoeff();
  for (int i = 0; i < eta_rho; i++)
    for (int j = 0; j < xi_rho; j++) {
      double eLon = GrdArr.GrdArrRho.LON(i, j);
      double eVal = (eLon - MinLon) / (MaxLon - MinLon);
      F(i, j) = eVal;
    }
  RecVar eRecVar;
  eRecVar.RecS.strAll = "trivialarrayplot";
  eRecVar.RecS.VarName1 = "Track";
  eRecVar.RecS.VarName2 = "Track";
  eRecVar.RecS.minval = 2;
  eRecVar.RecS.maxval = 3;
  eRecVar.RecS.Unit = "nondim.";
  eRecVar.F = F;
  return eRecVar;
}

void Set_iTime_eTimeDay(RecVar &eRecVar, int iTime, double eTimeDay) {
  eRecVar.RecS.iTime = iTime;
  eRecVar.RecS.eTimeDay = eTimeDay;
  eRecVar.RecS.strPres = DATE_ConvertMjd2mystringPres(eTimeDay);
  eRecVar.RecS.strFile = DATE_ConvertMjd2mystringFile(eTimeDay);
}

// clang-format off
#endif  // SRC_OCEAN_MODEL_DATA_LOADING_H_
// clang-format on
