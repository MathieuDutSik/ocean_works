#ifndef SRC_OCEAN_STATISTICS_H_
#define SRC_OCEAN_STATISTICS_H_

#include "Basic_string.h"
#include "MAT_Matrix.h"
#include "MAT_Tensor.h"
#include <algorithm>
#include <limits>
#include <string>
#include <utility>
#include <vector>

struct T_stat {
  int nbMeas;
  double MaxMeas;
  double MinMeas;
  double MaxModel;
  double MinModel;
  double MeanMeas;
  double MeanModel;
  double MeanError;
  double AbsoluteError;
  double RMSE;
  double CenteredRMSE;
  double Correlation;
  double ScatterIndex;
  double CenteredScatterIndex;
  double Slope;
};

struct PairMM {
  double Meas;
  double Model;
};

struct T_statString {
  std::string strMaxMeas;
  std::string strMinMeas;
  std::string strMaxModel;
  std::string strMinModel;
  std::string strMeanMeas;
  std::string strMeanModel;
  std::string strMeanError;
  std::string strAbsoluteError;
  std::string strRMSE;
  std::string strCenteredRMSE;
  std::string strCorrelation;
  std::string strScatterIndex;
  std::string strCenteredScatterIndex;
  std::string strSlope;
  std::string strNature = "ME    AE    RMSE CRMSE  CORR   SCI   CSCI";
  std::string str;
};

T_stat ComputeStatistics_Pair(std::vector<PairMM> const &eVect) {
  //  std::cerr << "step 1\n";
  T_stat eStat;
  int nbMeas = 0;
  double SumAbs = 0;
  double SumSqr = 0;
  double eSum1 = 0;
  double eSum2 = 0;
  double eSum11 = 0;
  double eSum12 = 0;
  double eSum22 = 0;
  double MaxMeas = -10 ^ (31);
  double MaxModel = -10 ^ (31);
  double MinMeas = 10 ^ (31);
  double MinModel = 10 ^ (31);
  //  std::cerr << "step 2\n";
  for (auto &ePair : eVect) {
    nbMeas++;
    double eMeas = ePair.Meas;
    double eModel = ePair.Model;
    //    std::cerr << "eMeas=" << eMeas << " eModel=" << eModel << "\n";
    MaxMeas = std::max(MaxMeas, eMeas);
    MaxModel = std::max(MaxModel, eModel);
    MinMeas = std::min(MinMeas, eMeas);
    MinModel = std::min(MinModel, eModel);
    eSum1 += eMeas;
    eSum2 += eModel;
    eSum11 += eMeas * eMeas;
    eSum12 += eMeas * eModel;
    eSum22 += eModel * eModel;
    SumAbs += fabs(eMeas - eModel);
    double eDiff = eMeas - eModel;
    SumSqr += eDiff * eDiff;
  }
  //  std::cerr << "step 3\n";
  double eME = (eSum2 - eSum1) / double(nbMeas);
  double eRMSE = sqrt(SumSqr / double(nbMeas));
  double eCentRMSE = sqrt(eRMSE * eRMSE - eME * eME);
  double eAE = SumAbs / double(nbMeas);
  double avgSum1 = eSum1 / double(nbMeas);
  double avgSum2 = eSum2 / double(nbMeas);
  double avgSum11 = eSum11 / double(nbMeas);
  double avgSum12 = eSum12 / double(nbMeas);
  double avgSum22 = eSum22 / double(nbMeas);
  double eProd11 = avgSum11 - avgSum1 * avgSum1;
  double eProd12 = avgSum12 - avgSum1 * avgSum2;
  double eProd22 = avgSum22 - avgSum2 * avgSum2;
  double TheCorr = eProd12 / sqrt(eProd11 * eProd22);
  double eScat = eRMSE / avgSum1;
  double eCentScat = eCentRMSE / avgSum1;
  double eSlope = eSum12 / eSum11;
  //  std::cerr << "step 4\n";
  //
  eStat.nbMeas = nbMeas;
  eStat.MaxMeas = MaxMeas;
  eStat.MinMeas = MinMeas;
  eStat.MaxModel = MaxModel;
  eStat.MinModel = MinModel;
  eStat.MeanMeas = avgSum1;
  eStat.MeanModel = avgSum2;
  eStat.MeanError = eME;
  eStat.AbsoluteError = eAE;
  eStat.RMSE = eRMSE;
  eStat.CenteredRMSE = eCentRMSE;
  eStat.Correlation = TheCorr;
  eStat.ScatterIndex = eScat;
  eStat.CenteredScatterIndex = eCentScat;
  eStat.Slope = eSlope;
  /*
  std::cerr << "step 5\n";
  std::cerr << "MaxMeas=" << MaxMeas << "\n";
  std::cerr << "MinMeas=" << MinMeas << "\n";
  std::cerr << "MaxModel=" << MaxModel << "\n";
  std::cerr << "MinModel=" << MinModel << "\n";

  std::cerr << "avgSum1=" << avgSum1 << "\n";
  std::cerr << "avgSum2=" << avgSum2 << "\n";
  std::cerr << "eME=" << eME << "\n";
  std::cerr << "eAE=" << eAE << "\n";
  std::cerr << "eRMSE=" << eRMSE << "\n";
  std::cerr << "eCentRMSE=" << eCentRMSE << "\n";
  std::cerr << "TheCorr=" << TheCorr << "\n";
  std::cerr << "eScat=" << eScat << "\n";
  std::cerr << "eCentScat=" << eCentScat << "\n";
  std::cerr << "eSlope=" << eSlope << "\n";*/
  //
  return eStat;
}

T_statString ComputeStatisticString_from_Statistics(T_stat const &eStat,
                                                    std::string const &opt) {
  auto fctstring = [&opt](double const &x) -> std::string {
    if (opt == "double") {
      return DoubleToString(x);
    }
    if (opt == "4dot2f") {
      return DoubleTo4dot2f(x);
    }
    std::cerr << "Allowed options are double or 4dot2f\n";
    std::cerr << "opt=" << opt << "\n";
    std::cerr << "So, wrong statement here\n";
    throw TerminalException{1};
  };
  T_statString eStatStr;
  eStatStr.strMaxMeas = fctstring(eStat.MaxMeas);
  //  std::cerr << "  print step 1\n";
  eStatStr.strMinMeas = fctstring(eStat.MinMeas);
  //  std::cerr << "  print step 2\n";
  eStatStr.strMaxModel = fctstring(eStat.MaxModel);
  //  std::cerr << "  print step 3\n";
  eStatStr.strMinModel = fctstring(eStat.MinModel);
  //  std::cerr << "  print step 4\n";
  eStatStr.strMeanMeas = fctstring(eStat.MeanMeas);
  //  std::cerr << "  print step 5\n";
  eStatStr.strMeanModel = fctstring(eStat.MeanModel);
  //  std::cerr << "  print step 6\n";
  eStatStr.strMeanError = fctstring(eStat.MeanError);
  //  std::cerr << "  print step 7\n";
  eStatStr.strAbsoluteError = fctstring(eStat.AbsoluteError);
  //  std::cerr << "  print step 8\n";
  eStatStr.strRMSE = fctstring(eStat.RMSE);
  //  std::cerr << "  print step 9\n";
  eStatStr.strCenteredRMSE = fctstring(eStat.CenteredRMSE);
  //  std::cerr << "  print step 10\n";
  eStatStr.strCorrelation = fctstring(eStat.Correlation);
  //  std::cerr << "  print step 11\n";
  eStatStr.strScatterIndex = fctstring(eStat.ScatterIndex);
  //  std::cerr << "  print step 12\n";
  eStatStr.strCenteredScatterIndex = fctstring(eStat.CenteredScatterIndex);
  //  std::cerr << "  print step 13\n";
  eStatStr.strSlope = fctstring(eStat.Slope);
  //  std::cerr << "step 6\n";
  //
  eStatStr.str = eStatStr.strMeanError + " " + eStatStr.strAbsoluteError + " " +
                 eStatStr.strRMSE + " " + eStatStr.strCenteredRMSE + " " +
                 eStatStr.strCorrelation + " " + eStatStr.strScatterIndex +
                 " " + eStatStr.strCenteredScatterIndex;
  return eStatStr;
}

template <typename T1, typename T2>
void ComputeStatisticCheckSizes(T1 const &v1, T2 const &v2) {
  size_t siz1 = v1.size();
  size_t siz2 = v2.size();
  if (siz1 != siz2) {
    std::cerr << "|ListMeas|=" << siz1 << " |ListModel|=" << siz2 << "\n";
    std::cerr << "Error in ComputeStatistics_vector\n";
    std::cerr << "Discrepancy in number of measurements\n";
    std::cerr << "Please solve the problem\n";
    throw TerminalException{1};
  }
}

template <typename Fset> T_stat ComputeStatistics_F(size_t nbEnt, Fset f_set) {
  std::vector<PairMM> ListPair(nbEnt);
  for (size_t iEnt = 0; iEnt < nbEnt; iEnt++)
    ListPair[iEnt] = f_set(iEnt);
  return ComputeStatistics_Pair(ListPair);
}

T_stat ComputeStatistics_vector(std::vector<double> const &ListMeas,
                                std::vector<double> const &ListModel) {
  ComputeStatisticCheckSizes(ListMeas, ListModel);
  auto f_set = [&](size_t iEnt) -> PairMM {
    return {ListMeas[iEnt], ListModel[iEnt]};
  };
  return ComputeStatistics_F(ListMeas.size(), f_set);
}

T_stat ComputeStatistics_MyVector(MyVector<double> const &ListMeas,
                                  MyVector<double> const &ListModel) {
  ComputeStatisticCheckSizes(ListMeas, ListModel);
  auto f_set = [&](size_t iEnt) -> PairMM {
    return {ListMeas(iEnt), ListModel(iEnt)};
  };
  return ComputeStatistics_F(ListMeas.size(), f_set);
}

T_stat ComputeStatistics_stdpair(
    std::vector<std::pair<double, double>> const &ListPairMM) {
  std::vector<PairMM> ListPair;
  for (auto &ePair : ListPairMM)
    ListPair.push_back({ePair.first, ePair.second});
  std::cerr << "After ListPair creation\n";
  return ComputeStatistics_Pair(ListPair);
}

void Print_Down_Statistics(std::ostream &os, std::string const &eName,
                           T_stat const &eStat) {
  int nbMeas = eStat.nbMeas;
  os << "        eName=" << eName << " nbMeas=" << nbMeas << "\n";
  if (nbMeas > 0) {
    os << "      MaxMeas=" << eStat.MaxMeas << "\n";
    os << "      MinMeas=" << eStat.MinMeas << "\n";
    os << "     MaxModel=" << eStat.MaxModel << "\n";
    os << "     MinModel=" << eStat.MinModel << "\n";
    os << "     MeanMeas=" << eStat.MeanMeas << "\n";
    os << "    MeanModel=" << eStat.MeanModel << "\n";
    os << "    MeanError=" << eStat.MeanError << "\n";
    os << "AbsoluteError=" << eStat.AbsoluteError << "\n";
    os << "         RMSE=" << eStat.RMSE << "\n";
    os << " CenteredRMSE=" << eStat.CenteredRMSE << "\n";
    os << "  Correlation=" << eStat.Correlation << "\n";
    os << " ScatterIndex=" << eStat.ScatterIndex << "\n";
    os << "  CenteredSci=" << eStat.CenteredScatterIndex << "\n";
  }
}

void PrintMMA_FCT(MyMatrix<double> const &F, MyMatrix<uint8_t> const &MSK,
                  std::string const &VarName, std::string const &UnitName) {
  int eta = F.rows();
  int xi = F.cols();
  double minval = 0, maxval = 0;
  double sum = 0;
  int nb = 0;
  bool IsFirst = true;
  int nWet = 0;
  for (int i = 0; i < eta; i++)
    for (int j = 0; j < xi; j++)
      if (MSK(i, j) == 1) {
        nWet++;
        double eVal = F(i, j);
        sum += eVal;
        nb++;
        if (IsFirst) {
          IsFirst = false;
          maxval = eVal;
          minval = eVal;
        } else {
          if (eVal > maxval)
            maxval = eVal;
          if (eVal < minval)
            minval = eVal;
        }
      }
  double eMean = sum / double(nb);
  std::cerr << "  " << VarName << " min=" << minval << " max=" << maxval
            << " avg=" << eMean << " " << UnitName << " nWet=" << nWet
            << " eta=" << eta << " xi=" << xi << "\n";
}

void IdentifyLinearInterpolationParts(MyVector<double> &ListVal,
                                      MyVector<double> const &ListTime,
                                      double const &MissingValue,
                                      double const &minAbsSlope,
                                      double const &MaxChangeV) {
  int len = ListVal.size();
  std::vector<int> ListStatus(len, 1);
  for (int i = 2; i < len; i++) {
    double slope =
        (ListVal(i) - ListVal(i - 2)) / (ListTime(i) - ListTime(i - 2));
    if (fabs(slope) > minAbsSlope) {
      double interpolVal =
          ListVal(i - 2) + (ListTime(i - 1) - ListTime(i - 2)) * slope;
      if (fabs(interpolVal - ListVal(i - 1)) < MaxChangeV)
        ListStatus[i - 1] = 0;
    }
  }
  for (int i = 0; i < len; i++)
    if (ListStatus[i] == 0)
      ListVal(i) = MissingValue;
}

std::vector<double> GetMinMaxAvg(Eigen::Tensor<double, 3> const &eTens) {
  auto LDim = eTens.dimensions();
  int dim0 = LDim[0];
  int dim1 = LDim[1];
  int dim2 = LDim[2];
  double maxval = std::numeric_limits<double>::min();
  double minval = std::numeric_limits<double>::max();
  double sumval = 0;
  for (int i0 = 0; i0 < dim0; i0++)
    for (int i1 = 0; i1 < dim1; i1++)
      for (int i2 = 0; i2 < dim2; i2++) {
        double val = eTens(i0, i1, i2);
        if (val < minval)
          minval = val;
        if (val > maxval)
          maxval = val;
        sumval += val;
      }
  double avgval = sumval / (dim0 * dim1 * dim2);
  return {minval, maxval, avgval};
}

// clang-format off
#endif // SRC_OCEAN_STATISTICS_H_
// clang-format on
