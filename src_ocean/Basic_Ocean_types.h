#ifndef SRC_OCEAN_BASIC_OCEAN_TYPES_H_
#define SRC_OCEAN_BASIC_OCEAN_TYPES_H_

#include "MAT_Matrix.h"
#include "MAT_Tensor.h"
#include "Timings.h"
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <unordered_map>

struct PairLL {
  double eLon;
  double eLat;
};

bool operator<(PairLL const &x, PairLL const &y) {
  if (x.eLon < y.eLon)
    return true;
  if (x.eLon > y.eLon)
    return false;
  return x.eLat < y.eLat;
}

struct RecSymbolic {
  double eTimeDay;
  int iTime;
  std::string strPres;
  std::string strFile;
  std::string strAll;
  //
  std::string FullVarName;
  std::string VarName1;
  std::string VarName2;
  std::string CFshortName;
  double minval;
  double maxval;
  double mindiff;
  double maxdiff;
  std::string Unit;
  std::string VarNature;
  std::string nameU, nameV;
  std::string strTime_ROMS;
  std::optional<std::string> varName_ROMS;
  std::optional<std::string> varName_ROMS_U, varName_ROMS_V;
  std::optional<std::string> varName_GRIB;
  std::optional<std::string> varName_GRIB_U, varName_GRIB_V;
};

struct NEMO_vars {
  std::vector<std::string> List2D_vars;
  std::vector<std::string> List3D_vars;
};

struct RecVar {
  RecSymbolic RecS;
  MyMatrix<double> U;
  MyMatrix<double> V;
  MyMatrix<double> F;
  Eigen::Tensor<double, 3> Uthree;
  Eigen::Tensor<double, 3> Vthree;
  Eigen::Tensor<double, 3> Tens3;
};

struct CoordGridArrayFD {
  MyMatrix<uint8_t> MSK;
  MyMatrix<double> LON, LAT, ANG;
  std::optional<MyMatrix<double>> DEP;
  //  For ROMS model pm and pn are needed
  MyMatrix<double> pm, pn;
  int nbWet;
  std::vector<int> Idx, Jdx;
};

const MyMatrix<double> &GetDEP(CoordGridArrayFD const &arr) {
  if (arr.DEP)
    return *arr.DEP;
  std::cerr << "The bathymetry has not been assigned\n";
  throw TerminalException{1};
}

struct ARVDtyp {
  bool IsAssigned;
  std::string ModelName;
  bool Zcoordinate;
  // N: number of vertical levels
  int N;
  // ListZ_r / ListZ_w are used for Z-coordinates models
  MyVector<double> ListZ_r;
  MyVector<double> ListZ_w;
  //  Sometimes, we have a vertical mask for Z-coordinates models
  std::optional<Eigen::Tensor<uint8_t, 3>> TensMSKvert;
  double Tcline;
  double hc;
  double theta_s;
  double theta_b;
  int Vtransform;
  int Vstretching;
  MyVector<double> Cs_r;
  MyVector<double> Cs_w;
  MyVector<double> sc_r;
  MyVector<double> sc_w;
};

//  IsSpherical = T for coordinates done in LON/LAT ( deg)
//  if = F then coordinates are different (maybe meter)
//  but we use the same names of LON/LAT, just meaning is different.
struct GridArray {
  std::string ModelName;
  int IsFE;
  bool IsSpherical;
  CoordGridArrayFD GrdArrRho, GrdArrU, GrdArrV, GrdArrPsi;
  MyMatrix<int> INE;
  MyVector<int> IOBP;
  bool L_IndexSelect;
  std::vector<int> I_IndexSelect;
  ARVDtyp ARVD;
};

struct GRIB_MessageInfo {
  std::string shortName;
  std::string cfVarName;
  std::string cfVarNameECMWF;
  std::string name;
  std::string units;
  int idx;
  double time;
  double timeStart;
  int stepRange;
  std::string FileName;
};

struct TripleXYZ {
  MyMatrix<double> X;
  MyMatrix<double> Y;
  MyMatrix<double> Z;
};

struct AnalyticalAlgorithm {
  std::vector<std::string> ListNameVariables;
  std::vector<double> ListConstantValuesRho;
  std::vector<double> ListConstantValuesU;
  std::vector<double> ListConstantValuesV;
};

struct MeasurementSingPos {
  double x;
  double y;
  double z;
  double dep;
  double year;
  double day;
};

struct MeasurementData {
  std::vector<MeasurementSingPos> l_pos;
  std::vector<double> l_data;
};

struct CompleteMeasurementData {
  std::unordered_map<std::string, MeasurementData> arr_meas;
};

struct ArrayHistory {
  std::string KindArchive;
  int nbFile, nbTime;
  double FirstTime, LastTime;
  std::string FirstTimeStr, LastTimeStr;
  std::vector<std::string> ListFileNames;
  std::vector<std::vector<GRIB_MessageInfo>> ListListMessages;
  std::unordered_map<std::string, std::vector<int>> MatchingByVariable;
  std::vector<GRIB_MessageInfo> ListAllMessage;
  std::vector<std::string> RawVarNames;
  std::vector<int> ListITime;
  std::vector<double> ListStartTime;
  std::vector<double> ListEndTime;
  std::vector<std::vector<int>> ListListIMesg;
  std::vector<int> ListIStartTime;
  std::vector<int> ListIFile;
  std::vector<int> ListIRec;
  std::vector<double> ListTime;
  std::string TimeSteppingInfo;
  std::string HisPrefix;
  std::map<std::string,
           std::vector<std::pair<double, std::vector<GRIB_MessageInfo>>>>
      FullOrganizedInfo;
  //  For NEMO, data is distributed into many prefix
  std::map<std::string, std::string> NEMO_vars_to_postfix;
  double SeparationTime;
  int nbDigit;
  int nbRecBegin;
  int nbRecMiddle;
  bool AppendVarName;
  std::string eModelName;
};

struct TotalArrGetData {
  GridArray GrdArr;
  ArrayHistory eArr;
};

// Same as in ROMS, ranges: Hz (0..N-1), z_w (0..N), z_r(0..N-1)
struct VerticalInfo {
  MyVector<double> Hz;
  MyVector<double> z_w;
  MyVector<double> z_r;
};

struct PlotBound {
  bool VariableRange;
  std::vector<std::string> BoundSingle_var;
  std::vector<double> BoundSingle_min;
  std::vector<double> BoundSingle_max;
  std::vector<std::string> BoundDiff_var;
  std::vector<double> BoundDiff_min;
  std::vector<double> BoundDiff_max;
};

struct PairRecVar {
  RecVar RecVar1;
  RecVar RecVar2;
};

struct QuadArray {
  double MinLon;
  double MaxLon;
  double MinLat;
  double MaxLat;
};

struct SingleArrayRegionAveraging {
  std::vector<std::vector<std::pair<int, int>>> ListListEtaXi;
};

struct SingleArrayInterpolation {
  int eta_out, xi_out;
  int eta_in, xi_in;
  GridArray GrdArrOut;
  ARVDtyp ARVDin;
  MySparseMatrix<double> SpMat;
  MyMatrix<double> DEPinInterp;
};

struct SingleArrayInterpolationGen {
  SingleArrayInterpolation e_arr;
  std::vector<SingleArrayInterpolation> l_arr;
};

struct TransectInformation {
  std::vector<PairLL> ListPairLL;
  MyVector<double> ListDimVar;
  std::vector<SingleArrayInterpolationGen> ListRec;
};

struct TransectInformation_3D {
  std::vector<PairLL> ListPairLL;
  MyVector<double> ListDimVar;
  SingleArrayInterpolationGen eRecInterp;
  //
  MyMatrix<double> matCoeffM;
  MyMatrix<double> matCoeffP;
  MyMatrix<int> matIdxM;
  MyMatrix<int> matIdxP;
  MyMatrix<uint8_t> MSK;
  int Ntotal;
  MyVector<double> ListVertPos;
  std::vector<int> FirstWetIndex;
  //
  double normU, normV;
};

#endif  // SRC_OCEAN_BASIC_OCEAN_TYPES_H_
