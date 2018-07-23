#ifndef BASIC_OCEAN_TYPES_DEFINE
#define BASIC_OCEAN_TYPES_DEFINE

#include "MAT_Matrix.h"
#include "MAT_Tensor.h"

struct PairLL {
  double eLon;
  double eLat;
};


bool operator<(PairLL const& x, PairLL const& y)
{
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
  std::string varName_ROMS;
  std::string varName_ROMS_U, varName_ROMS_V;
  std::string varName_GRIB;
  std::string varName_GRIB_U, varName_GRIB_V;
};



struct RecVar {
  RecSymbolic RecS;
  MyMatrix<double> U;
  MyMatrix<double> V;
  MyMatrix<double> F;
  Eigen::Tensor<double,3> Uthree;
  Eigen::Tensor<double,3> Vthree;
  Eigen::Tensor<double,3> Tens3;
};

struct CoordGridArrayFD {
  MyMatrix<int> MSK;
  MyMatrix<double> LON, LAT, DEP, ANG;
  int nbWet;
  bool HaveDEP;
  std::vector<int> Idx, Jdx;
};


struct ARVDtyp {
  bool IsAssigned;
  std::string ModelName;
  bool Zcoordinate;
  int N; // number of vertical levels
  MyVector<double> ListZ_r; // only used for Z-coordinates models
  MyVector<double> ListZ_w; // only used for Z-coordinates models
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

// IsSpherical = T for coordinates done in LON/LAT ( deg)
// if = F then coordinates are different (maybe meter)
// but we use the same names of LON/LAT, just meaning is different.
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



struct ArrayHistory {
  std::string KindArchive;
  int nbFile, nbTime;
  double FirstTime, LastTime;
  std::string FirstTimeStr, LastTimeStr;
  std::vector<std::string> ListFileNames;
  std::vector<std::vector<GRIB_MessageInfo>> ListListMessages;
  std::map<std::string, std::vector<int> > MatchingByVariable;
  std::vector<GRIB_MessageInfo> ListAllMessage;
  std::vector<std::string> RawVarNames;
  std::vector<int> ListITime;
  std::vector<double> ListStartTime;
  std::vector<int> ListIStartTime;
  std::vector<int> ListIFile;
  std::vector<int> ListIRec;
  std::vector<double> ListTime;
  std::string TimeSteppingInfo;
  std::string HisPrefix;
  std::map<std::string, std::vector<std::pair<double, std::vector<GRIB_MessageInfo>>>> FullOrganizedInfo;
  double SeparationTime;
  int nbDigit;
  int nbRecBegin;
  int nbRecMiddle;
  bool AppendVarName;
};




struct VerticalInfo {
  MyVector<double> Hz;  // range is 0..N-1
  MyVector<double> z_w; // range is 0..N
  MyVector<double> z_r; // range is 0..N-1
};

struct TotalArrGetData {
  GridArray GrdArr;
  ArrayHistory eArr;
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


struct SingleArrayInterpolation {
  int eta_out, xi_out;
  int eta_in, xi_in;
  GridArray GrdArrOut;
  ARVDtyp ARVDin;
  MySparseMatrix<double> SpMat;
};


struct TransectInformation {
  std::vector<PairLL> ListPairLL;
  MyVector<double> ListDimVar;
  std::vector<SingleArrayInterpolation> ListRec;
};

struct TransectInformation_3D {
  std::vector<PairLL> ListPairLL;
  MyVector<double> ListDimVar;
  SingleArrayInterpolation eRecInterp;
  //
  MyMatrix<double> matCoeffM;
  MyMatrix<double> matCoeffP;
  MyMatrix<int> matIdxM;
  MyMatrix<int> matIdxP;
  MyMatrix<int> MSK;
  int Ntotal;
  MyVector<double> ListVertPos;
  std::vector<int> FirstWetIndex;
  //
  double normU, normV;
};




#endif
