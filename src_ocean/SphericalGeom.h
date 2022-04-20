#ifndef SRC_OCEAN_SPHERICALGEOM_H_
#define SRC_OCEAN_SPHERICALGEOM_H_

#include "Basic_Ocean_types.h"
#include "Basic_string.h"
#include "Temp_common.h"
#include <string>
#include <vector>

bool IsPointInQuadArray(QuadArray const &eQuad, double const &eLon,
                        double const &eLat) {
  if (eLon < eQuad.MinLon)
    return false;
  if (eLon > eQuad.MaxLon)
    return false;
  if (eLat < eQuad.MinLat)
    return false;
  if (eLat > eQuad.MaxLat)
    return false;
  return true;
}

// Other attempt:
// http://stackoverflow.com/questions/11716268/point-in-polygon-algorithm

bool IsPointInside_Point(PairLL const &ePt, std::vector<PairLL> const &ListPt) {
  int nvert = ListPt.size();
  int i, j;
  bool c = false;
  for (i = 0, j = nvert - 1; i < nvert; j = i++) {
    if (((ListPt[i].eLat > ePt.eLat) != (ListPt[j].eLat > ePt.eLat)) &&
        (ePt.eLon < (ListPt[j].eLon - ListPt[i].eLon) *
                            (ePt.eLat - ListPt[i].eLat) /
                            (ListPt[j].eLat - ListPt[i].eLat) +
                        ListPt[i].eLon))
      c = !c;
  }
  return c;
}

bool IsPointInside(double const &testx, double const &testy,
                   std::vector<double> const &vertx,
                   std::vector<double> const &verty) {
  int nvert = vertx.size();
  int i, j;
  bool c = false;
  for (i = 0, j = nvert - 1; i < nvert; j = i++) {
    if (((verty[i] > testy) != (verty[j] > testy)) &&
        (testx <
         (vertx[j] - vertx[i]) * (testy - verty[i]) / (verty[j] - verty[i]) +
             vertx[i]))
      c = !c;
  }
  return c;
}

std::string ConvertLLDoubleToDegMinSec(double const &eLL) {
  double eDeg = floor(eLL);
  double res1 = eLL - eDeg;
  double eMin = floor(res1 * double(60));
  double res2 = 60 * res1 - eMin;
  double eSec = floor(res2 * double(60));
  std::string retStr =
      IntToString(eDeg) + "." + IntToString(eMin) + "." + IntToString(eSec);
  return retStr;
}

std::vector<PairLL> GetGridBoundary(MyMatrix<double> const &LON,
                                    MyMatrix<double> const &LAT,
                                    int const &iStart, int const &iEnd,
                                    int const &jStart, int const &jEnd) {
  int len = 2 * (iEnd - iStart) + 2 * (jEnd - jStart);
  std::vector<PairLL> ListPt(len);
  int idx = 0;
  for (int iEta = iStart; iEta <= iEnd - 1; iEta++) {
    int iXi = jStart;
    PairLL ePt{LON(iEta, iXi), LAT(iEta, iXi)};
    ListPt[idx] = ePt;
    idx++;
  }
  for (int iXi = jStart; iXi <= jEnd - 1; iXi++) {
    int iEta = iEnd;
    PairLL ePt{LON(iEta, iXi), LAT(iEta, iXi)};
    ListPt[idx] = ePt;
    idx++;
  }
  for (int iEta = iEnd; iEta >= iStart + 1; iEta--) {
    int iXi = jEnd;
    PairLL ePt{LON(iEta, iXi), LAT(iEta, iXi)};
    ListPt[idx] = ePt;
    idx++;
  }
  for (int iXi = jEnd; iXi >= jStart + 1; iXi--) {
    int iEta = iStart;
    PairLL ePt{LON(iEta, iXi), LAT(iEta, iXi)};
    ListPt[idx] = ePt;
    idx++;
  }
  return ListPt;
}

std::vector<PairLL> GetGridBoundaryTot(MyMatrix<double> const &LON,
                                       MyMatrix<double> const &LAT) {
  int nbRow = LON.rows();
  int nbCol = LON.cols();
  return GetGridBoundary(LON, LAT, 0, nbRow - 1, 0, nbCol - 1);
}

struct PairCoord {
  int i;
  int j;
};

PairCoord FindContaining(PairLL const &ePt, MyMatrix<double> const &LON,
                         MyMatrix<double> const &LAT) {
  int eta_rho = LON.rows();
  int xi_rho = LON.cols();
  int iStart, iEnd, jStart, jEnd;
  iStart = 0;
  iEnd = eta_rho - 1;
  jStart = 0;
  jEnd = xi_rho - 1;
  while (true) {
    std::vector<PairLL> ListPt1 =
        GetGridBoundary(LON, LAT, iStart, iEnd, jStart, jEnd);
    bool test1 = IsPointInside_Point(ePt, ListPt1);
    if (!test1)
      return {-1, -1};
    int iDiff = iEnd - iStart;
    int jDiff = jEnd - jStart;
    if (iDiff == 1 && jDiff == 1)
      break;
    if (iDiff > 1) {
      int iMid = int(roundl((float(iStart) + float(iEnd)) / double(2)));
      std::vector<PairLL> ListPt2 =
          GetGridBoundary(LON, LAT, iStart, iMid, jStart, jEnd);
      bool test2 = IsPointInside_Point(ePt, ListPt2);
      if (test2) {
        iEnd = iMid;
      } else {
        iStart = iMid;
      }
    }
    if (jDiff > 1) {
      int jMid = int(roundl((float(jStart) + float(jEnd)) / double(2)));
      std::vector<PairLL> ListPt3 =
          GetGridBoundary(LON, LAT, iStart, iEnd, jStart, jMid);
      bool test3 = IsPointInside_Point(ePt, ListPt3);
      if (test3) {
        jEnd = jMid;
      } else {
        jStart = jMid;
      }
    }
  }
  return {iStart, jStart};
}

std::vector<double> GetLLcoordinateXYZ(std::vector<double> const &eXYZ) {
  if (eXYZ.size() != 3) {
    std::cerr << "eXYZ should have size 3\n";
    throw TerminalException{1};
  }
  double eX = eXYZ[0];
  double eY = eXYZ[1];
  double eZ = eXYZ[2];
  double lat = asin(eZ);
  double eX2 = eX / cos(lat);
  double eY2 = eY / cos(lat);
  double lon = atan2(eY2, eX2);
  double pi = 3.141592653589792;
  double coef = double(180) / pi;
  double lonDeg = lon * coef;
  double latDeg = lat * coef;
  return {lonDeg, latDeg};
}

std::vector<double> GetXYZcoordinateLL(double const &LonDeg,
                                       double const &LatDeg) {
  double pi = 3.141592653589792;
  double lon = pi * LonDeg / double(180);
  double lat = pi * LatDeg / double(180);
  double x = cos(lon) * cos(lat);
  double y = sin(lon) * cos(lat);
  double z = sin(lat);
  return {x, y, z};
}

TripleXYZ ComputeTripleXYZ(MyMatrix<double> const &LONdeg,
                           MyMatrix<double> const &LATdeg) {
  double pi = 3.141592653589792;
  int eta_rho = LONdeg.rows();
  int xi_rho = LONdeg.cols();
  MyMatrix<double> X(eta_rho, xi_rho);
  MyMatrix<double> Y(eta_rho, xi_rho);
  MyMatrix<double> Z(eta_rho, xi_rho);
  for (int iEta = 0; iEta < eta_rho; iEta++)
    for (int iXi = 0; iXi < xi_rho; iXi++) {
      double lon = pi * LONdeg(iEta, iXi) / double(180);
      double lat = pi * LATdeg(iEta, iXi) / double(180);
      double x = cos(lon) * cos(lat);
      double y = sin(lon) * cos(lat);
      double z = sin(lat);
      X(iEta, iXi) = x;
      Y(iEta, iXi) = y;
      Z(iEta, iXi) = z;
    }
  return {X, Y, Z};
}

std::vector<double> GetXYZaverage(std::vector<double> const &sXYZ,
                                  std::vector<double> const &eXYZ,
                                  double const &sCoeff, double const &eCoeff) {
  int len = sXYZ.size();
  std::vector<double> retXYZ(len);
  double sum = 0;
  for (int i = 0; i < len; i++) {
    double eCoord = sXYZ[i] * sCoeff + eXYZ[i] * eCoeff;
    sum += eCoord * eCoord;
    retXYZ[i] = eCoord;
  }
  double norm = sqrt(sum);
  for (int i = 0; i < len; i++)
    retXYZ[i] /= norm;
  return retXYZ;
}

double GeodesicDistance(double const &LonDeg1, double const &LatDeg1,
                        double const &LonDeg2, double const &LatDeg2) {
  double pi = 3.141592653589792;
  double lon1 = pi * LonDeg1 / double(180);
  double lat1 = pi * LatDeg1 / double(180);
  double x1 = cos(lon1) * cos(lat1);
  double y1 = sin(lon1) * cos(lat1);
  double z1 = sin(lat1);

  double lon2 = pi * LonDeg2 / double(180);
  double lat2 = pi * LatDeg2 / double(180);
  double x2 = cos(lon2) * cos(lat2);
  double y2 = sin(lon2) * cos(lat2);
  double z2 = sin(lat2);
  double scalprod = x1 * x2 + y1 * y2 + z1 * z2;
  if (scalprod > 1)
    return 0;
  else
    return acos(scalprod);
}

double GeodesicDistance_pair(PairLL const &ePair1, PairLL const &ePair2) {
  return GeodesicDistance(ePair1.eLon, ePair1.eLat, ePair2.eLon, ePair2.eLat);
}

double GeodesicDistanceKM(double const &LonDeg1, double const &LatDeg1,
                          double const &LonDeg2, double const &LatDeg2) {
  double EarthRadius = 6370;
  return EarthRadius * GeodesicDistance(LonDeg1, LatDeg1, LonDeg2, LatDeg2);
}

double GeodesicDistanceKM_pair(PairLL const &ePair1, PairLL const &ePair2) {
  double EarthRadius = 6370;
  return EarthRadius * GeodesicDistance_pair(ePair1, ePair2);
}

double GeodesicDistanceM_General(double const &LonDeg1, double const &LatDeg1,
                                 double const &LonDeg2, double const &LatDeg2,
                                 bool const &IsSpherical) {
  if (IsSpherical)
    return GeodesicDistanceKM(LonDeg1, LatDeg1, LonDeg2, LatDeg2) * 1000;
  double DeltaLon = LonDeg1 - LonDeg2;
  double DeltaLat = LatDeg1 - LatDeg2;
  return sqrt(DeltaLon * DeltaLon + DeltaLat * DeltaLat);
}

double SphericalCoordinateArea(double const &eLon1, double const &eLon2,
                               double const &eLon3, double const &eLat1,
                               double const &eLat2, double const &eLat3) {
  double dist1 = GeodesicDistance(eLon1, eLat1, eLon2, eLat2);
  double dist2 = GeodesicDistance(eLon2, eLat2, eLon3, eLat3);
  double dist3 = GeodesicDistance(eLon3, eLat3, eLon1, eLat1);
  double DistS = (dist1 + dist2 + dist3) / double(2);
  double eTan1 = tan(DistS / double(2));
  double eTan2 = tan((DistS - dist1) / double(2));
  double eTan3 = tan((DistS - dist2) / double(2));
  double eTan4 = tan((DistS - dist3) / double(2));
  double eProd = eTan1 * eTan2 * eTan3 * eTan4;
  double sqrtProd = sqrt(eProd);
  double area = double(4) * atan(sqrtProd);
  return area;
}

double SphericalCoordinateAreaKM(double const &eLon1, double const &eLon2,
                                 double const &eLon3, double const &eLat1,
                                 double const &eLat2, double const &eLat3) {
  double EarthRadius = 6370;
  double area =
      SphericalCoordinateArea(eLon1, eLon2, eLon3, eLat1, eLat2, eLat3);
  return EarthRadius * EarthRadius * area;
}

double SphericalCoordinateAreaKM_pair(PairLL const& pair1, PairLL const& pair1, PairLL const& pair1) {
  return SphericalCoordinateAreaKM_pair(pair1.eLon, pair2.eLon, pair3.eLon, pair1.eLat, pair2.eLat, pair3.eLat);
}

template <typename T>
bool get_line_intersection(T const &p0_x, T const &p0_y, T const &p1_x,
                           T const &p1_y, T const &p2_x, T const &p2_y,
                           T const &p3_x, T const &p3_y, T *i_x, T *i_y) {
  T s1_x = p1_x - p0_x;
  T s1_y = p1_y - p0_y;
  T s2_x = p3_x - p2_x;
  T s2_y = p3_y - p2_y;

  T s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) /
        (-s2_x * s1_y + s1_x * s2_y);
  T t = (s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) /
        (-s2_x * s1_y + s1_x * s2_y);

  if (s >= 0 && s <= 1 && t >= 0 && t <= 1) {
    *i_x = p0_x + (t * s1_x);
    *i_y = p0_y + (t * s1_y);
    return true;
  }
  return false;
}

// The LonDomain1 / LatDomain1 is the big domain.
// The question is whether LonDomain2 / LatDomain2 intersects with Domain1 or is
// contained in it.
bool IsContainedInDomain(std::vector<double> const &LonDomain1,
                         std::vector<double> const &LatDomain1,
                         std::vector<double> const &LonDomain2,
                         std::vector<double> const &LatDomain2) {
  double deltaLat = VectorMax(LatDomain2) - VectorMin(LatDomain2);
  double deltaLon = VectorMax(LonDomain2) - VectorMin(LonDomain2);
  std::cerr << "deltaLat=" << deltaLat << " deltaLon=" << deltaLon << "\n";
  if (deltaLon > 300)
    return false;
  //
  // Fast heuristic method. Works in many cases
  //
  auto fastConclusionExclusion = [](std::vector<double> const &Field1,
                                    std::vector<double> const &Field2) -> bool {
    double min1 = VectorMin(Field1);
    double max1 = VectorMax(Field1);
    double min2 = VectorMin(Field2);
    double max2 = VectorMax(Field2);
    if (max2 < min1)
      return true;
    if (max1 < min2)
      return true;
    return false;
  };
  if (fastConclusionExclusion(LonDomain1, LonDomain2))
    return false;
  if (fastConclusionExclusion(LatDomain1, LatDomain2))
    return false;
  //
  // Now testing inclusions
  //
  int siz1 = LonDomain1.size();
  std::vector<PairLL> ListPt1(siz1);
  for (int i1 = 0; i1 < siz1; i1++) {
    PairLL ePt1{LonDomain1[i1], LatDomain1[i1]};
    ListPt1[i1] = ePt1;
  }
  int siz2 = LonDomain2.size();
  for (int i2 = 0; i2 < siz2; i2++) {
    PairLL ePt2{LonDomain2[i2], LatDomain2[i2]};
    bool test = IsPointInside_Point(ePt2, ListPt1);
    if (test)
      return true;
  }
  for (int i1 = 0; i1 < siz1; i1++) {
    int i1N = NextIdx(siz1, i1);
    double p0_x = LonDomain1[i1];
    double p0_y = LatDomain1[i1];
    double p1_x = LonDomain1[i1N];
    double p1_y = LatDomain1[i1N];
    for (int i2 = 0; i2 < siz2; i2++) {
      int i2N = NextIdx(siz2, i2);
      double p2_x = LonDomain2[i2];
      double p2_y = LatDomain2[i2];
      double p3_x = LonDomain2[i2N];
      double p3_y = LatDomain2[i2N];
      double int_x, int_y;
      bool test = get_line_intersection(p0_x, p0_y, p1_x, p1_y, p2_x, p2_y,
                                        p3_x, p3_y, &int_x, &int_y);
      if (test)
        return true;
    }
  }
  return false;
}

struct DiscInfo {
  PairLL eSample;
  PairLL avgPoint;
  double SpreadLon;
  double SpreadLat;
  double MaxSpread;
};

struct KTreeElt {
  std::vector<PairLL> ListPt;
  DiscInfo eDisc;
  //
  bool IsSplit;
  int iSub1;
  int iSub2;
};

MyMatrix<double>
GRID_GetMaxRadiusInfluence_kernel(MyMatrix<double> const &LON_rho,
                                  MyMatrix<double> const &LAT_rho,
                                  std::string const &strDistType) {
  int eta_rho = LON_rho.rows();
  int xi_rho = LON_rho.cols();
  int eChoice = -1;
  if (strDistType == "spherical")
    eChoice = 1;
  if (strDistType == "euclidean")
    eChoice = 2;
  if (eChoice == -1) {
    std::cerr << "Possible choices are spherical and euclidean\n";
    throw TerminalException{1};
  }
  MyMatrix<int> ListDir(4, 2);
  ListDir(0, 0) = 1;
  ListDir(0, 1) = 1;
  ListDir(1, 0) = 1;
  ListDir(1, 1) = -1;
  ListDir(2, 0) = -1;
  ListDir(2, 1) = 1;
  ListDir(3, 0) = -1;
  ListDir(3, 1) = -1;
  MyMatrix<double> RadiusInflMatrix(eta_rho, xi_rho);
  for (int iEta = 0; iEta < eta_rho; iEta++)
    for (int iXi = 0; iXi < xi_rho; iXi++) {
      double lon1 = LON_rho(iEta, iXi);
      double lat1 = LAT_rho(iEta, iXi);
      double MaxDist = 0;
      for (int iN = 0; iN < 4; iN++) {
        int iEtaN = iEta + ListDir(iN, 0);
        int iXiN = iXi + ListDir(iN, 1);
        if (iEtaN >= 0 && iEtaN < eta_rho && iXiN >= 0 && iXiN < xi_rho) {
          double lon2 = LON_rho(iEtaN, iXiN);
          double lat2 = LAT_rho(iEtaN, iXiN);
          double eDist;
          if (eChoice == 1) {
            eDist = GeodesicDistanceKM(lon1, lat1, lon2, lat2);
          } else {
            double dx = lon1 - lon2;
            double dy = lat1 - lat2;
            eDist = sqrt(dx * dx + dy * dy);
          }
          if (eDist > MaxDist)
            MaxDist = eDist;
        }
      }
      RadiusInflMatrix(iEta, iXi) = MaxDist;
    }
  return RadiusInflMatrix;
}

DiscInfo KTree_ComputeDisc(std::vector<PairLL> const &ListPt) {
  double MinLon = ListPt[0].eLon;
  double MaxLon = ListPt[0].eLon;
  double MinLat = ListPt[0].eLat;
  double MaxLat = ListPt[0].eLat;
  int siz = ListPt.size();
  double SumLon = 0;
  double SumLat = 0;
  for (int i = 0; i < siz; i++) {
    double eLon = ListPt[i].eLon;
    double eLat = ListPt[i].eLat;
    SumLon += eLon;
    SumLat += eLat;
    if (eLon < MinLon)
      MinLon = eLon;
    if (eLon > MaxLon)
      MaxLon = eLon;
    if (eLat < MinLat)
      MinLat = eLat;
    if (eLat > MaxLat)
      MaxLat = eLat;
  }
  double avgLon = SumLon / double(siz);
  double avgLat = SumLat / double(siz);
  double dist12 = GeodesicDistanceKM(MinLon, MinLat, MinLon, MaxLat);
  double dist23 = GeodesicDistanceKM(MinLon, MaxLat, MaxLon, MaxLat);
  double dist34 = GeodesicDistanceKM(MaxLon, MaxLat, MaxLon, MinLat);
  double dist41 = GeodesicDistanceKM(MaxLon, MinLat, MinLon, MinLat);
  //
  double SpreadLat = (dist12 + dist34) / double(2);
  double SpreadLon = (dist23 + dist41) / double(2);
  double MinDistKM = 20000;
  int idxMin = -1;
  for (int i = 0; i < siz; i++) {
    double eLon = ListPt[i].eLon;
    double eLat = ListPt[i].eLat;
    double dist = GeodesicDistanceKM(eLon, eLat, avgLon, avgLat);
    if (dist < MinDistKM) {
      MinDistKM = dist;
      idxMin = i;
    }
  }
  PairLL eSample = ListPt[idxMin];
  double eSampleLon = eSample.eLon;
  double eSampleLat = eSample.eLat;
  double MaxSpread = 0;
  for (int i = 0; i < siz; i++) {
    double eLon = ListPt[i].eLon;
    double eLat = ListPt[i].eLat;
    double dist = GeodesicDistanceKM(eLon, eLat, eSampleLon, eSampleLat);
    if (dist > MaxSpread)
      MaxSpread = dist;
  }
  PairLL avgPoint{avgLon, avgLat};
  return {eSample, avgPoint, SpreadLon, SpreadLat, MaxSpread};
}

std::vector<KTreeElt> SplitKTreeElt(KTreeElt const &eKD) {
  std::vector<PairLL> NewListPt1;
  std::vector<PairLL> NewListPt2;
  if (eKD.eDisc.SpreadLon > eKD.eDisc.SpreadLat) {
    for (auto &ePt : eKD.ListPt) {
      if (ePt.eLon < eKD.eDisc.avgPoint.eLon)
        NewListPt1.push_back(ePt);
      else
        NewListPt2.push_back(ePt);
    }
  } else {
    for (auto &ePt : eKD.ListPt) {
      if (ePt.eLat < eKD.eDisc.avgPoint.eLat)
        NewListPt1.push_back(ePt);
      else
        NewListPt2.push_back(ePt);
    }
  }
  //
  KTreeElt eKD1{{}, eKD.eDisc, true, -2, -2};
  KTreeElt eKD2{NewListPt1, KTree_ComputeDisc(NewListPt1), false, -1, -1};
  KTreeElt eKD3{NewListPt2, KTree_ComputeDisc(NewListPt2), false, -1, -1};
  return {eKD1, eKD2, eKD3};
}

std::vector<KTreeElt>
KTree_GetDecomposition(std::vector<PairLL> const &ListPtCoast) {
  std::vector<KTreeElt> TList;
  std::vector<int> IsTreated;
  int MaxNumberPerCell = 100;
  auto SplitComponent = [&](int const &iComp) -> void {
    int len = TList.size();
    std::vector<KTreeElt> NList = SplitKTreeElt(TList[iComp]);
    TList[iComp] = NList[0];
    TList[iComp].iSub1 = len;
    TList[iComp].iSub2 = len + 1;
    TList.push_back(NList[1]);
    TList.push_back(NList[2]);
    IsTreated.push_back(0);
    IsTreated.push_back(0);
  };
  KTreeElt eElt{ListPtCoast, KTree_ComputeDisc(ListPtCoast), false, -1, -1};
  TList.push_back(eElt);
  IsTreated.push_back(0);
  while (true) {
    int siz = TList.size();
    //    std::cerr << "siz=" << siz << "\n";
    bool IsFinished = true;
    for (int i = 0; i < siz; i++)
      if (IsTreated[i] == 0) {
        IsTreated[i] = 1;
        IsFinished = false;
        int len = TList[i].ListPt.size();
        if (len > MaxNumberPerCell)
          SplitComponent(i);
      }
    if (IsFinished)
      break;
  }
  //  for (int i=0; i<siz; i++) {
  //    std::cerr << "i=" << i << " split=" << TList[i].IsSplit << " |ListPt|="
  //    << TList[i].ListPt.size() << " MaxSpread=" << TList[i].eDisc.MaxSpread
  //    << "\n";
  //  }
  return TList;
}

double ShortestDistance(std::vector<KTreeElt> const &ListKT, PairLL const &ePt,
                        double &UpperEstimate) {
  double TheDist = UpperEstimate;
  double eLon = ePt.eLon;
  double eLat = ePt.eLat;
  //  std::cerr << "|ListKT|=" << ListKT.size() << "\n";
  std::vector<int> ListIdx{0};
  while (true) {
    if (ListIdx.size() == 0)
      break;
    //    std::cerr << "Before creating NewListIdx\n";
    std::vector<int> NewListIdx;
    for (int &eVal : ListIdx) {
      if (!ListKT[eVal].IsSplit) {
        int len = ListKT[eVal].ListPt.size();
        for (int i = 0; i < len; i++) {
          double nLon = ListKT[eVal].ListPt[i].eLon;
          double nLat = ListKT[eVal].ListPt[i].eLat;
          double nDist = GeodesicDistanceKM(nLon, nLat, eLon, eLat);
          if (nDist < TheDist)
            TheDist = nDist;
        }
      } else {
        double SampleLon = ListKT[eVal].eDisc.eSample.eLon;
        double SampleLat = ListKT[eVal].eDisc.eSample.eLat;
        double nDist = GeodesicDistanceKM(SampleLon, SampleLat, eLon, eLat);
        double LowerBound = nDist - ListKT[eVal].eDisc.MaxSpread;
        if (LowerBound < TheDist) {
          NewListIdx.push_back(ListKT[eVal].iSub1);
          NewListIdx.push_back(ListKT[eVal].iSub2);
        }
      }
    }
    ListIdx = NewListIdx;
  }
  return TheDist;
}

std::vector<double>
GetUpperEstimateMinDist(std::vector<PairLL> const &ListPt1,
                        std::vector<PairLL> const &ListPt2) {
  int nbPt2 = ListPt2.size();
  int nbPt1 = ListPt1.size();
  if (nbPt1 == 0) {
    std::cerr << "The list ListPt1 should not be empty\n";
    std::cerr << "nbPt1=" << nbPt1 << "\n";
    throw TerminalException{1};
  }
  //  std::cerr << "nbPt1=" << nbPt1 << " nbPt2=" << nbPt2 << "\n";
  auto fDist = [&](int const &i1, int const &i2) -> double {
    double eLon1 = ListPt1[i1].eLon;
    double eLat1 = ListPt1[i1].eLat;
    double eLon2 = ListPt2[i2].eLon;
    double eLat2 = ListPt2[i2].eLat;
    return GeodesicDistanceKM(eLon1, eLat1, eLon2, eLat2);
  };
  int i1 = 0;
  std::vector<double> ListUpperEst(nbPt2);
  for (int i2 = 0; i2 < nbPt2; i2++) {
    double eDist = fDist(i1, i2);
    //    std::cerr << "i2=" << i2 << "\n";
    for (int iter = 0; iter < 4; iter++) {
      int i1New = rand() % nbPt1;
      double nDist = fDist(i1New, i2);
      if (nDist < eDist) {
        eDist = nDist;
        i1 = i1New;
        // std::cerr << "  i1New=" << i1 << " eDist=" << eDist << "\n";
      }
    }
    while (true) {
      bool DoSomething = false;
      //      std::cerr << "  i1=" << i1 << " eDist=" << eDist << "\n";
      if (i1 > 0) {
        double nDist = fDist(i1 - 1, i2);
        if (nDist < eDist) {
          eDist = nDist;
          i1 = i1 - 1;
          DoSomething = true;
        }
      }
      if (i1 < nbPt1 - 1) {
        double nDist = fDist(i1 + 1, i2);
        if (nDist < eDist) {
          eDist = nDist;
          i1 = i1 + 1;
          DoSomething = true;
        }
      }
      if (!DoSomething)
        break;
    }
    ListUpperEst[i2] = eDist;
  }
  //  std::cerr << "Now leaving returning ListUpperEst\n";
  return ListUpperEst;
}

std::vector<double>
GetListMinimalDistances(std::vector<PairLL> const &ListPtCoast,
                        std::vector<PairLL> const &ListPt) {
  int nbPt = ListPt.size();
  std::vector<double> ListUpperEst =
      GetUpperEstimateMinDist(ListPtCoast, ListPt);
  std::vector<KTreeElt> ListKT = KTree_GetDecomposition(ListPtCoast);
  std::vector<double> ListShortest(nbPt);
  double TotalDefect = 0;
  for (int iPt = 0; iPt < nbPt; iPt++) {
    double eEst = ListUpperEst[iPt];
    PairLL ePt = ListPt[iPt];
    double eMinDist = ShortestDistance(ListKT, ePt, eEst);
    ListShortest[iPt] = eMinDist;
    double eDefect = eEst - eMinDist;
    TotalDefect += eDefect;
    //    std::cerr << "iPt=" << iPt << " eMinDist=" << eMinDist << " eEst=" <<
    //    eEst << "\n";
  }
  //  std::cerr << "TotalDefect=" << TotalDefect << "\n";
  return ListShortest;
}

void TwoPiNormalization(double &TheAng) {
  double ThePi = 3.141592653589792;
  if (TheAng < -ThePi) {
    TheAng += double(2) * ThePi;
  }
  if (TheAng > ThePi) {
    TheAng -= double(2) * ThePi;
  }
}

int MySign(double &TheVal) {
  if (TheVal > 0)
    return 1;
  if (TheVal < 0)
    return -1;
  return 0;
}

MyMatrix<double> CreateAngleMatrix(MyMatrix<double> const &LON_rho,
                                   MyMatrix<double> const &LAT_rho) {
  int eta_rho = LON_rho.rows();
  int xi_rho = LON_rho.cols();
  int eta_v = eta_rho - 1;
  int xi_v = xi_rho;
  MyMatrix<double> LONrad_v(eta_v, xi_v);
  MyMatrix<double> LATrad_v(eta_v, xi_v);
  MyMatrix<double> azim(eta_v - 1, xi_v);
  double ThePi = 3.141592653589792;
  double DegTwoRad = ThePi / double(180);
  for (int iEta = 0; iEta < eta_v; iEta++)
    for (int iXi = 0; iXi < xi_v; iXi++) {
      double eLon = (LON_rho(iEta, iXi) + LON_rho(iEta + 1, iXi)) / double(2);
      double eLat = (LAT_rho(iEta, iXi) + LAT_rho(iEta + 1, iXi)) / double(2);
      LONrad_v(iEta, iXi) = eLon * DegTwoRad;
      LATrad_v(iEta, iXi) = eLat * DegTwoRad;
    }
  for (int iEta = 0; iEta < eta_v - 1; iEta++)
    for (int iXi = 0; iXi < xi_v; iXi++) {
      double phi1 = LATrad_v(iEta, iXi);
      double xlam1 = LONrad_v(iEta, iXi);
      double phi2 = LATrad_v(iEta + 1, iXi);
      double xlam2 = LONrad_v(iEta + 1, iXi);
      double TPSI2 = tan(phi2);
      double dlam = xlam2 - xlam1;
      TwoPiNormalization(dlam);
      double cta12 = (cos(phi1) * TPSI2 - sin(phi1) * cos(dlam)) / sin(dlam);
      double eAzim = atan(double(1) / cta12);
      int signAzim = MySign(eAzim);
      int signDlam = MySign(dlam);
      int eFact2;
      if (signDlam != signAzim)
        eFact2 = 1;
      else
        eFact2 = 0;
      int eFact1 = -signAzim;
      double fAzim = eAzim + ThePi * eFact1 * eFact2;
      azim(iEta, iXi) = fAzim;
    }
  MyMatrix<double> ANG_rho(eta_rho, xi_rho);
  for (int iEta = 1; iEta < eta_v; iEta++)
    for (int iXi = 0; iXi < xi_v; iXi++)
      ANG_rho(iEta, iXi) = ThePi / double(2) - azim(iEta - 1, iXi);
  for (int iXi = 0; iXi < xi_v; iXi++) {
    ANG_rho(0, iXi) = ANG_rho(1, iXi);
    ANG_rho(eta_rho - 1, iXi) = ANG_rho(eta_rho - 2, iXi);
  }
  return ANG_rho;
}

constexpr double GetPI() { return 3.1415926535; }

#endif  // SRC_OCEAN_SPHERICALGEOM_H_
