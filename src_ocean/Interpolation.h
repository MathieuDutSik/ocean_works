// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_OCEAN_INTERPOLATION_H_
#define SRC_OCEAN_INTERPOLATION_H_

#include "SphericalGeom.h"
#include "Triangulations.h"
#include <algorithm>
#include <utility>
#include <vector>

#ifdef USE_OPENCV_LIBRARY
#include <opencv2/flann.hpp>
#include <opencv2/opencv.hpp>
#endif

struct SinglePartInterp {
  int eEta, eXi;
  double eCoeff;
};

struct SingleRecInterp {
  bool status;
  std::vector<SinglePartInterp> LPart;
};

struct QuadCoordinate {
  double MinLon;
  double MaxLon;
  double MinLat;
  double MaxLat;
};

std::vector<double> DetermineCoefficient(std::vector<double> const &X,
                                         std::vector<double> const &Y,
                                         double const &eX, double const &eY) {
  MyMatrix<double> A(3, 3);
  for (int i = 0; i < 3; i++) {
    A(0, i) = 1;
    A(1, i) = X[i];
    A(2, i) = Y[i];
  }
  MyVector<double> V(3);
  V(0) = 1;
  V(1) = eX;
  V(2) = eY;
  MyMatrix<double> eInv = A.inverse().eval();
  MyVector<double> eProduct = eInv * V;
  std::vector<double> LCoeff(3);
  for (int i = 0; i < 3; i++)
    LCoeff[i] = eProduct(i);
  /*
  double deltaX=eX;
  double deltaY=eY;
  for (int i=0; i<3; i++) {
    deltaX=deltaX - LCoeff[i]*X[i];
    deltaY=deltaY - LCoeff[i]*Y[i];
  }
//  std::cerr << "deltaX=" << deltaX << " deltaY=" << deltaY << "\n"; */
  return LCoeff;
}

bool TestFeasibilityByQuad(QuadCoordinate const &eQuad, double const &eLon,
                           double const &eLat) {
  //  std::cerr << "eLon=" << eLon << " eLat=" << eLat << "\n";
  //  std::cerr << "LON(min/max)=" << eQuad.MinLon << " / " << eQuad.MaxLon <<
  //  "\n"; std::cerr << "LAT(min/max)=" << eQuad.MinLat << " / " <<
  //  eQuad.MaxLat << "\n";
  if (eLon > eQuad.MinLon && eLon < eQuad.MaxLon && eLat > eQuad.MinLat &&
      eLat < eQuad.MaxLat)
    return true;
  return false;
}

std::vector<SingleRecInterp>
TRIG_FIND_ELE_Kernel(MyMatrix<int> const &INE, MyMatrix<double> const &X,
                     MyMatrix<double> const &Y, QuadCoordinate const &eQuad,
                     MyMatrix<double> const &ListXY,
                     bool const &AllowExtrapolation) {
  double THR = 1e-10;
  int mnp = X.rows();
  //  for (int i=0; i<mnp; i++)
  //    std::cerr << "i=" << i << " x=" << X(i) << " y=" << Y(i) << "\n";
  auto IsCorrect = [&](int const &ie, double const &Xp,
                       double const &Yp) -> bool {
    int ki = INE(ie, 0);
    int kj = INE(ie, 1);
    int kk = INE(ie, 2);
    double xi = X(ki);
    double yi = Y(ki);
    double xj = X(kj);
    double yj = Y(kj);
    double xk = X(kk);
    double yk = Y(kk);
    /*
    double area=xi*(yj-yk) + xj*(yk-yi) + xk*(yi-yj);
    std::cerr << "-------------------------------------\n";
    std::cerr << "ie=" << ie << " area=" << area << "\n";
    std::cerr << "kijk=" << ki << "," << kj << "," << kk << "\n";
    std::cerr << "x(ijk)=" << xi << "," << xj << "," << xk << "\n";
    std::cerr << "y(ijk)=" << yi << "," << yj << "," << yk << "\n";
    std::cerr << "-------------------------------------\n";*/
    double f1, f2, f3;
    f1 = xi * (yj - Yp) + xj * (Yp - yi) + Xp * (yi - yj);
    f2 = xj * (yk - Yp) + xk * (Yp - yj) + Xp * (yj - yk);
    f3 = xk * (yi - Yp) + xi * (Yp - yk) + Xp * (yk - yi);
    //    double sumF=f1 + f2 + f3;
    //    std::cerr << "sumF=" << sumF << "\n";
    if (f1 > -THR && f2 > -THR && f3 > -THR)
      return true;
    return false;
  };
  double dx = 0;
  double dy = 0;
  int nbEle = INE.rows();
  for (int ie = 0; ie < nbEle; ie++) {
    int ki = INE(ie, 0);
    int kj = INE(ie, 1);
    int kk = INE(ie, 2);
    double xi = X(ki);
    double yi = Y(ki);
    double xj = X(kj);
    double yj = Y(kj);
    double xk = X(kk);
    double yk = Y(kk);
    dx = std::max(dx, std::abs(xi - xj));
    dx = std::max(dx, std::abs(xi - xk));
    dx = std::max(dx, std::abs(xj - xk));
    dy = std::max(dy, std::abs(yi - yj));
    dy = std::max(dy, std::abs(yi - yk));
    dy = std::max(dy, std::abs(yj - yk));
  }
  std::vector<int> ListAdjTrig =
      GetUnstructuredTriangleAdjInfo_vectint(INE, mnp);
  auto SearchElement_V1 = [&](double const &eX, double const &eY) -> int {
    for (int iele = 0; iele < nbEle; iele++)
      if (IsCorrect(iele, eX, eY))
        return iele;
    return -1;
  };
  auto DistCentTriangle = [&](double const &eX, double const &eY,
                              int const &iEle) -> double {
    int i1 = INE(iEle, 0);
    int i2 = INE(iEle, 1);
    int i3 = INE(iEle, 2);
    double eXcent = (X(i1) + X(i2) + X(i3)) / static_cast<double>(3);
    double eYcent = (X(i1) + X(i2) + X(i3)) / static_cast<double>(3);
    double delX = eXcent - eX;
    double delY = eYcent - eY;
    return sqrt(delX * delX + delY * delY);
  };
  auto SearchElement = [&](double const &eX, double const &eY,
                           int const &iEltStart) -> int {
    if (IsCorrect(iEltStart, eX, eY))
      return iEltStart;
    int iEleWork = iEltStart;
    double distCurr = DistCentTriangle(eX, eY, iEleWork);
    //    std::cerr << "distCurr=" << distCurr << "\n";
    int nbIter = 0;
    while (true) {
      bool DoSomething = false;
      nbIter++;
      for (int i = 0; i < 3; i++) {
        int iEleAdj = ListAdjTrig[3 * iEleWork + i];
        if (iEleAdj != -1 && !DoSomething) {
          double eDist = DistCentTriangle(eX, eY, iEleAdj);
          if (eDist < distCurr) {
            iEleWork = iEleAdj;
            distCurr = eDist;
            DoSomething = true;
            if (IsCorrect(iEleWork, eX, eY)) {
              return iEleWork;
            }
          }
        }
      }
      if (!DoSomething)
        break;
    }
    return SearchElement_V1(eX, eY);
  };
  int nbPoint = ListXY.cols();
  int ielePrev = 0;
  std::vector<SingleRecInterp> LRec(nbPoint);
  int nbExtrapolation = 0;
  for (int iPoint = 0; iPoint < nbPoint; iPoint++) {
    double Xp = ListXY(0, iPoint);
    double Yp = ListXY(1, iPoint);
    int eElt = SearchElement(Xp, Yp, ielePrev);
    if (eElt >= 0)
      ielePrev = eElt;
    //    std::cerr << "iPoint=" << iPoint << " eElt=" << eElt << "\n";
    SingleRecInterp eRec;
    if (eElt == -1) {
      if (AllowExtrapolation) {
        nbExtrapolation++;
        bool IsFirst = true;
        double MinDist = 0;
        int idx = -1;
        for (int ip = 0; ip < mnp; ip++) {
          double dx = X(ip) - Xp;
          double dy = Y(ip) - Yp;
          double dist_sqr = dx * dx + dy * dy;
          if (IsFirst) {
            idx = ip;
            MinDist = dist_sqr;
            IsFirst = false;
          } else {
            if (dist_sqr < MinDist) {
              idx = ip;
              MinDist = dist_sqr;
            }
          }
        }
        SinglePartInterp ePart = {idx, 0, static_cast<double>(1)};
        eRec = {true, {ePart}};
      } else {
        eRec = {false, {}};
      }
    } else {
      std::vector<int> LEta(3);
      std::vector<double> Xcall(3), Ycall(3);
      for (int i = 0; i < 3; i++) {
        int IP = INE(eElt, i);
        LEta[i] = IP;
        Xcall[i] = X(IP);
        Ycall[i] = Y(IP);
      }
      std::vector<double> LCoeff = DetermineCoefficient(Xcall, Ycall, Xp, Yp);
      std::vector<SinglePartInterp> LPart(3);
      for (int i = 0; i < 3; i++) {
        LPart[i] = {LEta[i], 0, LCoeff[i]};
      }
      eRec = {true, LPart};
    };
    LRec[iPoint] = eRec;
  }
  std::cerr << "TRIG_FIND_ELE_Kernel nbExtrapolation=" << nbExtrapolation
            << "\n";
  return LRec;
}

std::vector<SingleRecInterp>
TRIG_FIND_ELE(MyMatrix<int> const &INE, MyMatrix<double> const &X,
              MyMatrix<double> const &Y, QuadCoordinate const &eQuad,
              MyMatrix<double> const &ListXY, bool const &AllowExtrapolation) {
  int nbElt = ListXY.cols();
  std::vector<PairLL> ListElt(nbElt);
  for (int iElt = 0; iElt < nbElt; iElt++) {
    double eLon = ListXY(0, iElt);
    double eLat = ListXY(1, iElt);
    PairLL ePair{eLon, eLat};
    ListElt[iElt] = ePair;
  }
  std::cerr << "We have ListElt\n";
  std::pair<std::vector<int>, std::vector<int>> ePairSort =
      SortingLists(ListElt);
  std::cerr << "We have ePairSort\n";
  std::vector<int> SortedToOrig = ePairSort.first;
  std::vector<int> OrigToSorted = ePairSort.second;
  MyMatrix<double> ListXYsort(2, nbElt);
  for (int iElt = 0; iElt < nbElt; iElt++) {
    int iEltOrig = SortedToOrig[iElt];
    for (int u = 0; u < 2; u++)
      ListXYsort(u, iElt) = ListXY(u, iEltOrig);
  }
  std::cerr << "We have ListXYsort\n";
  std::vector<SingleRecInterp> eVect =
      TRIG_FIND_ELE_Kernel(INE, X, Y, eQuad, ListXYsort, AllowExtrapolation);
  std::cerr << "We have eVect\n";
  std::vector<SingleRecInterp> retVect(nbElt);
  for (int iEltOrig = 0; iEltOrig < nbElt; iEltOrig++) {
    int iEltSort = OrigToSorted[iEltOrig];
    retVect[iEltOrig] = eVect[iEltSort];
  }
  std::cerr << "We have retVect\n";
  return retVect;
}

MyMatrix<int> GetDirection() {
  MyMatrix<int> eMat(2, 4);
  eMat(0, 0) = 0;
  eMat(1, 0) = 0;
  eMat(0, 1) = 1;
  eMat(1, 1) = 0;
  eMat(0, 2) = 0;
  eMat(1, 2) = 1;
  eMat(0, 3) = 1;
  eMat(1, 3) = 1;
  return eMat;
}

std::vector<SingleRecInterp> FD_FIND_ELE(CoordGridArrayFD const &CoordGridArr,
                                         QuadCoordinate const &eQuad,
                                         MyMatrix<double> const &ListXY,
                                         bool const &AllowExtrapolation) {
  double THR = 1e-10;
  MyMatrix<int> MatDir = GetDirection();
  int nbEta = CoordGridArr.LON.rows();
  int nbXi = CoordGridArr.LON.cols();
  struct DataStr {
    int i;
    int j;
    double lon;
    double lat;
  };
  std::vector<DataStr> ListWetEntry;
  for (int i = 0; i < nbEta; i++)
    for (int j = 0; j < nbXi; j++)
      if (CoordGridArr.MSK(i, j) == 1) {
        double lon = CoordGridArr.LON(i, j);
        double lat = CoordGridArr.LAT(i, j);
        ListWetEntry.push_back({i, j, lon, lat});
      }
  auto FindCoefficient = [&](int const &eEta, int const &eXi, double const &eX,
                             double const &eY) -> std::vector<double> {
    std::vector<double> X, Y, LCoeff;
    double TheMin;
    X = {CoordGridArr.LON(eEta, eXi), CoordGridArr.LON(eEta + 1, eXi),
         CoordGridArr.LON(eEta, eXi + 1)};
    Y = {CoordGridArr.LAT(eEta, eXi), CoordGridArr.LAT(eEta + 1, eXi),
         CoordGridArr.LAT(eEta, eXi + 1)};
    //    std::cerr << "Before DetermineCoefficient 1\n";
    LCoeff = DetermineCoefficient(X, Y, eX, eY);
    //    std::cerr << "After DetermineCoefficient 1\n";
    TheMin = VectorMin(LCoeff);
    if (TheMin > -THR) {
      //      std::cerr << "1: LCoeff=[" << LCoeff[0] << "," << LCoeff[1] << ","
      //      << LCoeff[2] << "]\n";
      double lambda1 = LCoeff[1];
      double lambda2 = LCoeff[2];
      double eCoeff00 = (1 - lambda1) * (1 - lambda2);
      double eCoeff10 = lambda1 * (1 - lambda2);
      double eCoeff01 = (1 - lambda1) * lambda2;
      double eCoeff11 = lambda1 * lambda2;
      //      std::cerr << "eX=" << eX << " eY=" << eY << "\n";
      //      std::cerr << "LON=" << CoordGridArr.LON(eEta,eXi) << "," <<
      //      CoordGridArr.LON(eEta+1,eXi) << "," <<
      //      CoordGridArr.LON(eEta,eXi+1) << "," <<
      //      CoordGridArr.LON(eEta+1,eXi+1) << "\n"; std::cerr << "LAT=" <<
      //      CoordGridArr.LAT(eEta,eXi) << "," << CoordGridArr.LAT(eEta+1,eXi)
      //      << "," << CoordGridArr.LAT(eEta,eXi+1) << "," <<
      //      CoordGridArr.LAT(eEta+1,eXi+1) << "\n"; std::cerr << "1:
      //      eCoeff00=" << eCoeff00 << " 10=" << eCoeff10 << " 01=" << eCoeff01
      //      << " 11=" << eCoeff11 << "\n";
      return {eCoeff00, eCoeff10, eCoeff01, eCoeff11};
      //
    }
    X = {CoordGridArr.LON(eEta + 1, eXi + 1), CoordGridArr.LON(eEta + 1, eXi),
         CoordGridArr.LON(eEta, eXi + 1)};
    Y = {CoordGridArr.LAT(eEta + 1, eXi + 1), CoordGridArr.LAT(eEta + 1, eXi),
         CoordGridArr.LAT(eEta, eXi + 1)};
    //    std::cerr << "Before DetermineCoefficient 2\n";
    LCoeff = DetermineCoefficient(X, Y, eX, eY);
    //    std::cerr << "After DetermineCoefficient 2\n";
    TheMin = VectorMin(LCoeff);
    if (TheMin > -THR) {
      //      std::cerr << "2: LCoeff=[" << LCoeff[0] << "," << LCoeff[1] << ","
      //      << LCoeff[2] << "]\n";
      double lambda1 = LCoeff[2];
      double lambda2 = LCoeff[1];
      double eCoeff00 = lambda1 * lambda2;
      double eCoeff10 = (1 - lambda1) * lambda2;
      double eCoeff01 = lambda1 * (1 - lambda2);
      double eCoeff11 = (1 - lambda1) * (1 - lambda2);
      //      std::cerr << "2: eCoeff00=" << eCoeff00 << " 10=" << eCoeff10 << "
      //      01=" << eCoeff01 << " 11=" << eCoeff11 << "\n";
      return {eCoeff00, eCoeff10, eCoeff01, eCoeff11};
    }
    return {-1, -1, -1, -1};
  };
  auto FindRecordArray = [&](int const &eEta, int const &eXi, double const &eX,
                             double const &eY, SingleRecInterp &eRec) -> bool {
    std::vector<double> LCoeff = FindCoefficient(eEta, eXi, eX, eY);
    if (LCoeff[0] != -1) {
      std::vector<SinglePartInterp> LPart;
      int sumMSK = 0;
      double sumWeight = 0;
      for (int i = 0; i < 4; i++) {
        int nEta = eEta + MatDir(0, i);
        int nXi = eXi + MatDir(1, i);
        int eMSK = CoordGridArr.MSK(nEta, nXi);
        sumMSK += eMSK;
        if (eMSK == 1)
          sumWeight += LCoeff[i];
      }
      if (sumMSK == 0)
        return false;
      double deltaX = eX;
      double deltaY = eY;
      for (int i = 0; i < 4; i++) {
        int nEta = eEta + MatDir(0, i);
        int nXi = eXi + MatDir(1, i);
        if (CoordGridArr.MSK(nEta, nXi) == 1) {
          double alpha = LCoeff[i] / sumWeight;
          SinglePartInterp ePart = {nEta, nXi, alpha};
          deltaX = deltaX - LCoeff[i] * CoordGridArr.LON(nEta, nXi);
          deltaY = deltaY - LCoeff[i] * CoordGridArr.LAT(nEta, nXi);
          LPart.push_back(ePart);
        }
      }
      eRec = {true, LPart};
      return true;
    }
    return false;
  };
  auto AdmissibleEtaXi = [&](int const &eEta, int const &eXi) -> bool {
    if (eEta >= 0 && eEta < nbEta - 1 && eXi >= 0 && eXi < nbXi - 1)
      return true;
    return false;
  };
  int nbExtrapolationPoint = 0;
  auto GetExtrapolationRecord = [&](double const &eX,
                                    double const &eY) -> SingleRecInterp {
    if (!AllowExtrapolation)
      return {false, {}};
    nbExtrapolationPoint++;
    // Now considering all entries
    bool IsFirst = true;
    double MinNorm;
    int iselect = -1;
    int jselect = -1;
    for (auto &eWetEntry : ListWetEntry) {
      double dx = eWetEntry.lon - eX;
      double dy = eWetEntry.lat - eY;
      double eNorm = sqrt(dx * dx + dy * dy);
      if (IsFirst) {
        IsFirst = false;
        MinNorm = eNorm;
        iselect = eWetEntry.i;
        jselect = eWetEntry.j;
      } else {
        if (eNorm < MinNorm) {
          MinNorm = eNorm;
          iselect = eWetEntry.i;
          jselect = eWetEntry.j;
        }
      }
    }
    if (IsFirst) {
      return {false, {}};
    } else {
      return {true, {{iselect, jselect, static_cast<double>(1)}}};
    }
  };
  auto FindRecord = [&](int const &eEta, int const &eXi, double const &eX,
                        double const &eY) -> SingleRecInterp {
    //    std::cerr << "FindRecord, step 1\n";
    bool testQuad = TestFeasibilityByQuad(eQuad, eX, eY);
    //    std::cerr << "FindRecord, step 1\n";
    if (!testQuad)
      return GetExtrapolationRecord(eX, eY);
    //    std::cerr << "FindRecord, step 2\n";
    auto fEvaluateCorrectness = [&](int const &fEta, int const &fXi,
                                    bool &DoSomething,
                                    SingleRecInterp &eRec) -> bool {
      if (AdmissibleEtaXi(fEta, fXi)) {
        DoSomething = true;
        bool test = FindRecordArray(fEta, fXi, eX, eY, eRec);
        if (test)
          return true;
      }
      return false;
    };
    //    std::cerr << "FindRecord, step 3\n";
    SingleRecInterp eRec;
    int sizExp = 3;
    //    std::cerr << "FindRecord, step 4\n";
    for (int siz = 0; siz < sizExp; siz++) {
      std::vector<int> ePair;
      bool DoSomething = false;
      if (siz == 0) {
        if (fEvaluateCorrectness(eEta, eXi, DoSomething, eRec))
          return eRec;
      }
      for (int i = -siz; i < siz; i++) {
        if (fEvaluateCorrectness(eEta - siz, eXi + i, DoSomething, eRec))
          return eRec;
        if (fEvaluateCorrectness(eEta + i, eXi + siz, DoSomething, eRec))
          return eRec;
        if (fEvaluateCorrectness(eEta + siz, eXi - i, DoSomething, eRec))
          return eRec;
        if (fEvaluateCorrectness(eEta - i, eXi - siz, DoSomething, eRec))
          return eRec;
      }
      if (!DoSomething)
        return GetExtrapolationRecord(eX, eY);
    }
    //    std::cerr << "FindRecord, step 5\n";
    PairLL ePt{eX, eY};
    PairCoord ePair = FindContaining(ePt, CoordGridArr.LON, CoordGridArr.LAT);
    //    std::cerr << "FindRecord, step 6\n";
    if (ePair.i == -1)
      return GetExtrapolationRecord(eX, eY);
    //    std::cerr << "FindRecord, step 7\n";
    bool test = FindRecordArray(ePair.i, ePair.j, eX, eY, eRec);
    //    std::cerr << "FindRecord, step 8\n";
    if (test)
      return eRec;
    //    std::cerr << "FindRecord, step 9\n";
    return GetExtrapolationRecord(eX, eY);
  };
  int nbPoint = ListXY.cols();
  std::cerr << "nbPoint=" << nbPoint << "\n";
  int iEtaPrev = 0;
  int iXiPrev = 0;
  std::vector<SingleRecInterp> LRec(nbPoint);
  for (int iPoint = 0; iPoint < nbPoint; iPoint++) {
    //    std::cerr << "iPoint=" << iPoint << " / " << nbPoint << "\n";
    double Xp = ListXY(0, iPoint);
    double Yp = ListXY(1, iPoint);
    //    std::cerr << "Xp=" << Xp << " Yp=" << Yp << "\n";
    SingleRecInterp eEnt = FindRecord(iEtaPrev, iXiPrev, Xp, Yp);
    /*
    std::cerr << "iPoint=" << iPoint << " eEnt=";
    for (auto & eP : eEnt.LPart)
      std::cerr << "(" << eP.eEta << "," << eP.eXi << ")";
    std::cerr << "\n";
    */
    LRec[iPoint] = eEnt;
    //    std::cerr << "iPoint=" << iPoint << " / " << nbPoint << " x/y=" << Xp
    //    << "/" << Yp << " status=" << eEnt.status << "\n";
    if (eEnt.status) {
      iEtaPrev = eEnt.LPart[0].eEta;
      iXiPrev = eEnt.LPart[0].eXi;
    }
  }
  std::cerr << "FD_FIND_ELE, nbExtrapolationPoint = " << nbExtrapolationPoint
            << "\n";
  return LRec;
}

void Print_InterpolationError(std::vector<SingleRecInterp> const &LRec,
                              GridArray const &GrdArr,
                              MyMatrix<double> const &ListXY) {
  int nbPoint = LRec.size();
  double TotalErr = 0;
  int nbPointInside = 0;
  MyMatrix<double> const &LON = GrdArr.GrdArrRho.LON;
  MyMatrix<double> const &LAT = GrdArr.GrdArrRho.LAT;
  for (int iPoint = 0; iPoint < nbPoint; iPoint++) {
    double deltaX = -ListXY(0, iPoint);
    double deltaY = -ListXY(1, iPoint);
    SingleRecInterp eSing = LRec[iPoint];
    //    std::cerr << "iPoint=" << iPoint << " / " << nbPoint << "\n";
    if (eSing.status) {
      nbPointInside++;
      for (auto &ePart : eSing.LPart) {
        int eEta = ePart.eEta;
        int eXi = ePart.eXi;
        double eCoeff = ePart.eCoeff;
        //        std::cerr << "eEta=" << eEta << " eXi=" << eXi << " |LON|=" <<
        //        LON.rows() << " / " << LON.cols() << " |LAT|=" << LAT.rows()
        //        << " / " << LAT.cols() << "\n";
        deltaX += eCoeff * LON(eEta, eXi);
        deltaY += eCoeff * LAT(eEta, eXi);
      }
      double eErr = std::abs(deltaX) + std::abs(deltaY);
      //      std::cerr << "iPoint=" << iPoint << " eErr=" << eErr << "\n";
      TotalErr = TotalErr + eErr;
    }
  }
  std::cerr << "Total Interpolation error = " << TotalErr << "\n";
  std::cerr << "nbPoint=" << nbPoint << " nbPointInside=" << nbPointInside
            << "\n";
}

QuadCoordinate Get_QuadCoordinate(GridArray const &GrdArr) {
  double MinLon = 0;
  double MaxLon = 0;
  double MinLat = 0;
  double MaxLat = 0;
  bool IsFirst = true;
  int eta_rho = GrdArr.GrdArrRho.LON.rows();
  int xi_rho = GrdArr.GrdArrRho.LON.cols();
  for (int i = 0; i < eta_rho; i++)
    for (int j = 0; j < xi_rho; j++) {
      double eLon = GrdArr.GrdArrRho.LON(i, j);
      double eLat = GrdArr.GrdArrRho.LAT(i, j);
      if (IsFirst) {
        MinLon = eLon;
        MaxLon = eLon;
        MinLat = eLat;
        MaxLat = eLat;
      } else {
        if (eLon > MaxLon)
          MaxLon = eLon;
        if (eLon < MinLon)
          MinLon = eLon;
        if (eLat > MaxLat)
          MaxLat = eLat;
        if (eLat < MinLat)
          MinLat = eLat;
      }
      IsFirst = false;
    }
  return {MinLon, MaxLon, MinLat, MaxLat};
}

std::vector<SingleRecInterp>
General_FindInterpolationWeight(GridArray const &GrdArr,
                                MyMatrix<double> const &ListXY,
                                bool const &AllowExtrapolation) {
  if (ListXY.rows() != 2) {
    std::cerr << "The number of rows needs to be equal to 2\n";
    throw TerminalException{1};
  }
  std::cerr << "General_FindInterpolationWeight : AllowExtrapolation="
            << AllowExtrapolation << "\n";
  std::vector<SingleRecInterp> LRec;
  QuadCoordinate eQuad = Get_QuadCoordinate(GrdArr);
  if (GrdArr.IsFE == 0) {
    std::cerr << "Before FD_FIND_ELE\n";
    LRec = FD_FIND_ELE(GrdArr.GrdArrRho, eQuad, ListXY, AllowExtrapolation);
  } else {
    std::cerr << "Before TRIG_FIND_ELE\n";
    LRec = TRIG_FIND_ELE(GrdArr.INE, GrdArr.GrdArrRho.LON, GrdArr.GrdArrRho.LAT,
                         eQuad, ListXY, AllowExtrapolation);
  }
  Print_InterpolationError(LRec, GrdArr, ListXY);
  return LRec;
}

#ifdef USE_OPENCV_LIBRARY
std::vector<SingleRecInterp>
NearestInterpolation_FindWeight_LonLat(MyMatrix<double> const &LON,
                                       MyMatrix<double> const &LAT,
                                       MyMatrix<double> const &ListXY) {
  cv::Mat_<double> features(0, 2);
  int siz = LON.rows();
  for (int i = 0; i < siz; i++) {
    cv::Mat row = (cv::Mat_<double>(1, 2) << LON(i, 0), LAT(i, 0));
    features.push_back(row);
  }
  cv::flann::Index flann_index(features, cv::flann::KDTreeIndexParams(1));
  int nbPoint = ListXY.cols();
  std::vector<SingleRecInterp> LRec(nbPoint);
  for (int iPoint = 0; iPoint < nbPoint; iPoint++) {
    cv::Mat indices, dists; // neither assume type nor size here !
    int knn = 1;
    double eLon = ListXY(0, iPoint);
    double eLat = ListXY(1, iPoint);
    cv::Mat query = (cv::Mat_<double>(1, 2) << eLon, eLat);
    flann_index.knnSearch(query, indices, dists, knn,
                          cv::flann::SearchParams(32));
    if (indices.rows != 1) {
      std::cerr << "The size is wrong here\n";
      throw TerminalException{1};
    }
    int eEta = indices.at<int>(0);
    SingleRecInterp eRec{true, {{eEta, 0, static_cast<double>(1)}}};
    LRec[iPoint] = eRec;
  }
  return LRec;
}

std::vector<SingleRecInterp>
NearestInterpolation_FindWeight_FD(GridArray const &GrdArr,
                                   MyMatrix<double> const &ListXY) {
  int nbWetRho = GrdArr.GrdArrRho.nbWet;
  std::vector<int> ListI(nbWetRho), ListJ(nbWetRho);
  MyMatrix<double> LON(nbWetRho, 1), LAT(nbWetRho, 1);
  int pos = 0;
  int etaRho = GrdArr.GrdArrRho.LON.rows();
  int xiRho = GrdArr.GrdArrRho.LON.cols();
  for (int i = 0; i < etaRho; i++)
    for (int j = 0; j < xiRho; j++)
      if (GrdArr.GrdArrRho.MSK(i, j) == 1) {
        ListI[pos] = i;
        ListJ[pos] = j;
        LON(pos, 0) = GrdArr.GrdArrRho.LON(i, j);
        LAT(pos, 0) = GrdArr.GrdArrRho.LAT(i, j);
        pos++;
      }
  std::vector<SingleRecInterp> LRec =
      NearestInterpolation_FindWeight_LonLat(LON, LAT, ListXY);
  std::vector<SingleRecInterp> LRecReturn;
  for (auto &eRec : LRec) {
    std::vector<SinglePartInterp> LPartNew;
    for (auto &ePart : eRec.LPart) {
      int eEta = ePart.eEta;
      int eEtaNew = ListI[eEta];
      int eXiNew = ListJ[eEta];
      SinglePartInterp ePartNew{eEtaNew, eXiNew, ePart.eCoeff};
      LPartNew.push_back(ePartNew);
    }
    SingleRecInterp eRecNew{eRec.status, LPartNew};
    LRecReturn.push_back(eRecNew);
  }
  return LRecReturn;
}

std::vector<SingleRecInterp>
NearestInterpolation_FindWeight_FE(GridArray const &GrdArr,
                                   MyMatrix<double> const &ListXY) {
  return NearestInterpolation_FindWeight_LonLat(GrdArr.GrdArrRho.LON,
                                                GrdArr.GrdArrRho.LAT, ListXY);
}

std::vector<SingleRecInterp>
NearestInterpolation_FindWeight(GridArray const &GrdArr,
                                MyMatrix<double> const &ListXY) {
  if (GrdArr.IsFE == 0) {
    return NearestInterpolation_FindWeight_FD(GrdArr, ListXY);
  } else {
    return NearestInterpolation_FindWeight_FE(GrdArr, ListXY);
  }
}
#endif

// clang-format off
#endif  // SRC_OCEAN_INTERPOLATION_H_
// clang-format on
