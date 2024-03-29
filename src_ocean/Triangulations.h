// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_OCEAN_TRIANGULATIONS_H_
#define SRC_OCEAN_TRIANGULATIONS_H_

#include "Basic_Ocean_types.h"
#include "GRAPH_GraphicalBasic.h"
#include "SphericalGeom.h"
#include <limits>
#include <set>
#include <vector>

GraphSparseImmutable GetUnstructuredVertexAdjInfo(MyMatrix<int> const &INE,
                                                  size_t nbNode) {
  size_t nbEle = INE.rows();
  std::vector<size_t> ListNbEnt(nbNode, 0);
  for (size_t iEle = 0; iEle < nbEle; iEle++)
    for (size_t i = 0; i < 3; i++) {
      size_t eVert = INE(iEle, i);
      ListNbEnt[eVert] += 2;
    }
  std::cerr << "nbNode=" << nbNode << "\n";
  size_t TotalSum_unrefined = 6 * nbEle;
  size_t miss_val = std::numeric_limits<size_t>::max();
  std::vector<size_t> ListStart_unrefined(nbNode + 1);
  ListStart_unrefined[0] = 0;
  for (size_t iNode = 0; iNode < nbNode; iNode++)
    ListStart_unrefined[iNode + 1] =
        ListStart_unrefined[iNode] + ListNbEnt[iNode];
  std::vector<size_t> ListListAdj_unrefined(TotalSum_unrefined, miss_val);
  std::vector<size_t> ListIndexPos(nbNode, 0);
  auto fInsert = [&](size_t const &eVert, size_t const &eVertAdj) -> void {
    size_t eStart = ListStart_unrefined[eVert];
    size_t eEnd = ListStart_unrefined[eVert] + ListIndexPos[eVert];
    if (ListListAdj_unrefined[eEnd] != miss_val) {
      std::cerr << "Logical error in the code\n";
      throw TerminalException{1};
    }
    for (size_t i = eStart; i < eEnd; i++)
      if (ListListAdj_unrefined[i] == eVertAdj)
        return;
    ListIndexPos[eVert]++;
    ListListAdj_unrefined[eEnd] = eVertAdj;
  };
  for (size_t iEle = 0; iEle < nbEle; iEle++)
    for (int i = 0; i < 3; i++) {
      int iNext = NextIdx(3, i);
      int iPrev = PrevIdx(3, i);
      size_t eVert = INE(iEle, i);
      size_t eVertP = INE(iEle, iPrev);
      size_t eVertN = INE(iEle, iNext);
      fInsert(eVert, eVertP);
      fInsert(eVert, eVertN);
    }
  std::vector<size_t> ListStart(nbNode + 1, 0);
  for (size_t iNode = 0; iNode < nbNode; iNode++)
    ListStart[iNode + 1] = ListStart[iNode] + ListIndexPos[iNode];
  size_t TotalSum = ListStart[nbNode];
  std::vector<size_t> ListListAdj(TotalSum);
  for (size_t iNode = 0; iNode < nbNode; iNode++) {
    size_t eStart = ListStart[iNode];
    size_t eStart_unrefined = ListStart_unrefined[iNode];
    size_t siz = ListIndexPos[iNode];
    for (size_t i = 0; i < siz; i++)
      ListListAdj[eStart + i] = ListListAdj_unrefined[eStart_unrefined + i];
  }
  return GraphSparseImmutable(nbNode, ListStart, ListListAdj);
}

std::vector<int> GetListNeighbor(MyMatrix<int> const &INE, int nbNode) {
  std::vector<int> Status(nbNode, 0), Neighbor(nbNode, -1);
  std::vector<int> PrevVert(nbNode), NextVert(nbNode), Collected(nbNode);
  int mne = INE.rows();
  int sizMesh = INE.cols();
  for (int ie = 0; ie < mne; ie++) {
    for (int i = 0; i < sizMesh; i++) {
      int iPrev = PrevIdx(sizMesh, i);
      int iNext = NextIdx(sizMesh, i);
      int ip = INE(ie, i);
      int ipnext = INE(ie, iNext);
      int ipprev = INE(ie, iPrev);
      if (Status[ip] == 0) {
        Status[ip] = 1;
        PrevVert[ip] = ipprev;
        NextVert[ip] = ipnext;
      }
    }
  }
  for (int iNode=0; iNode<nbNode; iNode++) {
    if (Status[iNode] == 0) {
      std::cerr << "Node iNode=" << iNode << " was not matched  by the process\n";
      std::cerr << "This is possible if the node is not contained in any triangle\n";
      throw TerminalException{1};
    }
  }
  for (int i = 0; i < nbNode; i++)
    Status[i] = 0;
  while (true) {
    for (int i = 0; i < nbNode; i++)
      Collected[i] = 0;
    // The principle is that if a node is not saturated,
    // then we can do some increment.
    for (int ie = 0; ie < mne; ie++) {
      for (int i = 0; i < sizMesh; i++) {
        int iPrev = PrevIdx(sizMesh, i);
        int iNext = NextIdx(sizMesh, i);
        int ip = INE(ie, i);
        int ipnext = INE(ie, iNext);
        int ipprev = INE(ie, iPrev);
        if (Status[ip] == 0) {
          int zNext = NextVert[ip];
          if (zNext == ipprev) {
            Collected[ip] = 1;
            NextVert[ip] = ipnext;
            if (NextVert[ip] == PrevVert[ip])
              Status[ip] = 1;
          }
        }
      }
    }
    bool IsFinished = true;
    for (int i = 0; i < nbNode; i++) {
      if (Collected[i] == 0 && Status[i] == 0) {
        Status[i] = -1;
        Neighbor[i] = NextVert[i];
      }
      if (Status[i] == 0)
        IsFinished = false;
    }
    if (IsFinished)
      break;
  }
  return Neighbor;
}



// Get the points that are contained in more than 1 triangle and for which the
// structure is discontinuous
std::vector<int> GetPointsInDiscontinuedTriangles(MyMatrix<int> const &INE, int nbNode) {
  std::vector<int> ListNbTriangleContained(nbNode,0);
  int mne = INE.rows();
  for (int ie=0; ie<mne; ie++) {
    for (int i=0; i<3; i++) {
      int ip = INE(ie,i);
      ListNbTriangleContained[ip]++;
    }
  }
  int n_max = 0;
  for (int ip=0; ip<nbNode; ip++) {
    if (ListNbTriangleContained[ip] > n_max)
      n_max = ListNbTriangleContained[ip];
  }
  std::vector<int> TotalAtt(nbNode,0);
  int upper_lim = 2 * n_max;
  MyMatrix<int> MatNeighbor = ZeroMatrix<int>(nbNode, upper_lim);
  MyMatrix<int> MatAtt = ZeroMatrix<int>(nbNode, upper_lim);
  auto f_insert=[&](int const& ip, int const& jp) -> void {
    for (int u=0; u<TotalAtt[ip]; u++) {
      if (MatNeighbor(ip, u) == jp) {
        MatAtt(ip, u) += 1;
        return;
      }
    }
    MatNeighbor(ip, TotalAtt[ip]) = jp;
    MatAtt(ip, TotalAtt[ip]) = 1;
    TotalAtt[ip] += 1;
  };
  for (int ie=0; ie<mne; ie++) {
    for (int i=0; i<3; i++) {
      int ip = INE(ie, i);
      for (int j=0; j<3; j++) {
        if (i != j) {
          int jp = INE(ie, j);
          f_insert(ip, jp);
        }
      }
    }
  }
  if (false) {
    int TestNode = 125755;
    std::map<int,size_t> map_test;
    for (int ie=0; ie<mne; ie++) {
      bool IsMatching = false;
      for (int i=0; i<3; i++) {
        if (INE(ie,i) == TestNode) {
          IsMatching = true;
        }
      }
      if (IsMatching) {
        std::cerr << "trig ie=" << ie << " INE =";
        for (int i=0; i<3; i++) {
          int ePt = INE(ie,i);
          std::cerr << " " << ePt;
          if (ePt != TestNode)
            map_test[ePt]++;
        }
        std::cerr << "\n";
      }
    }
    for (auto & kv : map_test) {
      std::cerr << "ePt=" << kv.first << " mult=" << kv.second << "\n";
    }
  }
  int n_inner = 0;
  int n_bnd = 0;
  std::ofstream os("ListBnd1");
  auto status_discont=[&](int ip) -> int {
    int n1 = 0;
    for (int u=0; u<TotalAtt[ip]; u++) {
      if (MatAtt(ip, u) == 1)
        n1++;
    }
    if (n1 == 0)
      n_inner++;
    if (n1 == 2) {
      n_bnd++;
      os << ip << "\n";
    }
    if (n1 == 0 || n1 == 2)
      return 1;
    return 0;
  };
  std::vector<int> ListStatus(nbNode);
  int n_discont = 0;
  for (int ip=0; ip<nbNode; ip++) {
    int status = status_discont(ip);
    ListStatus[ip] = status;
    if (status == 0)
      n_discont++;
  }
  std::cerr << "n_discont=" << n_discont << " n_inner=" << n_inner << " n_bnd=" << n_bnd << "\n";
  return ListStatus;
}

// Computation of the boundary status
// 0: should not happen in return
// 1: should be inside node
// -1: should be boundary node
MyVector<int> GetBoundaryStatus(MyMatrix<int> const &INE, int nbNode) {
  std::vector<int> Neighbor = GetListNeighbor(INE, nbNode);
  MyVector<int> Status(nbNode);
  for (int iNode = 0; iNode < nbNode; iNode++) {
    if (Neighbor[iNode] == -1) {
      Status(iNode) = 1;
    } else {
      Status(iNode) = -1;
    }
  }
  return Status;
}

// There is a more advanced version of this code that handle
// pathological cases in WWM code in wwm_netcdf.F90 in the section
// See also interpol.F90 in polymesh code.
//
// Here we use the simpler code.
std::vector<std::vector<int>> GetListBoundaryCycles(MyMatrix<int> const &INE,
                                                    int nbNode) {
  std::vector<int> Neighbor = GetListNeighbor(INE, nbNode);
  int nbBound = 0;
  std::ofstream os("ListBnd2");
  for (int iVert = 0; iVert < nbNode; iVert++)
    if (Neighbor[iVert] >= 0) {
      nbBound++;
      os << iVert << "\n";
    }
  std::vector<int> ListMap(nbNode, -1);
  std::vector<int> ListMapRev(nbBound);
  std::cerr << "nbNode=" << nbNode << " nbBound=" << nbBound << "\n";
  int idx = 0;
  for (int iVert = 0; iVert < nbNode; iVert++)
    if (Neighbor[iVert] >= 0) {
      ListMap[iVert] = idx;
      ListMapRev[idx] = iVert;
      idx++;
    }
  std::vector<int> NeighborRed(nbBound);
  for (int iB = 0; iB < nbBound; iB++) {
    int iVert = ListMapRev[iB];
    int iVertImg = Neighbor[iVert];
    int iBimg = ListMap[iVertImg];
    NeighborRed[iB] = iBimg;
  }
  std::vector<int> Status(nbBound, -1);
  auto GetOneUnset = [&]() -> int {
    for (int iB = 0; iB < nbBound; iB++)
      if (Status[iB] == -1)
        return iB;
    return -1;
  };
  std::vector<std::vector<int>> eListList;
  while (true) {
    int eFirst = GetOneUnset();
    if (eFirst == -1)
      break;
    std::vector<int> eList;
    std::vector<int> eList_idxB;
    auto f_terminate=[&](std::string const& strerror) -> void {
      std::cerr << "Error: " << strerror << "\n";
      std::map<int,size_t> map;
      for (auto & eEnt : eList) {
        map[eEnt]++;
      }
      std::map<size_t,size_t> mult;
      for (auto & kv : map) {
        mult[kv.second]++;
      }
      for (auto & kv : mult) {
        std::cerr << "multiplicity " << kv.first << " attained " << kv.second << " times\n";
      }
      std::cerr << "nbNode=" << nbNode << " nbBound=" << nbBound << "\n";
      throw TerminalException{1};
    };
    int idxB = eFirst;
    while (true) {
      eList_idxB.push_back(idxB);
      int eVert = ListMapRev[idxB];
      eList.push_back(eVert);
      if (int(eList.size()) > nbNode) {
        f_terminate("eList is clearly too large");
      }
      if (Status[idxB] == 0) {
        f_terminate("already passed point");
      }
      Status[idxB] = 0;
      idxB = NeighborRed[idxB];
      if (idxB == eFirst)
        break;
    }
    eListList.push_back(eList);
  }
  return eListList;
}

MyMatrix<int> GetEdgeSet(MyMatrix<int> const &INE, int nbNode) {
  std::vector<int> ListDegree(nbNode, 0);
  int mne = INE.rows();
  for (int ie = 0; ie < mne; ie++)
    for (int i = 0; i < 3; i++) {
      int ip = INE(ie, i);
      ListDegree[ip] += 2;
    }
  int MaxDeg = 0;
  for (int iNode = 0; iNode < nbNode; iNode++) {
    int eDeg = ListDegree[iNode];
    if (eDeg > MaxDeg)
      MaxDeg = eDeg;
  }
  for (int iNode = 0; iNode < nbNode; iNode++)
    ListDegree[iNode] = 0;
  MyMatrix<int> ListAdjacency(nbNode, MaxDeg);
  for (int ie = 0; ie < mne; ie++) {
    int i1 = INE(ie, 0);
    int i2 = INE(ie, 1);
    int i3 = INE(ie, 2);
    int eDeg1 = ListDegree[i1];
    int eDeg2 = ListDegree[i2];
    int eDeg3 = ListDegree[i3];
    ListAdjacency(i1, eDeg1) = i2;
    ListAdjacency(i1, eDeg1 + 1) = i3;
    ListAdjacency(i2, eDeg2) = i1;
    ListAdjacency(i2, eDeg2 + 1) = i3;
    ListAdjacency(i3, eDeg3) = i1;
    ListAdjacency(i3, eDeg3 + 1) = i2;
    ListDegree[i1] = eDeg1 + 2;
    ListDegree[i2] = eDeg2 + 2;
    ListDegree[i3] = eDeg3 + 2;
  }
  int nbEdge = 0;
  for (int iNode = 0; iNode < nbNode; iNode++) {
    std::set<int> eSet;
    int eDeg = ListDegree[iNode];
    for (int iAdj = 0; iAdj < eDeg; iAdj++) {
      int eAdj = ListAdjacency(iNode, iAdj);
      if (eAdj > iNode)
        eSet.insert(eAdj);
    }
    nbEdge += eSet.size();
  }
  MyMatrix<int> ListEdges(nbEdge, 2);
  int iEdge = 0;
  for (int iNode = 0; iNode < nbNode; iNode++) {
    std::set<int> eSet;
    int eDeg = ListDegree[iNode];
    for (int iAdj = 0; iAdj < eDeg; iAdj++) {
      int eAdj = ListAdjacency(iNode, iAdj);
      if (eAdj > iNode)
        eSet.insert(eAdj);
    }
    for (auto eAdj : eSet) {
      ListEdges(iEdge, 0) = iNode;
      ListEdges(iEdge, 1) = eAdj;
      iEdge++;
    }
  }
  return ListEdges;
}

std::set<std::pair<int,int>> GetEdgeSet_set(MyMatrix<int> const &INE, int nbNode) {
  MyMatrix<int> TheMat = GetEdgeSet(INE, nbNode);
  std::set<std::pair<int,int>> TheSet;
  int nbRow = TheMat.rows();
  for (int iRow=0; iRow<nbRow; iRow++) {
    int eVert = TheMat(iRow,0);
    int fVert = TheMat(iRow,1);
    std::pair<int,int> epair = {eVert,fVert};
    TheSet.insert(epair);
  }
  return TheSet;
}



MyMatrix<int> GetFaceEdgeConnectivity(int const &nbNode,
                                      MyMatrix<int> const &LEdge,
                                      MyMatrix<int> const &INE) {
  int nbEdge = LEdge.rows();
  std::vector<int> eVect(2 * nbEdge, 0);
  std::vector<int> eVectNB(nbNode, 0);
  for (int iEdge = 0; iEdge < nbEdge; iEdge++) {
    for (int i = 0; i < 2; i++) {
      int ePt = LEdge(iEdge, i);
      eVectNB[ePt]++;
    }
  }
  std::vector<int> ListShift(nbNode, 0);
  for (int iNode = 1; iNode < nbNode; iNode++)
    ListShift[iNode] = ListShift[iNode - 1] + eVectNB[iNode - 1];
  int NbIncidence = ListShift[nbNode - 1] + eVectNB[nbNode - 1];
  MyMatrix<int> PairInfo(NbIncidence, 2);
  std::vector<int> eVectNB_att(nbNode, 0);
  for (int iEdge = 0; iEdge < nbEdge; iEdge++) {
    for (int i = 0; i < 2; i++) {
      int j = 1 - i;
      int ePt1 = LEdge(iEdge, i);
      int ePt2 = LEdge(iEdge, j);
      int eShift = ListShift[ePt1];
      int pos = eVectNB_att[ePt1];
      PairInfo(eShift + pos, 0) = ePt2;
      PairInfo(eShift + pos, 1) = iEdge;
      eVectNB_att[ePt1] = pos + 1;
    }
  }
  auto GetIedge = [&](int const &ePt1, int const &ePt2) -> int {
    int eShift = ListShift[ePt1];
    int siz = eVectNB[ePt1];
    for (int i = 0; i < siz; i++) {
      if (PairInfo(eShift + i, 0) == ePt2)
        return PairInfo(eShift + i, 1);
    }
    return -1; // should not happen at all
  };
  MyMatrix<int> FaceEdgeConn(nbEdge, 2);
  std::vector<int> nbEdgeAtt(nbEdge, 0);
  int nbFace = INE.rows();
  for (int iFace = 0; iFace < nbFace; iFace++) {
    for (int i = 0; i < 3; i++) {
      int j = (i + 1) % 3;
      int ePt1 = INE(iFace, i);
      int ePt2 = INE(iFace, j);
      int iEdge = GetIedge(ePt1, ePt2);
      int pos = nbEdgeAtt[iEdge];
      FaceEdgeConn(iEdge, pos) = iFace;
      nbEdgeAtt[iEdge] = pos + 1;
    }
  }
  for (int iEdge = 0; iEdge < nbEdge; iEdge++) {
    for (int pos = nbEdgeAtt[iEdge] + 1; pos < 2; pos++)
      FaceEdgeConn(iEdge, pos) = -999;
  }
  return FaceEdgeConn;
}

std::vector<int>
GetUnstructuredTriangleAdjInfo_vectint(MyMatrix<int> const &INE, int nbNode) {
  std::vector<int> ListDegree(nbNode, 0);
  int mne = INE.rows();
  for (int ie = 0; ie < mne; ie++)
    for (int i = 0; i < 3; i++) {
      int ip = INE(ie, i);
      ListDegree[ip]++;
    }
  int MaxDeg = 0;
  for (int iNode = 0; iNode < nbNode; iNode++) {
    int eDeg = ListDegree[iNode];
    if (eDeg > MaxDeg)
      MaxDeg = eDeg;
  }
  for (int iNode = 0; iNode < nbNode; iNode++)
    ListDegree[iNode] = 0;
  MyMatrix<int> IncidenceTrigVert(nbNode, MaxDeg);
  for (int ie = 0; ie < mne; ie++) {
    for (int i = 0; i < 3; i++) {
      int i1 = INE(ie, i);
      int eDeg = ListDegree[i1];
      IncidenceTrigVert(i1, eDeg) = ie;
      ListDegree[i1] = eDeg + 1;
    }
  }
  std::vector<int> ListAdj(3 * mne);
  auto HasSpecifiedVertex = [&](int const &iTrig, int const &eVert) -> bool {
    for (int i = 0; i < 3; i++)
      if (INE(iTrig, i) == eVert)
        return true;
    return false;
  };
  for (int ie = 0; ie < mne; ie++)
    for (int i = 0; i < 3; i++) {
      int iNext = NextIdx(3, i);
      int ip1 = INE(ie, i);
      int ip2 = INE(ie, iNext);
      int eAdj = -1;
      for (int iTrig = 0; iTrig < ListDegree[ip1]; iTrig++) {
        int eTrig = IncidenceTrigVert(ip1, iTrig);
        if (HasSpecifiedVertex(eTrig, ip2) && eTrig != ie)
          eAdj = eTrig;
      }
      ListAdj[3 * ie + i] = eAdj;
    }
  return ListAdj;
}

std::vector<double> GetVectorAreas(GridArray const &GrdArr) {
  int mne = GrdArr.INE.rows();
  std::vector<double> ListArea(mne);
  for (int ie = 0; ie < mne; ie++) {
    int i1 = GrdArr.INE(ie, 0);
    int i2 = GrdArr.INE(ie, 1);
    int i3 = GrdArr.INE(ie, 2);
    double eLon1 = GrdArr.GrdArrRho.LON(i1, 0);
    double eLon2 = GrdArr.GrdArrRho.LON(i2, 0);
    double eLon3 = GrdArr.GrdArrRho.LON(i3, 0);
    double eLat1 = GrdArr.GrdArrRho.LAT(i1, 0);
    double eLat2 = GrdArr.GrdArrRho.LAT(i2, 0);
    double eLat3 = GrdArr.GrdArrRho.LAT(i3, 0);
    double deltaLON12 = eLon2 - eLon1;
    double deltaLAT12 = eLat2 - eLat1;
    double deltaLON13 = eLon3 - eLon1;
    double deltaLAT13 = eLat3 - eLat1;
    double eArea = deltaLON13 * deltaLAT12 - deltaLON12 * deltaLAT13;
    ListArea[ie] = eArea;
  }
  return ListArea;
}

void CHECK_COORDINATE_ORIENTATION(GridArray const &GrdArr) {
  int mnp = GrdArr.GrdArrRho.LON.rows();
  int mne = GrdArr.INE.rows();
  int nbPlus = 0;
  int nbMinus = 0;
  std::vector<double> ListArea = GetVectorAreas(GrdArr);
  for (int ie = 0; ie < mne; ie++) {
    double eArea = ListArea[ie];
    if (eArea > 0) {
      nbPlus++;
    } else {
      nbMinus++;
    }
    for (int i = 0; i < 3; i++) {
      int IP = GrdArr.INE(ie, i);
      if (IP < 0 || IP >= mnp) {
        std::cerr << "Error in the unstructured grid\n";
        std::cerr << "mnp=" << mnp << "  mne=" << mne << "\n";
        std::cerr << "ie=" << ie << "\n";
        std::cerr << "INE=[" << GrdArr.INE(ie, 0) << " , " << GrdArr.INE(ie, 1)
                  << " , " << GrdArr.INE(ie, 2) << "]\n";
        throw TerminalException{1};
      }
    }
  }
  if (nbPlus > 0 && nbMinus > 0) {
    std::cerr << "Orientation error\n";
    std::cerr << "nbPlus =" << nbPlus << "\n";
    std::cerr << "nbMinus=" << nbMinus << "\n";
    throw TerminalException{1};
  }
  std::cerr << "nbPlus = " << nbPlus << "  nbMinus = " << nbMinus << "\n";
}

void CHECK_UnstructuredGrid(GridArray const &GrdArr) {
  int mnp = GrdArr.GrdArrRho.LON.rows();
  int mne = GrdArr.INE.rows();
  std::cerr << "mne=" << mne << "\n";
  CHECK_COORDINATE_ORIENTATION(GrdArr);
  MyVector<int> Status = GetBoundaryStatus(GrdArr.INE, mnp);
  int nbStatusNormal = 0;
  int nbStatusBound = 0;
  for (int i = 0; i < mnp; i++) {
    if (Status[i] == 1)
      nbStatusNormal++;
    if (Status[i] == -1)
      nbStatusBound++;
  }
  std::cerr << "nbStatusNormal=" << nbStatusNormal
            << " nbStatusBound=" << nbStatusBound << "\n";
  double MinDistKM = 10000000;
  double MaxDistKM = 0;
  double MinAreaSqrKM = 10000000000000;
  double MaxAreaSqrKM = 0;
  auto GetPairLL = [&](int pt) -> PairLL {
    return {GrdArr.GrdArrRho.LON(pt), GrdArr.GrdArrRho.LAT(pt)};
  };
  for (int ie = 0; ie < mne; ie++) {
    for (int i1 = 0; i1 < 3; i1++) {
      int i2;
      if (i1 == 2) {
        i2 = 0;
      } else {
        i2 = i1 + 1;
      }
      int pt1 = GrdArr.INE(ie, i1);
      int pt2 = GrdArr.INE(ie, i2);
      PairLL pair1 = GetPairLL(pt1);
      PairLL pair2 = GetPairLL(pt2);
      double eDistKM = GeodesicDistanceKM_pair(pair1, pair2);
      if (eDistKM < MinDistKM)
        MinDistKM = eDistKM;
      if (eDistKM > MaxDistKM)
        MaxDistKM = eDistKM;
    }
    int pt1 = GrdArr.INE(ie, 0);
    int pt2 = GrdArr.INE(ie, 1);
    int pt3 = GrdArr.INE(ie, 2);
    PairLL pair1 = GetPairLL(pt1);
    PairLL pair2 = GetPairLL(pt2);
    PairLL pair3 = GetPairLL(pt3);
    double areaKM = SphericalCoordinateAreaKM_pair(pair1, pair2, pair3);
    if (areaKM < MinAreaSqrKM)
      MinAreaSqrKM = areaKM;
    if (areaKM > MaxAreaSqrKM)
      MaxAreaSqrKM = areaKM;
  }
  std::cerr << "dist km  : min=" << MinDistKM << " max=" << MaxDistKM << "\n";
  std::cerr << "area km2 : min=" << MinAreaSqrKM << " max=" << MaxAreaSqrKM
            << "\n";
  //
  (void)GetPointsInDiscontinuedTriangles(GrdArr.INE, mnp);

  std::vector<std::vector<int>> ListConn =
      GetListBoundaryCycles(GrdArr.INE, mnp);
  int nbConn = ListConn.size();
  int nbIsland = nbConn - 1;
  std::cerr << "number of island resolved = " << nbIsland << "\n";
}

std::vector<int> GetLongestCycle(std::vector<std::vector<int>> const& ListCyc) {
  int iCycSel = -1;
  size_t maxSize = 0;
  int n_cycle = ListCyc.size();
  for (int iCyc=0; iCyc<n_cycle; iCyc++) {
    size_t len = ListCyc[iCyc].size();
    if (len > maxSize) {
      maxSize = len;
      iCycSel = iCyc;
    }
  }
  return ListCyc[iCycSel];
}

std::vector<int> GetShortestSegment(std::vector<int> const& eCyc, int const& pos0, int const& pos1) {
  auto get_index=[&](int pos) -> int {
    int len = eCyc.size();
    for (int u=0; u<len; u++)
      if (eCyc[u] == pos)
        return u;
    return -1;
  };
  int len = eCyc.size();
  std::cerr << "eCyc =";
  for (auto & eVal : eCyc)
    std::cerr << " " << eVal;
  std::cerr << "\n";
  auto get_segment=[&](int const& pos_start, int const& pos_end) -> std::vector<int> {
    std::vector<int> the_list{eCyc[pos_start]};
    int pos = pos_start;
    while(true) {
      int pos_next;
      if (pos == len - 1) {
        pos_next = 0;
      } else {
        pos_next = pos + 1;
      }
      pos = pos_next;
      the_list.push_back(eCyc[pos]);
      if (pos == pos_end)
        break;
    }
    return the_list;
  };
  int index0 = get_index(pos0);
  int index1 = get_index(pos1);
  std::cerr << "index0=" << index0 << " index1=" << index1 << "\n";
  std::vector<int> segmentA = get_segment(index0, index1);
  std::vector<int> segmentB = get_segment(index1, index0);
  if (segmentA.size() < segmentB.size()) {
    return segmentA;
  } else {
    return segmentB;
  }
}



void CHECK_CombinatorialGrid(GridArray const &GrdArr) {
  int mnp = GrdArr.GrdArrRho.LON.rows();
  int mne = GrdArr.INE.rows();
  std::vector<int> CCON(mnp, 0);
  int POS_TRICK[3][2];
  POS_TRICK[0][0] = 1;
  POS_TRICK[1][0] = 2;
  POS_TRICK[2][0] = 0;
  POS_TRICK[0][1] = 2;
  POS_TRICK[1][1] = 0;
  POS_TRICK[2][1] = 1;
  for (int ie = 0; ie < mne; ie++)
    for (int i = 0; i < 3; i++) {
      int ip = GrdArr.INE(ie, i);
      CCON[ip]++;
    }
  int MAXMNECON = 0;
  for (int ip = 0; ip < mnp; ip++) {
    int eCon = CCON[ip];
    if (eCon > MAXMNECON)
      MAXMNECON = eCon;
  }
  std::vector<int> CHILF(mnp, 0);
  Eigen::Tensor<int, 3> CELLVERTEX(mnp, MAXMNECON, 2);
  for (int ie = 0; ie < mne; ie++)
    for (int j = 0; j < 3; j++) {
      int i = GrdArr.INE(ie, j);
      CELLVERTEX(i, CHILF[i], 0) = ie;
      CELLVERTEX(i, CHILF[i], 1) = j;
      CHILF[i]++;
    }
  MyMatrix<int> IE_CELL2(mnp, MAXMNECON);
  MyMatrix<int> POS_CELL2(mnp, MAXMNECON);
  for (int ip = 0; ip < mnp; ip++)
    for (int i = 0; i < CCON[ip]; i++) {
      IE_CELL2(ip, i) = CELLVERTEX(ip, i, 0);
      POS_CELL2(ip, i) = CELLVERTEX(ip, i, 1);
    }
  for (int ie = 0; ie < mne; ie++)
    for (int i = 0; i < 3; i++) {
      int INEXT = POS_TRICK[i][0];
      int ip = GrdArr.INE(ie, i);
      int IP_NEXT = GrdArr.INE(ie, INEXT);
      int nbMatch = 0;
      std::vector<int> Lmatch;
      for (int icon = 0; icon < CCON[ip]; icon++) {
        int ie2 = IE_CELL2(ip, icon);
        if (ie != ie2) {
          int POS = POS_CELL2(ip, icon);
          int POS_NEXT = POS_TRICK[POS][0];
          int IP_ADJ_NEXT = GrdArr.INE(ie2, POS_NEXT);
          if (IP_ADJ_NEXT == IP_NEXT) {
            std::cerr << "Combinatorial orientability problem\n";
            std::cerr << "IE=" << ie << " IE2=" << ie2 << "\n";
            std::cerr << "IP=" << ip << " IP_NEXT=" << IP_NEXT << "\n";
            throw TerminalException{1};
          }
          int POS_PREV = POS_TRICK[POS][1];
          int IP_ADJ_PREV = GrdArr.INE(ie2, POS_PREV);
          if (IP_ADJ_PREV == IP_NEXT) {
            nbMatch++;
            Lmatch.push_back(ie2);
          }
        }
      }
      if (nbMatch > 1) {
        std::cerr << "nbMatch is too large.\n";
        std::cerr << "Should be 0 for boundary edge\n";
        std::cerr << "Should be 1 for interior edges\n";
        std::cerr << "ie=" << ie << " i=" << i << "\n";
        std::cerr << " ip=" << ip << " IP_NEXT=" << IP_NEXT << "\n";
        std::cerr << "nbMatch=" << nbMatch << "\n";
        for (int iMatch = 0; iMatch < nbMatch; iMatch++) {
          int iem = Lmatch[iMatch];
          std::cerr << "  iMatch=" << iMatch << " ie=" << iem << "\n";
          std::cerr << "     ine=[" << GrdArr.INE(iem, 0) << ","
                    << GrdArr.INE(iem, 1) << "," << GrdArr.INE(iem, 2) << "]\n";
        }
        throw TerminalException{1};
      }
    }
  std::cerr << "Now leaving the combinatorial check\n";
}

GridArray MergeNeighboringVertices(GridArray const& GrdArr, double const& CritDistKM) {
  std::unordered_set<std::pair<int,int>> SetEdges;
  auto f_insert=[&](int const& u, int const& v) -> void {
    if (u < v)
      SetEdges.insert({u,v});
    else
      SetEdges.insert({v,u});
  };
  int mne = GrdArr.INE.rows();
  int mnp = GrdArr.GrdArrRho.LON.rows();
  double MinDistKM = std::numeric_limits<double>::max();
  for (int ie=0; ie<mne; ie++) {
    for (int i=0; i<3; i++) {
      int j = (i + 1) % 3;
      int ip = GrdArr.INE(ie, i);
      int jp = GrdArr.INE(ie, j);
      double eLon = GrdArr.GrdArrRho.LON(ip,0);
      double eLat = GrdArr.GrdArrRho.LAT(ip,0);
      double fLon = GrdArr.GrdArrRho.LON(jp,0);
      double fLat = GrdArr.GrdArrRho.LAT(jp,0);
      double distKM = GeodesicDistanceKM(eLon, eLat, fLon, fLat);
      if (distKM < CritDistKM)
        f_insert(ip, jp);
      if (distKM < MinDistKM)
        MinDistKM = distKM;
    }
  }
  int n_edge = SetEdges.size();
  std::cerr << "n_edge=" << n_edge << "\n";
  MyMatrix<size_t> ListEdge(n_edge,2);
  int pos = 0;
  for (auto & eEdge : SetEdges) {
    ListEdge(pos,0) = eEdge.first;
    ListEdge(pos,1) = eEdge.second;
    pos++;
  }
  GraphListAdj GR(ListEdge, mnp);
  std::vector<size_t> ListStatus = ConnectedComponents_vector(GR);
  int nbConn = VectorMax(ListStatus) + 1;
  std::cerr << "nbConn=" << nbConn << "\n";
  int mne_red = 0;
  for (int ie=0; ie<mne; ie++) {
    int i0 = GrdArr.INE(ie,0);
    int i1 = GrdArr.INE(ie,1);
    int i2 = GrdArr.INE(ie,2);
    int stat0 = ListStatus[i0];
    int stat1 = ListStatus[i1];
    int stat2 = ListStatus[i2];
    if (stat0 != stat1 && stat0 != stat2 && stat1 != stat2) {
      mne_red++;
    }
  }
  std::cerr << "mne=" << mne << " mne_red=" << mne_red << "\n";
  MyMatrix<int> INEred(mne_red,3);
  pos = 0;
  for (int ie=0; ie<mne; ie++) {
    int i0 = GrdArr.INE(ie,0);
    int i1 = GrdArr.INE(ie,1);
    int i2 = GrdArr.INE(ie,2);
    int stat0 = ListStatus[i0];
    int stat1 = ListStatus[i1];
    int stat2 = ListStatus[i2];
    if (stat0 != stat1 && stat0 != stat2 && stat1 != stat2) {
      INEred(pos,0) = stat0;
      INEred(pos,1) = stat1;
      INEred(pos,2) = stat2;
      pos++;
    }
  }
  MyMatrix<double> const& DEP = *GrdArr.GrdArrRho.DEP;
  MyMatrix<double> LONred(nbConn,1), LATred(nbConn,1), DEPred(nbConn,1);
  std::vector<size_t> Noccur(nbConn,0);
  for (int i=0; i<mnp; i++) {
    int stat = ListStatus[i];
    Noccur[stat]++;
    LONred(stat,0) += GrdArr.GrdArrRho.LON(i,0);
    LATred(stat,0) += GrdArr.GrdArrRho.LAT(i,0);
    DEPred(stat,0) += DEP(i,0);
  }
  for (int i_conn=0; i_conn < nbConn; i_conn++) {
    double div = 1.0 / static_cast<double>(Noccur[i_conn]);
    LONred(i_conn,0) *= div;
    LATred(i_conn,0) *= div;
    DEPred(i_conn,0) *= div;
  }
  GridArray GrdArrRet = GrdArr;
  GrdArrRet.INE = INEred;
  GrdArrRet.GrdArrRho.LON = LONred;
  GrdArrRet.GrdArrRho.LAT = LATred;
  GrdArrRet.GrdArrRho.DEP = DEPred;
  CHECK_UnstructuredGrid(GrdArrRet);
  return GrdArrRet;
}

GridArray SelectSubsetVertices(GridArray const& GrdArr, std::vector<int> const& ListStatus) {
  // The vertices whose value in ListStatus is 1 are kept.
  int mnp = GrdArr.GrdArrRho.LON.rows();
  int mne = GrdArr.INE.rows();
  int mnp_red = 0;
  for (int iP=0; iP<mnp; iP++) {
    if (ListStatus[iP] == 1)
      mnp_red += 1;
  }
  std::cerr << "SelectSubsetVertices: mnp=" << mnp << " mnp_red=" << mnp_red << "\n";
  std::vector<int> MapNode(mnp, -1);
  std::vector<int> MapNodeRev(mnp_red,-1);
  int idx=0;
  for (int iP=0; iP<mnp; iP++) {
    if (ListStatus[iP] == 1) {
      MapNodeRev[idx] = iP;
      MapNode[iP] = idx;
      idx++;
    }
  }
  int mne_red = 0;
  for (int ie=0; ie<mne; ie++) {
    int i0 = GrdArr.INE(ie,0);
    int i1 = GrdArr.INE(ie,1);
    int i2 = GrdArr.INE(ie,2);
    int stat0 = ListStatus[i0];
    int stat1 = ListStatus[i1];
    int stat2 = ListStatus[i2];
    if (stat0 == 1 && stat1 == 1 && stat2 == 1) {
      mne_red++;
    }
  }
  std::vector<int> MapEle(mne, -1);
  std::vector<int> MapEleRev(mne_red,-1);
  std::cerr << "SelectSubsetVertices: mne=" << mne << " mne_red=" << mne_red << "\n";
  MyMatrix<int> INEred(mne_red, 3);
  int pos = 0;
  for (int ie=0; ie<mne; ie++) {
    int i0 = GrdArr.INE(ie, 0);
    int i1 = GrdArr.INE(ie, 1);
    int i2 = GrdArr.INE(ie, 2);
    int j0 = MapNode[i0];
    int j1 = MapNode[i1];
    int j2 = MapNode[i2];
    if (j0 != -1 && j1 != -1 && j2 != -1) {
      INEred(pos, 0) = j0;
      INEred(pos, 1) = j1;
      INEred(pos, 2) = j2;
      MapEle[ie] = pos;
      MapEleRev[pos] = ie;
      pos++;
    }
  }
  MyMatrix<double> const& DEP = *GrdArr.GrdArrRho.DEP;
  MyMatrix<double> LONred(mnp_red,1), LATred(mnp_red,1), DEPred(mnp_red,1);
  for (int i=0; i<mnp_red; i++) {
    int iP = MapNodeRev[i];
    LONred(i,0) += GrdArr.GrdArrRho.LON(iP,0);
    LATred(i,0) += GrdArr.GrdArrRho.LAT(iP,0);
    DEPred(i,0) += DEP(iP,0);
  }
  GridArray GrdArrRet = GrdArr;
  GrdArrRet.INE = INEred;
  GrdArrRet.GrdArrRho.LON = LONred;
  GrdArrRet.GrdArrRho.LAT = LATred;
  GrdArrRet.GrdArrRho.DEP = DEPred;
  // It is impossible to apply a check on output, since the produced grid
  // may be incorrect and need further treatment before being usable.
  // So do not call CHECK_UnstructuredGrid
  return GrdArrRet;
}


GridArray KeepLargestConnectedComponent(GridArray const& GrdArr) {
  std::unordered_set<std::pair<int,int>> SetEdges;
  auto f_insert=[&](int const& u, int const& v) -> void {
    if (u < v)
      SetEdges.insert({u,v});
    else
      SetEdges.insert({v,u});
  };
  int mne = GrdArr.INE.rows();
  int mnp = GrdArr.GrdArrRho.LON.rows();
  for (int ie=0; ie<mne; ie++) {
    for (int i=0; i<3; i++) {
      int j = (i + 1) % 3;
      int ip = GrdArr.INE(ie, i);
      int jp = GrdArr.INE(ie, j);
      f_insert(ip, jp);
    }
  }
  int n_edge = SetEdges.size();
  std::cerr << "n_edge=" << n_edge << "\n";
  MyMatrix<size_t> ListEdge(n_edge,2);
  int pos = 0;
  int max_list_edges = 0;
  for (auto & eEdge : SetEdges) {
    ListEdge(pos,0) = eEdge.first;
    ListEdge(pos,1) = eEdge.second;
    if (eEdge.first > max_list_edges)
      max_list_edges = eEdge.first;
    if (eEdge.second > max_list_edges)
      max_list_edges = eEdge.second;
    pos++;
  }
  std::cerr << "max_list_edges=" << max_list_edges << " mnp=" << mnp << "\n";
  GraphListAdj GR(ListEdge, mnp);
  std::vector<size_t> ListStatus = ConnectedComponents_vector(GR);
  int nbConn = VectorMax(ListStatus) + 1;
  std::vector<int> ListConnSize(nbConn, 0);
  for (int i=0; i<mnp; i++) {
    int iConn = ListStatus[i];
    ListConnSize[iConn] += 1;
  }
  size_t iConnMax = std::numeric_limits<size_t>::max();;
  int maxConnSize = -1;
  for (int iConn=0; iConn<nbConn; iConn++) {
    int ConnSize = ListConnSize[iConn];
    if (ConnSize > maxConnSize) {
      maxConnSize = ConnSize;
      iConnMax = iConn;
    }
  }
  std::vector<int> ListStatus_B(mnp,0);
  for (int iP=0; iP<mnp; iP++) {
    if (ListStatus[iP] == iConnMax)
      ListStatus_B[iP] = 1;
  }
  return SelectSubsetVertices(GrdArr, ListStatus_B);
}

GridArray RemoveIsolatedPoints(GridArray const& GrdArr) {
  int mne = GrdArr.INE.rows();
  int mnp = GrdArr.GrdArrRho.LON.rows();
  std::vector<int> ListStatus(mnp);
  for (int ie=0; ie<mne; ie++) {
    for (int i=0; i<3; i++) {
      int ip = GrdArr.INE(ie, i);
      ListStatus[ip] = 1;
    }
  }
  return SelectSubsetVertices(GrdArr, ListStatus);
}



// clang-format off
#endif // SRC_OCEAN_TRIANGULATIONS_H_
// clang-format on
