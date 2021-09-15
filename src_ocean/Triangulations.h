#ifndef TRIANGULATION_FCT_DEFINE
#define TRIANGULATION_FCT_DEFINE

#include "Basic_Ocean_types.h"
#include "GRAPH_GraphicalBasic.h"
#include "SphericalGeom.h"


GraphSparseImmutable GetUnstructuredVertexAdjInfo(MyMatrix<int> const& INE, size_t nbNode)
{
  size_t nbEle=INE.rows();
  std::vector<size_t> ListNbEnt(nbNode,0);
  for (size_t iEle=0; iEle<nbEle; iEle++)
    for (size_t i=0; i<3; i++) {
      size_t eVert=INE(iEle,i);
      ListNbEnt[eVert]+=2;
    }
  std::cerr << "nbNode=" << nbNode << "\n";
  size_t TotalSum_unrefined=6*nbEle;
  size_t miss_val=std::numeric_limits<size_t>::max();
  std::vector<size_t> ListStart_unrefined(nbNode+1);
  ListStart_unrefined[0] = 0;
  for (size_t iNode=0; iNode<nbNode; iNode++)
    ListStart_unrefined[iNode+1]=ListStart_unrefined[iNode] + ListNbEnt[iNode];
  std::vector<size_t> ListListAdj_unrefined(TotalSum_unrefined,miss_val);
  std::vector<size_t> ListIndexPos(nbNode,0);
  auto fInsert=[&](size_t const& eVert, size_t const& eVertAdj) -> void {
    size_t eStart=ListStart_unrefined[eVert];
    size_t eEnd=ListStart_unrefined[eVert] + ListIndexPos[eVert];
    if (ListListAdj_unrefined[eEnd] != miss_val) {
      std::cerr << "Logical error in the code\n";
      throw TerminalException{1};
    }
    for (size_t i=eStart; i<eEnd; i++)
      if (ListListAdj_unrefined[i] == eVertAdj)
	return;
    ListIndexPos[eVert]++;
    ListListAdj_unrefined[eEnd]=eVertAdj;
  };
  for (size_t iEle=0; iEle<nbEle; iEle++)
    for (int i=0; i<3; i++) {
      int iNext=NextIdx(3,i);
      int iPrev=PrevIdx(3,i);
      size_t eVert=INE(iEle,i);
      size_t eVertP=INE(iEle,iPrev);
      size_t eVertN=INE(iEle,iNext);
      fInsert(eVert, eVertP);
      fInsert(eVert, eVertN);
    }
  std::vector<size_t> ListStart(nbNode+1,0);
  for (size_t iNode=0; iNode<nbNode; iNode++)
    ListStart[iNode+1]=ListStart[iNode] + ListIndexPos[iNode];
  size_t TotalSum=ListStart[nbNode];
  std::vector<size_t> ListListAdj(TotalSum);
  for (size_t iNode=0; iNode<nbNode; iNode++) {
    size_t eStart=ListStart[iNode];
    size_t eStart_unrefined=ListStart_unrefined[iNode];
    size_t siz=ListIndexPos[iNode];
    for (size_t i=0; i<siz; i++)
      ListListAdj[eStart + i]=ListListAdj_unrefined[eStart_unrefined+i];
  }
  return GraphSparseImmutable(nbNode, ListStart, ListListAdj);
}

// Computation of the boundary status
// 0: should not happen in return
// 1: should be inside node
// -1: should be boundary node
MyVector<int> GetBoundaryStatus(MyMatrix<int> const& INE, int nbNode)
{
  MyVector<int> Status(nbNode);
  for (int i=0; i<nbNode; i++)
    Status(i)=0;
  std::vector<int> PrevVert(nbNode), NextVert(nbNode), Collected(nbNode);
  int mne=INE.rows();
  for (int ie=0; ie<mne; ie++) {
    for (int i=0; i<3; i++) {
      int iPrev=PrevIdx(3, i);
      int iNext=NextIdx(3, i);
      int ip=INE(ie,i);
      int ipnext=INE(ie,iNext);
      int ipprev=INE(ie,iPrev);
      if (Status(ip) == 0) {
	Status(ip)=1;
	PrevVert[ip]=ipprev;
	NextVert[ip]=ipnext;
      }
    }
  }
  for (int i=0; i<nbNode; i++)
    Status(i)=0;
  while(true) {
    for (int i=0; i<nbNode; i++)
      Collected[i]=0;
    for (int ie=0; ie<mne; ie++) {
      for (int i=0; i<3; i++) {
	int iPrev=PrevIdx(3, i);
	int iNext=NextIdx(3, i);
	int ip=INE(ie,i);
	int ipnext=INE(ie,iNext);
	int ipprev=INE(ie,iPrev);
	if (Status(ip) == 0) {
	  int zNext=NextVert[ip];
	  if (zNext == ipprev) {
	    Collected[ip]=1;
	    NextVert[ip]=ipnext;
	    if (NextVert[ip] == PrevVert[ip])
	      Status(ip)=1;
	  }
	}
      }
    }
    int IsFinished=1;
    for (int i=0; i<nbNode; i++) {
      if (Collected[i] == 0 && Status(i) == 0)
	Status(i)=-1;
      if (Status(i) == 0)
	IsFinished=0;
    }
    if (IsFinished == 1)
      break;
  }
  return Status;
}


std::vector<int> GetListNeighbor(MyMatrix<int> const& INE, int nbNode)
{
  std::vector<int> Status(nbNode,0), Neighbor(nbNode,-1);
  std::vector<int> PrevVert(nbNode), NextVert(nbNode), Collected(nbNode);
  int mne=INE.rows();
  int sizMesh=INE.cols();
  for (int ie=0; ie<mne; ie++) {
    for (int i=0; i<sizMesh; i++) {
      int iPrev=PrevIdx(sizMesh, i);
      int iNext=NextIdx(sizMesh, i);
      int ip=INE(ie,i);
      int ipnext=INE(ie,iNext);
      int ipprev=INE(ie,iPrev);
      if (Status[ip] == 0) {
	Status[ip]=1;
	PrevVert[ip]=ipprev;
	NextVert[ip]=ipnext;
      }
    }
  }
  for (int i=0; i<nbNode; i++)
    Status[i]=0;
  while(true) {
    for (int i=0; i<nbNode; i++)
      Collected[i]=0;
    for (int ie=0; ie<mne; ie++) {
      for (int i=0; i<sizMesh; i++) {
	int iPrev=PrevIdx(sizMesh, i);
	int iNext=NextIdx(sizMesh, i);
	int ip=INE(ie,i);
	int ipnext=INE(ie,iNext);
	int ipprev=INE(ie,iPrev);
	if (Status[ip] == 0) {
	  int zNext=NextVert[ip];
	  if (zNext == ipprev) {
	    Collected[ip]=1;
	    NextVert[ip]=ipnext;
	    if (NextVert[ip] == PrevVert[ip])
	      Status[ip]=1;
	  }
	}
      }
    }
    bool IsFinished=true;
    for (int i=0; i<nbNode; i++) {
      if (Collected[i] == 0 && Status[i] == 0) {
	Status[i]=-1;
	Neighbor[i]=NextVert[i];
      }
      if (Status[i] == 0)
	IsFinished=false;
    }
    if (IsFinished)
      break;
  }
  return Neighbor;
}







// There is a more advanced version of this code that handle
// pathological cases in WWM code in wwm_netcdf.F90 in the section
// See also interpol.F90 in polymesh code.
//
// Here we use the simpler code.
std::vector<std::vector<int>> GetListBoundaryCycles(MyMatrix<int> const& INE, int nbNode)
{
  std::vector<int> Neighbor=GetListNeighbor(INE, nbNode);
  int nbBound=0;
  for (int iVert=0; iVert<nbNode; iVert++)
    if (Neighbor[iVert] >= 0)
      nbBound++;
  std::vector<int> ListMap(nbNode,-1);
  std::vector<int> ListMapRev(nbBound);
  std::cerr << "nbNode=" << nbNode << " nbBound=" << nbBound << "\n";
  int idx=0;
  for (int iVert=0; iVert<nbNode; iVert++)
    if (Neighbor[iVert] >= 0) {
      ListMap[iVert]=idx;
      ListMapRev[idx]=iVert;
      idx++;
    }
  std::vector<int> NeighborRed(nbBound);
  for (int iB=0; iB<nbBound; iB++) {
    int iVert=ListMapRev[iB];
    int iVertImg=Neighbor[iVert];
    int iBimg=ListMap[iVertImg];
    NeighborRed[iB]=iBimg;
  }
  std::vector<int> Status(nbBound,-1);
  auto GetOneUnset=[&]() -> int {
    for (int iB=0; iB<nbBound; iB++)
      if (Status[iB] == -1)
	return iB;
    return -1;
  };
  std::vector<std::vector<int>> eListList;
  while(true) {
    int eFirst=GetOneUnset();
    if (eFirst == -1)
      break;
    std::vector<int> eList;
    int idxB=eFirst;
    while(true) {
      int eVert=ListMapRev[idxB];
      eList.push_back(eVert);
      if (int(eList.size()) > nbNode) {
	std::cerr << "eList is clearly too large\n";
	std::cerr << "Logical error\n";
	throw TerminalException{1};
      }
      Status[idxB]=0;
      idxB=NeighborRed[idxB];
      if (idxB == eFirst)
	break;
    }
    eListList.push_back(eList);
  }
  return eListList;
}



MyMatrix<int> GetEdgeSet(MyMatrix<int> const& INE, int nbNode)
{
  std::vector<int> ListDegree(nbNode, 0);
  int mne=INE.rows();
  for (int ie=0; ie<mne; ie++)
    for (int i=0; i<3; i++) {
      int ip=INE(ie,i);
      ListDegree[ip]+=2;
    }
  int MaxDeg=0;
  for (int iNode=0; iNode<nbNode; iNode++) {
    int eDeg=ListDegree[iNode];
    if (eDeg > MaxDeg)
      MaxDeg=eDeg;
  }
  for (int iNode=0; iNode<nbNode; iNode++)
    ListDegree[iNode]=0;
  MyMatrix<int> ListAdjacency(nbNode, MaxDeg);
  for (int ie=0; ie<mne; ie++) {
    int i1=INE(ie,0);
    int i2=INE(ie,1);
    int i3=INE(ie,2);
    int eDeg1=ListDegree[i1];
    int eDeg2=ListDegree[i2];
    int eDeg3=ListDegree[i3];
    ListAdjacency(i1, eDeg1  )=i2;
    ListAdjacency(i1, eDeg1+1)=i3;
    ListAdjacency(i2, eDeg2  )=i1;
    ListAdjacency(i2, eDeg2+1)=i3;
    ListAdjacency(i3, eDeg3  )=i1;
    ListAdjacency(i3, eDeg3+1)=i2;
    ListDegree[i1] = eDeg1 + 2;
    ListDegree[i2] = eDeg2 + 2;
    ListDegree[i3] = eDeg3 + 2;
  }
  int nbEdge=0;
  for (int iNode=0; iNode<nbNode; iNode++) {
    std::set<int> eSet;
    int eDeg=ListDegree[iNode];
    for (int iAdj=0; iAdj<eDeg; iAdj++) {
      int eAdj=ListAdjacency(iNode, iAdj);
      if (eAdj > iNode)
	eSet.insert(eAdj);
    }
    nbEdge += eSet.size();
  }
  MyMatrix<int> ListEdges(nbEdge,2);
  int iEdge=0;
  for (int iNode=0; iNode<nbNode; iNode++) {
    std::set<int> eSet;
    int eDeg=ListDegree[iNode];
    for (int iAdj=0; iAdj<eDeg; iAdj++) {
      int eAdj=ListAdjacency(iNode, iAdj);
      if (eAdj > iNode)
	eSet.insert(eAdj);
    }
    for (auto eAdj : eSet) {
      ListEdges(iEdge,0)=iNode;
      ListEdges(iEdge,1)=eAdj;
      iEdge++;
    }
  }
  return ListEdges;
}


MyMatrix<int> GetFaceEdgeConnectivity(int const& nbNode, MyMatrix<int> const& LEdge, MyMatrix<int> const& INE)
{
  int nbEdge=LEdge.rows();
  std::vector<int> eVect(2*nbEdge,0);
  std::vector<int> eVectNB(nbNode,0);
  for (int iEdge=0; iEdge<nbEdge; iEdge++) {
    for (int i=0; i<2; i++) {
      int ePt = LEdge(iEdge,i);
      eVectNB[ePt]++;
    }
  }
  std::vector<int> ListShift(nbNode,0);
  for (int iNode=1; iNode<nbNode; iNode++)
    ListShift[iNode] = ListShift[iNode-1] + eVectNB[iNode-1];
  int NbIncidence = ListShift[nbNode-1] + eVectNB[nbNode-1];
  MyMatrix<int> PairInfo(NbIncidence,2);
  std::vector<int> eVectNB_att(nbNode,0);
  for (int iEdge=0; iEdge<nbEdge; iEdge++) {
    for (int i=0; i<2; i++) {
      int j = 1 - i;
      int ePt1 = LEdge(iEdge,i);
      int ePt2 = LEdge(iEdge,j);
      int eShift = ListShift[ePt1];
      int pos = eVectNB_att[ePt1];
      PairInfo(eShift+pos, 0) = ePt2;
      PairInfo(eShift+pos, 1) = iEdge;
      eVectNB_att[ePt1] = pos+1;
    }
  }
  auto GetIedge=[&](int const& ePt1, int const& ePt2) -> int {
    int eShift = ListShift[ePt1];
    int siz    = eVectNB[ePt1];
    for (int i=0; i<siz; i++) {
      if (PairInfo(eShift + i, 0) == ePt2)
	return PairInfo(eShift + i, 1);
    }
    return -1; // should not happen at all
  };
  MyMatrix<int> FaceEdgeConn(nbEdge,2);
  std::vector<int> nbEdgeAtt(nbEdge,0);
  int nbFace=INE.rows();
  for (int iFace=0; iFace<nbFace; iFace++) {
    for (int i=0; i<3; i++) {
      int j = (i+1) % 3;
      int ePt1=INE(iFace,i);
      int ePt2=INE(iFace,j);
      int iEdge=GetIedge(ePt1, ePt2);
      int pos=nbEdgeAtt[iEdge];
      FaceEdgeConn(iEdge,pos) = iFace;
      nbEdgeAtt[iEdge] = pos + 1;
    }
  }
  for (int iEdge=0; iEdge<nbEdge; iEdge++) {
    for (int pos=nbEdgeAtt[iEdge]+1; pos<2; pos++)
      FaceEdgeConn(iEdge,pos) = -999;
  }
  return FaceEdgeConn;
}





std::vector<int> GetUnstructuredTriangleAdjInfo_vectint_V1(MyMatrix<int> const& INE, int nbNode)
{
  int nbEle=INE.rows();
  std::vector<int> ListDegree(nbNode, 0);
  int mne=INE.rows();
  for (int ie=0; ie<mne; ie++)
    for (int i=0; i<3; i++) {
      int ip=INE(ie,i);
      ListDegree[ip]+=2;
    }
  int MaxDeg=0;
  for (int iNode=0; iNode<nbNode; iNode++) {
    int eDeg=ListDegree[iNode];
    if (eDeg > MaxDeg)
      MaxDeg=eDeg;
  }
  for (int iNode=0; iNode<nbNode; iNode++)
    ListDegree[iNode]=0;
  MyMatrix<int> ListAdjacency(nbNode, MaxDeg);
  for (int ie=0; ie<mne; ie++) {
    for (int i=0; i<3; i++) {
      int iNext=NextIdx(3,i);
      int iPrev=PrevIdx(3,i);
      int i1=INE(ie,i);
      int i2=INE(ie,iNext);
      int i3=INE(ie,iPrev);
      int eDeg=ListDegree[i1];
      ListAdjacency(i1, eDeg  )=i2;
      ListAdjacency(i1, eDeg+1)=i3;
      ListDegree[i1] = eDeg + 2;
    }
  }
  int nbEdge=0;
  for (int iNode=0; iNode<nbNode; iNode++) {
    std::set<int> eSet;
    int eDeg=ListDegree[iNode];
    for (int iAdj=0; iAdj<eDeg; iAdj++) {
      int eAdj=ListAdjacency(iNode, iAdj);
      if (eAdj > iNode)
	eSet.insert(eAdj);
    }
    nbEdge += eSet.size();
  }
  MyMatrix<int> ListEdges(nbEdge,2);
  int posEdge=0;
  std::vector<int> IndexStart(nbNode);
  std::vector<int> IndexEnd  (nbNode);
  for (int iNode=0; iNode<nbNode; iNode++) {
    std::set<int> eSet;
    int eDeg=ListDegree[iNode];
    for (int iAdj=0; iAdj<eDeg; iAdj++) {
      int eAdj=ListAdjacency(iNode, iAdj);
      if (eAdj > iNode)
	eSet.insert(eAdj);
    }
    IndexStart[iNode]=posEdge;
    for (auto eAdj : eSet) {
      ListEdges(posEdge,0)=iNode;
      ListEdges(posEdge,1)=eAdj;
      posEdge++;
    }
    IndexEnd  [iNode]=posEdge;
  }
  auto GetIEdge=[&](int const& eVert1, int const& eVert2) -> int {
    int idxStart=IndexStart[eVert1];
    int idxEnd=IndexEnd[eVert1];
    for (int iEdge=idxStart; iEdge<idxEnd; iEdge++) {
      if (ListEdges(iEdge,0) != eVert1) {
	std::cerr << "Clear inconsistency in code\n";
	throw TerminalException{1};
      }
      if (ListEdges(iEdge,1) == eVert2)
	return iEdge;
    }
    std::cerr << "Failed to find the correct indexes iEdgeIter\n";
    throw TerminalException{1};
  };
  MyMatrix<int> LEdge=GetEdgeSet(INE, nbNode);
  std::vector<int> NumberMatch(nbEdge, 0);
  MyMatrix<int> IncidenceTrigEdge(nbEdge,2);
  for (int iEle=0; iEle<nbEle; iEle++) {
    for (int i=0; i<3; i++) {
      int iNext=NextIdx(3,i);
      int eVert=INE(iEle, i);
      int eNext=INE(iEle, iNext);
      int eVert1, eVert2;
      if (eVert < eNext) {
	eVert1=eVert;
	eVert2=eNext;
      }
      else {
	eVert1=eNext;
	eVert2=eVert;
      }
      int iEdge=GetIEdge(eVert1, eVert2);
      int pos=NumberMatch[iEdge];
      IncidenceTrigEdge(iEdge,pos)=iEle;
      NumberMatch[iEdge]=pos+1;
    }
  }
  std::vector<int> ListAdj(3*nbEle,-1);
  std::vector<int> DegreeTriangle(nbEle,0);
  for (int iEdge=0; iEdge<nbEdge; iEdge++)
    if (NumberMatch[iEdge] == 2) {
      int iEle1=IncidenceTrigEdge(iEdge,0);
      int iEle2=IncidenceTrigEdge(iEdge,1);
      int pos1=DegreeTriangle[iEle1];
      int pos2=DegreeTriangle[iEle2];
      ListAdj[3*iEle1 + pos1] = iEle2;
      ListAdj[3*iEle2 + pos2] = iEle1;
      DegreeTriangle[iEle1]=pos1+1;
      DegreeTriangle[iEle2]=pos2+1;
    }
  return ListAdj;
}




std::vector<int> GetUnstructuredTriangleAdjInfo_vectint(MyMatrix<int> const& INE, int nbNode)
{
  std::vector<int> ListDegree(nbNode, 0);
  int mne=INE.rows();
  for (int ie=0; ie<mne; ie++)
    for (int i=0; i<3; i++) {
      int ip=INE(ie,i);
      ListDegree[ip]++;
    }
  int MaxDeg=0;
  for (int iNode=0; iNode<nbNode; iNode++) {
    int eDeg=ListDegree[iNode];
    if (eDeg > MaxDeg)
      MaxDeg=eDeg;
  }
  for (int iNode=0; iNode<nbNode; iNode++)
    ListDegree[iNode]=0;
  MyMatrix<int> IncidenceTrigVert(nbNode, MaxDeg);
  for (int ie=0; ie<mne; ie++) {
    for (int i=0; i<3; i++) {
      int i1=INE(ie,i);
      int eDeg=ListDegree[i1];
      IncidenceTrigVert(i1, eDeg) = ie;
      ListDegree[i1] = eDeg + 1;
    }
  }
  std::vector<int> ListAdj(3*mne);
  auto HasSpecifiedVertex=[&](int const& iTrig, int const & eVert) -> bool {
    for (int i=0; i<3; i++)
      if (INE(iTrig,i) == eVert)
	return true;
    return false;
  };
  for (int ie=0; ie<mne; ie++)
    for (int i=0; i<3; i++) {
      int iNext=NextIdx(3,i);
      int ip1=INE(ie,i);
      int ip2=INE(ie,iNext);
      int eAdj=-1;
      for (int iTrig=0; iTrig<ListDegree[ip1]; iTrig++) {
	int eTrig=IncidenceTrigVert(ip1, iTrig);
	if (HasSpecifiedVertex(eTrig, ip2) && eTrig != ie)
	  eAdj=eTrig;
      }
      ListAdj[3*ie + i] = eAdj;
    }
  return ListAdj;
}










void CHECK_UnstructuredGrid(GridArray const& GrdArr)
{
  int mnp=GrdArr.GrdArrRho.LON.rows();
  int mne=GrdArr.INE.rows();
  std::cerr << "mne=" << mne << "\n";
  int nbPlus=0;
  int nbMinus=0;
  for (int ie=0; ie<mne; ie++) {
    int i1=GrdArr.INE(ie,0);
    int i2=GrdArr.INE(ie,1);
    int i3=GrdArr.INE(ie,2);
    if (i1 == i2 || i1 == i3 || i2 == i3) {
      std::cerr << "For ie=" << ie << "\n";
      std::cerr << "We have i123=" << i1 << "," << i2 << "," << i3 << "\n";
    }
    double xi = GrdArr.GrdArrRho.LON(i1);
    double yi = GrdArr.GrdArrRho.LAT(i1);
    double xj = GrdArr.GrdArrRho.LON(i2);
    double yj = GrdArr.GrdArrRho.LAT(i2);
    double xk = GrdArr.GrdArrRho.LON(i3);
    double yk = GrdArr.GrdArrRho.LAT(i3);
    double area=xi*(yj-yk) + xj*(yk-yi) + xk*(yi-yj);
    if (area > 0)
      nbPlus++;
    if (area < 0)
      nbMinus++;
    for (int i=0; i<3; i++) {
      int IP=GrdArr.INE(ie,i);
      if (IP < 0 || IP >= mnp) {
	std::cerr << "Error in the unstructured grid\n";
	std::cerr << "mnp=" << mnp << "  mne=" << mne << "\n";
	std::cerr << "ie=" << ie << "\n";
	std::cerr << "INE=[" << GrdArr.INE(ie,0) << " , " << GrdArr.INE(ie,1) << " , " << GrdArr.INE(ie,2) << "]\n";
	throw TerminalException{1};
      }
    }
  }
  if (nbPlus > 0 && nbMinus > 0) {
    std::cerr << "The grid is incorrectly oriented\n";
    std::cerr << "mne=" << mne << " : nbPlus=" << nbPlus << "  nbMinus=" << nbMinus << "\n";
    throw TerminalException{1};
  }
  MyVector<int> Status=GetBoundaryStatus(GrdArr.INE, mnp);
  int nbStatusNormal=0;
  int nbStatusBound=0;
  for (int i=0; i<mnp; i++) {
    if (Status[i] == 1)
      nbStatusNormal++;
    if (Status[i] == -1)
      nbStatusBound++;
  }
  std::cerr << "nbStatusNormal=" << nbStatusNormal << " nbStatusBound=" << nbStatusBound << "\n";
  double MinDistKM=10000000;
  double MaxDistKM=0;
  double MinAreaSqrKM=10000000000000;
  double MaxAreaSqrKM=0;
  for (int ie=0; ie<mne; ie++) {
    for (int i1=0; i1<3; i1++) {
      int i2;
      if (i1 == 2) {
	i2=0;
      }
      else {
	i2=i1+1;
      }
      int pt1=GrdArr.INE(ie,i1);
      int pt2=GrdArr.INE(ie,i2);
      double eLon1 = GrdArr.GrdArrRho.LON(pt1);
      double eLat1 = GrdArr.GrdArrRho.LAT(pt1);
      double eLon2 = GrdArr.GrdArrRho.LON(pt2);
      double eLat2 = GrdArr.GrdArrRho.LAT(pt2);
      double eDistKM=GeodesicDistanceKM(eLon1, eLat1, eLon2, eLat2);
      if (eDistKM < MinDistKM)
	MinDistKM=eDistKM;
      if (eDistKM > MaxDistKM)
	MaxDistKM=eDistKM;
    }
    int pt1=GrdArr.INE(ie,0);
    int pt2=GrdArr.INE(ie,1);
    int pt3=GrdArr.INE(ie,2);
    double eLon1 = GrdArr.GrdArrRho.LON(pt1);
    double eLat1 = GrdArr.GrdArrRho.LAT(pt1);
    double eLon2 = GrdArr.GrdArrRho.LON(pt2);
    double eLat2 = GrdArr.GrdArrRho.LAT(pt2);
    double eLon3 = GrdArr.GrdArrRho.LON(pt3);
    double eLat3 = GrdArr.GrdArrRho.LAT(pt3);
    double areaKM=SphericalCoordinateAreaKM(eLon1, eLon2, eLon3, eLat1, eLat2, eLat3);
    if (areaKM < MinAreaSqrKM)
      MinAreaSqrKM=areaKM;
    if (areaKM > MaxAreaSqrKM)
      MaxAreaSqrKM=areaKM;
  }
  std::cerr << "dist km  : min=" << MinDistKM << " max=" << MaxDistKM << "\n";
  std::cerr << "area km2 : min=" << MinAreaSqrKM << " max=" << MaxAreaSqrKM << "\n";
  //
  std::vector<std::vector<int>> ListConn=GetListBoundaryCycles(GrdArr.INE, mnp);
  int nbConn=ListConn.size();
  int nbIsland=nbConn-1;
  std::cerr << "number of island resolved = " << nbIsland << "\n";
}

void CHECK_CombinatorialGrid(GridArray const& GrdArr)
{
  int mnp=GrdArr.GrdArrRho.LON.rows();
  int mne=GrdArr.INE.rows();
  std::vector<int> CCON(mnp,0);
  int POS_TRICK[3][2];
  POS_TRICK[0][0] = 1;
  POS_TRICK[1][0] = 2;
  POS_TRICK[2][0] = 0;
  POS_TRICK[0][1] = 2;
  POS_TRICK[1][1] = 0;
  POS_TRICK[2][1] = 1;
  for (int ie=0; ie<mne; ie++)
    for (int i=0; i<3; i++) {
      int ip=GrdArr.INE(ie,i);
      CCON[ip]++;
    }
  int MAXMNECON=0;
  for (int ip=0; ip<mnp; ip++) {
    int eCon=CCON[ip];
    if (eCon > MAXMNECON)
      MAXMNECON=eCon;
  }
  std::vector<int> CHILF(mnp,0);
  Eigen::Tensor<int,3> CELLVERTEX(mnp,MAXMNECON,2);
  for (int ie=0; ie<mne; ie++)
    for (int j=0; j<3; j++) {
      int i=GrdArr.INE(ie,j);
      CELLVERTEX(i, CHILF[i] ,0) = ie;
      CELLVERTEX(i, CHILF[i] ,1) = j;
      CHILF[i]++;
    }
  int COUNT_MAX=0;
  for (int ip=0; ip<mnp; ip++)
    COUNT_MAX += CCON[ip];
  MyMatrix<int> IE_CELL2 (mnp,MAXMNECON);
  MyMatrix<int> POS_CELL2(mnp,MAXMNECON);
  int j=0;
  for (int ip=0; ip<mnp; ip++)
    for (int i=0; i<CCON[ip]; i++) {
      IE_CELL2 (ip,i) = CELLVERTEX(ip,i,0);
      POS_CELL2(ip,i) = CELLVERTEX(ip,i,1);
      j++;
    }
  for (int ie=0; ie<mne; ie++)
    for (int i=0; i<3; i++) {
      int INEXT=POS_TRICK[i][0];
      int ip=GrdArr.INE(ie,i);
      int IP_NEXT=GrdArr.INE(ie, INEXT);
      int nbMatch=0;
      std::vector<int> Lmatch;
      for (int icon=0; icon<CCON[ip]; icon++) {
	int ie2=IE_CELL2(ip,icon);
	if (ie != ie2) {
	  int POS=POS_CELL2(ip,icon);
	  int POS_NEXT=POS_TRICK[POS][0];
	  int IP_ADJ_NEXT=GrdArr.INE(ie2, POS_NEXT);
	  if (IP_ADJ_NEXT == IP_NEXT) {
	    std::cerr << "Combinatorial orientability problem\n";
	    std::cerr << "IE=" << ie << " IE2=" << ie2 << "\n";
	    std::cerr << "IP=" << ip << " IP_NEXT=" << IP_NEXT << "\n";
	    throw TerminalException{1};
	  }
	  int POS_PREV=POS_TRICK[POS][1];
	  int IP_ADJ_PREV=GrdArr.INE(ie2, POS_PREV);
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
	for (int iMatch=0; iMatch<nbMatch; iMatch++) {
	  int iem=Lmatch[iMatch];
	  std::cerr << "  iMatch=" << iMatch << " ie=" << iem << "\n";
	  std::cerr << "     ine=[" << GrdArr.INE(iem,0) << "," << GrdArr.INE(iem,1) << "," << GrdArr.INE(iem,2) << "]\n";
	}
	throw TerminalException{1};
      }
    }
  std::cerr << "Now leaving the combinatorial check\n";
}




void CHECK_COORDINATE_ORIENTATION(GridArray const& GrdArr)
{
  int mne=GrdArr.INE.rows();
  int nbPlus=0;
  int nbMinus=0;
  for (int ie=0; ie<mne; ie++) {
    int i1=GrdArr.INE(ie, 0);
    int i2=GrdArr.INE(ie, 1);
    int i3=GrdArr.INE(ie, 2);
    double eLon1=GrdArr.GrdArrRho.LON(i1,0);
    double eLon2=GrdArr.GrdArrRho.LON(i2,0);
    double eLon3=GrdArr.GrdArrRho.LON(i3,0);
    double eLat1=GrdArr.GrdArrRho.LAT(i1,0);
    double eLat2=GrdArr.GrdArrRho.LAT(i2,0);
    double eLat3=GrdArr.GrdArrRho.LAT(i3,0);
    double deltaLON12=eLon2 - eLon1;
    double deltaLAT12=eLat2 - eLat1;
    double deltaLON13=eLon3 - eLon1;
    double deltaLAT13=eLat3 - eLat1;
    double eArea=deltaLON13 * deltaLAT12 - deltaLON12 * deltaLAT13;
    if (eArea > 0) {
      nbPlus++;
    }
    else {
      nbMinus++;
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


#endif
