#ifndef KERNEL_TRANSECT_INCLUDE
#define KERNEL_TRANSECT_INCLUDE


#include "NamelistExampleOcean.h"
#include "SphericalGeom.h"
#include "Interpolation.h"
#include "Model_interpolation.h"
#include "Model_grids.h"
#include "CommonFuncModel.h"
#include "Basic_plot.h"


struct PointOutTrans {
  double lon;
  double lat;
  double depth;
  std::string name;
};


std::vector<PointOutTrans> ReadStationCoordinate_File(std::string const& eFile)
{
  std::ifstream is(eFile);
  std::cerr << "We have is\n";
  std::string line;
  bool IsFirst=true;
  std::vector<PointOutTrans> ListPoint;
  while(!is.eof()) {
    std::getline(is, line);
    int siz=line.size();
    //    std::cerr << "line=" << line << " siz=" << siz << "\n";
    if (!IsFirst && siz > 0) {
      std::vector<std::string> LStr=STRING_Split(line, "\t");
      if (LStr.size() != 4) {
	std::cerr << "The size of LStr must be 4\n";
	std::cerr << " |LStr|=" << LStr.size() << "\n";
	throw TerminalException{1};
      }
      //      std::cerr << "  After the split\n";
      std::string name=LStr[0];
      double lon, lat, depth;
      //      std::cerr << "LStr[1]=" << LStr[1] << "\n";
      //      std::cerr << "LStr[2]=" << LStr[2] << "\n";
      //      std::cerr << "LStr[3]=" << LStr[3] << "\n";
      std::string strLon=STRING_Replace(LStr[3], ",", ".");
      //      std::cerr << "strLon=" << strLon << "\n";
      std::string strLat=STRING_Replace(LStr[2], ",", ".");
      //      std::cerr << "strLat=" << strLat << "\n";
      std::string strDep=STRING_Replace(LStr[1], ",", ".");
      //      std::cerr << "strDep=" << strDep << "\n";
      std::istringstream(strDep) >> depth;
      std::istringstream(strLat) >> lat;
      std::istringstream(strLon) >> lon;
      //      std::cerr << "After the istringstreams\n";
      ListPoint.push_back({lon, lat, depth, name});
    }
    IsFirst=false;
  }
  return ListPoint;
}




MyVector<double> ObtainDimensionVariable(std::vector<PairLL> const& ListPair)
{
  int siz=ListPair.size();
  MyVector<double> ListLon(siz), ListLat(siz);
  for (int i=0; i<siz; i++) {
    ListLon(i)=ListPair[i].eLon;
    ListLat(i)=ListPair[i].eLat;
  }
  double RangeLat=ListLat.maxCoeff() - ListLat.minCoeff();
  double RangeLon=ListLon.maxCoeff() - ListLon.minCoeff();
  if (RangeLat > RangeLon) {
    return ListLat;
  }
  else {
    return ListLon;
  }
}


std::vector<bool> DetermineBelonging_ListXY(GridArray const& GrdArr, MyMatrix<double> const& ListXY)
{
  int nbPoint=ListXY.cols();
  std::vector<SingleRecInterp> LSingle=General_FindInterpolationWeight(GrdArr, ListXY);
  std::vector<bool> ListStatus(nbPoint,true);
  for (int i=0; i<nbPoint; i++) {
    if (!LSingle[i].status)
      ListStatus[i]=false;
    std::cerr << "DetermineBelonging i=" << i << " status=" << LSingle[i].status << "\n";
  }
  return ListStatus;
}


SingleArrayInterpolation ComputeArrayInterpolation_ListXY(GridArray const& GrdArr, MyMatrix<double> const& ListXY)
{
  int nbPoint=ListXY.cols();
  std::vector<SingleRecInterp> LSingle=General_FindInterpolationWeight(GrdArr, ListXY);
  int eta_in=GrdArr.GrdArrRho.LON.rows();
  int xi_in =GrdArr.GrdArrRho.LON.cols();
  int nbError=0;
  for (int i=0; i<nbPoint; i++) {
    if (!LSingle[i].status) {
      std::cerr << "Point i=" << i << " / " << nbPoint << "\n";
      std::cerr << "X=" << ListXY(0,i) << " Y=" << ListXY(1,i) << "\n";
      std::cerr << "is found outside the domain\n";
      nbError++;
    }
    std::cerr << "ComputeArrayInterpolation i=" << i << " status=" << LSingle[i].status << "\n";
  }
  if (nbError > 0) {
    std::cerr << "Some points are outside\n";
    PrintGridArray(std::cerr, GrdArr);
    throw TerminalException{1};
  }
  int eta_out=nbPoint;
  int xi_out=1;
  std::vector<int> LEta(eta_out), LXi(eta_out);
  for (int i=0; i<eta_out; i++) {
    LEta[i]=i;
    LXi[i]=0;
  }
  ARVDtyp ARVDin=GetTrivialARrayVerticalDescription();
  return ConvertToArrayInt(eta_out, xi_out, eta_in, xi_in, LEta, LXi, LSingle, GrdArr, ARVDin);
}





TransectInformation GetTransectInformation(std::vector<GridArray> const& ListGrdArr,
					   double const& eLonStart, double const& eLatStart, 
					   double const& eLonEnd, double const& eLatEnd,
					   double const& eResolKM)
{
  int nbGrid=ListGrdArr.size();
  bool IsSpherical=ListGrdArr[0].IsSpherical;
  for (int iGrid=0; iGrid<nbGrid; iGrid++)
    if (ListGrdArr[iGrid].IsSpherical != IsSpherical) {
      std::cerr << "We should have all the grids spherical or not\n";
      std::cerr << "and that seems not to be the case\n";
      throw TerminalException{1};
    }
  double eDistKM;
  if (IsSpherical) {
    eDistKM=GeodesicDistanceKM(eLonStart, eLatStart, eLonEnd, eLatEnd);
  }
  else {
    double diffLON=eLonStart - eLonEnd;
    double diffLAT=eLatStart - eLatEnd;
    eDistKM=sqrt(diffLON*diffLON + diffLAT*diffLAT);
  }
  int NbSubdi=int(round(eDistKM/eResolKM));
  std::cerr << "eDistKM=" << eDistKM << "\n";
  std::cerr << "eResolKM=" << eResolKM << "\n";
  std::cerr << "NbSubdi=" << NbSubdi << "\n";
  double deltaLON=(eLonEnd - eLonStart)/double(NbSubdi);
  double deltaLAT=(eLatEnd - eLatStart)/double(NbSubdi);
  std::vector<PairLL> ListPairLL(NbSubdi+1);
  MyMatrix<double> ListXY(2,NbSubdi+1);
  for (int i=0; i<=NbSubdi; i++) {
    double eLon=eLonStart + double(i)*deltaLON;
    double eLat=eLatStart + double(i)*deltaLAT;
    PairLL ePair{eLon, eLat};
    ListPairLL[i]=ePair;
    ListXY(0,i)=eLon;
    ListXY(1,i)=eLat;
  }
  TransectInformation eTransect;
  eTransect.ListPairLL=ListPairLL;
  eTransect.ListDimVar=ObtainDimensionVariable(ListPairLL);
  for (int iGrid=0; iGrid<nbGrid; iGrid++)
    eTransect.ListRec.push_back(ComputeArrayInterpolation_ListXY(ListGrdArr[iGrid], ListXY));
  return eTransect;
}


TransectInformation_3D GetTransectInformation_3D(TransectInformation const& eTrans, 
						 GridArray const& GrdArr,
						 Eigen::Tensor<double,3> const VertCoord, 
						 double const& VertResolM)
{
  TransectInformation_3D eTrans3;
  eTrans3.ListPairLL=eTrans.ListPairLL;
  eTrans3.ListDimVar=eTrans.ListDimVar;
  eTrans3.eRecInterp=eTrans.ListRec[0];
  //
  int nbPoint=eTrans.ListPairLL.size();
  auto LDim=VertCoord.dimensions();
  int Nvert=LDim[0];
  Eigen::Tensor<double,3> TheRes=SingleInterpolationOfField_3D(eTrans.ListRec[0], VertCoord);
  MyMatrix<double> VertCoordS=DimensionExtraction(TheRes, 2, 0);
  MyMatrix<double> DEPinterp=SingleInterpolationOfField_2D(eTrans.ListRec[0], GrdArr.GrdArrRho.DEP);
  double maxDep=DEPinterp.maxCoeff();
  int NbVert=int(round(maxDep / VertResolM)) + 1;
  double DeltaZ=maxDep / double(NbVert);
  MyVector<double> ListVertPos(NbVert+1);
  for (int i=0; i<=NbVert; i++) {
    double eVertPos = -maxDep + i*DeltaZ;
    ListVertPos(i)=eVertPos;
  }
  int Ntotal=NbVert+1;
  MyMatrix<double> matCoeffM(Ntotal,nbPoint);
  MyMatrix<double> matCoeffP(Ntotal,nbPoint);
  MyMatrix<int> matIdxM(Ntotal,nbPoint);
  MyMatrix<int> matIdxP(Ntotal,nbPoint);
  MyMatrix<int> MSK(Ntotal,nbPoint);
  std::vector<int> FirstWetIndex(Ntotal,-1);
  double eps=0.00001;
  for (int iPt=0; iPt<nbPoint; iPt++)
    for (int i=0; i<=NbVert; i++) {
      double eVert=ListVertPos(i);
      int idxM=-1;
      int idxP=-1;
      double coefP=0;
      double coefM=0;
      int eMSK=0;
      if (eVert >= -DEPinterp(iPt,0)) {
	eMSK=1;
	if (eVert < VertCoordS(0,iPt) + eps) {
	  idxM=0;
	  idxP=0;
	  coefP=1;
	  coefM=0;
	}
	else {
	  if (eVert > VertCoordS(Nvert-1,iPt) - eps) {
	    idxM=Nvert-1;
	    idxP=Nvert-1;
	    coefP=1;
	    coefM=0;
	  }
	  else {
	    bool IsMatch=false;
	    for (int iVert=0; iVert<Nvert-1; iVert++) {
	      double dep1=VertCoordS(iVert  ,iPt);
	      double dep2=VertCoordS(iVert+1,iPt);
	      if (eVert >= dep1 && eVert <= dep2) {
		IsMatch=true;
		double alpha1=(dep2 - eVert)/(dep2 - dep1);
		double alpha2=(eVert - dep1)/(dep2 - dep1);
		idxM=iVert;
		idxP=iVert+1;
		coefM=alpha1;
		coefP=alpha2;
	      }
	    }
	    if (!IsMatch) {
	      std::cerr << "Failed to find matching depth\n";
	      throw TerminalException{1};
	    }
	  }
	}
      }
      matCoeffM(i, iPt)=coefM;
      matCoeffP(i, iPt)=coefP;
      matIdxM(i, iPt)=idxM;
      matIdxP(i, iPt)=idxP;
      MSK(i,iPt)=eMSK;
      if (eMSK == 1 && FirstWetIndex[iPt] == -1)
	FirstWetIndex[iPt]=i;
    }
  for (int iPt=0; iPt<nbPoint; iPt++) {
    std::cerr << "iPt=" << iPt << " FirstWetIndex=" << FirstWetIndex[iPt] << "\n";
  }
  eTrans3.matCoeffM=matCoeffM;
  eTrans3.matCoeffP=matCoeffP;
  eTrans3.matIdxM=matIdxM;
  eTrans3.matIdxP=matIdxP;
  eTrans3.MSK=MSK;
  eTrans3.ListVertPos=ListVertPos;
  eTrans3.Ntotal=Ntotal;
  eTrans3.FirstWetIndex=FirstWetIndex;
  //
  double deltaLON=eTrans.ListPairLL[nbPoint-1].eLon - eTrans.ListPairLL[0].eLon;
  double deltaLAT=eTrans.ListPairLL[nbPoint-1].eLat - eTrans.ListPairLL[0].eLat;
  double norm=sqrt(pow(deltaLON, 2) + pow(deltaLAT, 2));
  double normLON=deltaLON / norm;
  double normLAT=deltaLAT / norm;
  eTrans3.normU = -normLAT;
  eTrans3.normV = normLON;
  //
  return eTrans3;
}




GridArray GetGridArrayFromTransect3(TransectInformation_3D const& eTrans3)
{
  GridArray GrdArr;
  GrdArr.IsFE=0;
  GrdArr.IsSpherical=false;
  GrdArr.GrdArrRho.MSK=eTrans3.MSK;
  int nbPoint=eTrans3.ListPairLL.size();
  int Ntotal=eTrans3.Ntotal;
  MyMatrix<double> LAT(Ntotal, nbPoint);
  MyMatrix<double> LON(Ntotal, nbPoint);
  for (int iPt=0; iPt<nbPoint; iPt++)
    for (int iVert=0; iVert<Ntotal; iVert++) {
      LAT(iVert,iPt)=eTrans3.ListVertPos(iVert);
      LON(iVert,iPt)=eTrans3.ListDimVar(iPt);
    }
  GrdArr.GrdArrRho.LAT=LAT;
  GrdArr.GrdArrRho.LON=LON;
  return GrdArr;
}




MyMatrix<double> TransectInterpolation_3D(TransectInformation_3D const& eTrans3, Eigen::Tensor<double,3> const& F)
{
  //  std::cerr << "max(F)=" << maxCoeff(F) << "\n";
  Eigen::Tensor<double,3> F2=SingleInterpolationOfField_3D(eTrans3.eRecInterp, F);
  //  std::cerr << "max(F2)=" << maxCoeff(F2) << "\n";
  
  MyMatrix<double> F3=DimensionExtraction(F2, 2, 0);
  //  std::cerr << "max(F3)=" << F3.maxCoeff() << "\n";
  int Ntotal=eTrans3.Ntotal;
  int nbPoint=eTrans3.ListPairLL.size();
  MyMatrix<double> F4(Ntotal, nbPoint);
  for (int iVert=0; iVert<Ntotal; iVert++)
    for (int iPt=0; iPt<nbPoint; iPt++) {
      double sum=0;
      int idxM=eTrans3.matIdxM(iVert, iPt);
      int idxP=eTrans3.matIdxP(iVert, iPt);
      if (idxM >= 0) {
	double coefM=eTrans3.matCoeffM(iVert,iPt);
	double coefP=eTrans3.matCoeffP(iVert,iPt);
	sum = coefM * F3(idxM,iPt) + coefP * F3(idxP,iPt);
      }
      F4(iVert,iPt)=sum;
    }
  //  std::cerr << "max(F4)=" << F4.maxCoeff() << "\n";
  return F4;
}











std::vector<TransectInformation> RetrieveListTransect(SingleBlock const& eBlPLOT, std::vector<GridArray> const& ListGrdArr)
{
  std::vector<double> ListLonStart=eBlPLOT.ListListDoubleValues.at("TransectLonStart");
  std::vector<double> ListLatStart=eBlPLOT.ListListDoubleValues.at("TransectLatStart");
  std::vector<double> ListLonEnd=eBlPLOT.ListListDoubleValues.at("TransectLonEnd");
  std::vector<double> ListLatEnd=eBlPLOT.ListListDoubleValues.at("TransectLatEnd");
  std::vector<double> ListResolKM=eBlPLOT.ListListDoubleValues.at("SpatialResolutionTransectKM");
  int nbTrans=ListResolKM.size();
  size_t nbTrans_t=ListResolKM.size();
  if (nbTrans_t != ListLonStart.size() || nbTrans_t != ListLatStart.size() || nbTrans_t != ListLonEnd.size() || nbTrans_t != ListLonEnd.size()) {
    std::cerr << "Error in the transect information input\n";
    std::cerr << "The number should all be the same for:\n";
    std::cerr << "|TransectLonStart|=" << ListLonStart.size() << "\n";
    std::cerr << "|TransectLatStart|=" << ListLatStart.size() << "\n";
    std::cerr << "  |TransectLonEnd|=" << ListLonEnd.size() << "\n";
    std::cerr << "  |TransectLatEnd|=" << ListLatEnd.size() << "\n";
    std::cerr << "  |SpatialResolutionTransectKM|=" << ListResolKM.size() << "\n";
    std::cerr << "TransectLonStart, TransectLatStart, TransectLonEnd, TransectLatEnd and SpatialResolutionTransectKM\n";
    throw TerminalException{1};
  }
  std::vector<TransectInformation> ListTransect(nbTrans);
  for (int iTrans=0; iTrans<nbTrans; iTrans++)
    ListTransect[iTrans]=GetTransectInformation(ListGrdArr,
						ListLonStart[iTrans], ListLatStart[iTrans],
						ListLonEnd[iTrans], ListLatEnd[iTrans],
						ListResolKM[iTrans]);
  return ListTransect;
}










#endif
