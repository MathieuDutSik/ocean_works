#ifndef INCLUDE_SMOOTHING_BATHYMETRY
#define INCLUDE_SMOOTHING_BATHYMETRY

#include "Namelist.h"
#include "Model_grids.h"
#include "Triangulations.h"
#include "Model_interpolation.h"
#include "MAT_Matrix.h"
#include "GRAPH_GraphicalBasic.h"
#include "POLY_LinearProgramming_GLPK.h"


MyMatrix<double> GetRoughnessFactor(MyMatrix<double> const& TheBathy, GridArray const& GrdArr)
{
  std::pair<GraphSparseImmutable, std::vector<std::pair<int,int>>> eGLP = GetGraphSparseVertexAdjacency(GrdArr);
  int nb_point = eGLP.second.size();
  MyVector<double> VectFrac(nb_point,0);
  for (int iPoint=0; iPoint<nb_point; iPoint++) {
    std::pair<int,int> ePair = eGLP.second[iPoint];
    std::vector<int> ListAdj=eGLP.first.Adjacency(iPoint);
    double dep1 = GrdArr.GrdArrRho.DEP(ePair.first, ePair.second);
    double maxR = 0;
    for (auto & eADJ : ListAdj) {
      std::pair<int,int> fPair = eGLP.second[eADJ];
      double dep2 = GrdArr.GrdArrRho.DEP(fPair.first, fPair.second);
      double rFrac=T_abs(dep1 - dep2) / (dep1 + dep2);
      if (rFrac > maxR)
	maxR = rFrac;
    }
    VectFrac(iPoint) = maxR;
  }
  int eta_rho=TheBathy.rows();
  int xi_rho=TheBathy.cols();
  MyMatrix<double> Fret = ZeroMatrix<double>(eta_rho,xi_rho);
  for (int iPoint=0; iPoint<nb_point; iPoint++) {
    std::pair<int,int> ePair = eGLP.second[iPoint];
    Fret(ePair.first, ePair.second) = VectFrac(iPoint);
  }
  return Fret;
}




MyVector<int> GetBadPoints(GridArray const& GrdArr, double const& rx0max, int const& NeighborLevel)
{
  std::pair<GraphSparseImmutable, std::vector<std::pair<int,int>>> eGLP = GetGraphSparseVertexAdjacency(GrdArr);
  int nb_point = eGLP.second.size();
  MyMatrix<double> RoughMat=GetRoughnessFactor(GrdArr.GrdArrRho.DEP, GrdArr);
  MyVector<int> ListBadPoint=ZeroVector<int>(nb_point);
  double MaxVal=0;
  int nbBad=0;
  for (int iPoint=0; iPoint<nb_point; iPoint++) {
    std::pair<int,int> ePair = eGLP.second[iPoint];
    double eVal=RoughMat(ePair.first, ePair.second);
    if (eVal > MaxVal)
      MaxVal = eVal;
    if (eVal > rx0max) {
      ListBadPoint(iPoint)=1;
      nbBad++;
    }
  }
  std::cerr << "rx0 = " << MaxVal << " nbBad = " << nbBad << "\n";
  for (int iter=0; iter<NeighborLevel; iter++) {
    MyVector<int> NewListBadPoint = ListBadPoint;
    for (int iPoint=0; iPoint<nb_point; iPoint++) {
      if (ListBadPoint(iPoint) == 1) {
	std::vector<int> ListAdj=eGLP.first.Adjacency(iPoint);
	for (auto & eVal : ListAdj)
	  NewListBadPoint(eVal)=1;
      }
    }
    ListBadPoint=NewListBadPoint;
  }
  return ListBadPoint;
}


// Adapted from the matlab programs GRID_LinProgGetIJS_rx0.m and related
MyMatrix<double> DoLinearProgrammingSmoothing(GridArray const& GrdArr, double const& rx0max, int const& NeighborLevel)
{
  std::pair<GraphSparseImmutable, std::vector<std::pair<int,int>>> eGLP = GetGraphSparseVertexAdjacency(GrdArr);
  int nb_point = eGLP.second.size();
  MyMatrix<int> ListBadPoint=GetBadPoints(GrdArr, rx0max, NeighborLevel);
  std::vector<int> eList;
  for (int iPoint=0; iPoint<nb_point; iPoint++)
    if (ListBadPoint(iPoint,0) == 1)
      eList.push_back(iPoint);
  std::cerr << "|eList|=" << eList.size() << "\n";
  GraphSparseImmutable eGI = InducedSubgraph<GraphSparseImmutable,GraphSparseImmutable>(eGLP.first, eList);
  std::vector<std::vector<int>> ListConn=ConnectedComponents_set(eGI);
  std::cerr << "|ListConn|=" << ListConn.size() << "\n";
  MyVector<double> DEPwork(nb_point);
  for (int iPoint=0; iPoint<nb_point; iPoint++) {
    std::pair<int,int> ePair = eGLP.second[iPoint];
    DEPwork(iPoint) = GrdArr.GrdArrRho.DEP(ePair.first, ePair.second);
  }
  for (auto & eConn : ListConn) {
    int sizConn=eConn.size();
    std::cerr << "-------------------------------------------------------------------\n";
    std::cerr << "sizConn=" << sizConn << "\n";
    std::vector<int> ListMap(nb_point,-1);
    for (int i=0; i<sizConn; i++) {
      int ePt=eList[eConn[i]];
      ListMap[ePt]=i;
    }
    std::vector<MyVector<int>> ListPair;
    for (int i=0; i<sizConn; i++) {
      int ePt=eList[eConn[i]];
      std::vector<int> ListAdj=eGLP.first.Adjacency(ePt);
      //      std::cerr << "ePt=" << ePt << " |ListAdj|=" << ListAdj.size() << "\n";
      for (auto & eAdj : ListAdj) {
	int j=ListMap[eAdj];
	if (j != -1) {
	  MyVector<int> eVect(2);
	  eVect(0)=ePt;
	  eVect(1)=eAdj;
	  ListPair.push_back(eVect);
	}
      }
    }
    int nbPair=ListPair.size();
    //    std::cerr << " nbPair=" << nbPair << "\n";
    typedef Eigen::Triplet<double> T2;
    std::vector<double> ListVal;
    std::vector<T2> ListTriplet;
    int idx=0;
    double r=rx0max;
    for (auto & ePair : ListPair) {
      int ePt=ePair(0);
      int fPt=ePair(1);
      double dep1=GrdArr.GrdArrRho.DEP(ePt,0);
      double dep2=GrdArr.GrdArrRho.DEP(fPt,0);
      double CST=(-1-r)*dep1 + (1-r)*dep2;
      ListVal.push_back(CST);
      int ePos=ListMap[ePt];
      int fPos=ListMap[fPt];
      T2 eTr(idx, ePos, 1+r);
      T2 fTr(idx, fPos, -1+r);
      ListTriplet.push_back(eTr);
      ListTriplet.push_back(fTr);
      idx++;
    }
    // second set of inequalities xi <= h_i and -x_i <= h_i
    for (int i=0; i<sizConn; i++) {
      ListVal.push_back(0);
      T2 eTr1(idx, sizConn+i,1);
      T2 eTr2(idx, i, -1);
      ListTriplet.push_back(eTr1);
      ListTriplet.push_back(eTr2);
      idx++;
      ListVal.push_back(0);
      T2 eTr3(idx, sizConn+i,1);
      T2 eTr4(idx, i, 1);
      ListTriplet.push_back(eTr3);
      ListTriplet.push_back(eTr4);
      idx++;
    }
    int dimProb=2*sizConn;
    int nbIneq=idx;
    std::cerr << "nbPair=" << nbPair << " dimProb=" << dimProb << " nbIneq=" << nbIneq << "\n";
    MySparseMatrix<double> Aspmat(nbIneq, dimProb);
    std::cerr << "We have Aspmat\n";
    Aspmat.setFromTriplets(ListTriplet.begin(), ListTriplet.end());
    std::cerr << "Aspmat insert ListTriplets\n";
    MyVector<double> ListAconst=VectorFromStdVector(ListVal);
    std::cerr << "We have ListAconst\n";
    //
    int nbEqua=0;
    MySparseMatrix<double> Bspmat(nbEqua, dimProb);
    MyVector<double> ListBconst(nbEqua);
    std::cerr << "We have ListBconst\n";
    //
    MyVector<double> ToBeMinimized(dimProb);
    for (int i=0; i<sizConn; i++)
      ToBeMinimized(i)=0;
    for (int i=0; i<sizConn; i++)
      ToBeMinimized(sizConn+i)=1;
    std::cerr << "We have ToBeMinimized\n";
    GLPKoption eGLPKoption;
    eGLPKoption.UseDouble=true;
    eGLPKoption.UseExact=false;
    eGLPKoption.UseXcheck=false;
    std::cerr << "Before GLPJ_LinearProgramming_Kernel\n";
    LpSolutionSimple<double> eSol = GLPK_LinearProgramming_Kernel_Sparse_PROC(Aspmat, ListAconst, Bspmat, ListBconst, ToBeMinimized, eGLPKoption);
    std::cerr << "After GLPJ_LinearProgramming_Kernel\n";
    int nbRow=eSol.DirectSolution.rows();
    if (nbRow == 0) {
      std::cerr << "number of rows in DirectSolution is not correct\n";
      throw TerminalException{1};
    }
    for (auto & ePair : ListPair) {
      int ePt=ePair(0);
      int fPt=ePair(1);
      int ePos=ListMap[ePt];
      int fPos=ListMap[fPt];
      double x1=eSol.DirectSolution(ePos);
      double x2=eSol.DirectSolution(fPos);
      double h1=GrdArr.GrdArrRho.DEP(ePt,0);
      double h2=GrdArr.GrdArrRho.DEP(fPt,0);
      double CST=(-1-r)*h1 + (1-r)*h2;
      double dep1=h1 + x1;
      double dep2=h2 + x2;
      double delta=(1+r)*dep1 + (-1+r)*dep2;
      if (delta < 0) {
	std::cerr << "Error in our construction at ePt=" << ePt << " fPt=" << fPt << "\n";
	std::cerr << "ePos=" << ePos << " fPos=" << fPos << "\n";
	std::cerr << "dep1=" << dep1 << " dep2=" << dep2 << "\n";
	std::cerr << "  x1=" << x1   << "   x2=" << x2 << "\n";
	std::cerr << "  h1=" << h1   << "   h2=" << h2 << "\n";
	std::cerr << "delta=" << delta << " CST=" << CST << "\n";
      }
    }

    for (int i=0; i<sizConn; i++) {
      int ePt=eList[eConn[i]];
      double eVal=eSol.DirectSolution(i);
      //      std::cerr << "i=" << i << " eVal=" << eVal << "\n";
      DEPwork(ePt) += eVal;
    }
    std::cerr << "After assignation\n";
  }
  MyMatrix<double> DEPret = GrdArr.GrdArrRho.DEP;
  for (int iPoint=0; iPoint<nb_point; iPoint++) {
    std::pair<int,int> ePair = eGLP.second[iPoint];
    DEPret(ePair.first, ePair.second) = DEPwork(iPoint);
  }
  return DEPret;
}


MyMatrix<double> DoMartinhoBatteenSmoothing(GridArray const& GrdArr, double const& rx0max)
{
  std::pair<GraphSparseImmutable, std::vector<std::pair<int,int>>> eGLP = GetGraphSparseVertexAdjacency(GrdArr);
  MyMatrix<double> DEPret=GrdArr.GrdArrRho.DEP;
  int nb_point = eGLP.second.size();
  double r=rx0max;
  double eFact=(1-r)/(1+r);
  double epsilon=0.000001;
  while(true) {
    int nbOperation=0;
    for (int iPoint=0; iPoint<nb_point; iPoint++) {
      std::pair<int,int> ePair = eGLP.second[iPoint];
      double minDepAllowed=DEPret(ePair.first, ePair.second)*eFact;
      std::vector<int> ListAdj=eGLP.first.Adjacency(iPoint);
      for (auto & eAdj : ListAdj) {
        std::pair<int,int> fPair = eGLP.second[eAdj];
	double dep=DEPret(fPair.first, fPair.second);
	if (dep < minDepAllowed - epsilon) {
	  nbOperation++;
	  DEPret(fPair.first, fPair.second)=minDepAllowed;
	}
      }
    }
    if (nbOperation == 0)
      break;
  }
  return DEPret;
}




FullNamelist NAMELIST_GetStandard_Bathymetry_Smoothing()
{
  std::map<std::string, SingleBlock> ListBlock;
  // PROC
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string> > ListListStringValues1;
  std::map<std::string, std::vector<int> > ListListIntValues1;
  std::map<std::string, std::vector<double> > ListListDoubleValues1;
  ListStringValues1["GridFileIn"]="unset";
  ListStringValues1["GridFileOut"]="unset";
  ListStringValues1["Method"]="unset";
  ListDoubleValues1["rx0max"]=1;
  ListIntValues1["NeighborLevel"]=5;
  SingleBlock BlockPROC;
  BlockPROC.ListIntValues=ListIntValues1;
  BlockPROC.ListBoolValues=ListBoolValues1;
  BlockPROC.ListDoubleValues=ListDoubleValues1;
  BlockPROC.ListStringValues=ListStringValues1;
  BlockPROC.ListListStringValues=ListListStringValues1;
  BlockPROC.ListListIntValues=ListListIntValues1;
  BlockPROC.ListListDoubleValues=ListListDoubleValues1;
  ListBlock["PROC"]=BlockPROC;
  // Final part
  return {ListBlock, "undefined"};
}





void DoFullSmoothing(FullNamelist const& eFull)
{
  //
  // Loading the grid
  //
  SingleBlock eBlMODEL=eFull.ListBlock.at("MODEL");
  std::string eModelName=eBlMODEL.ListStringValues.at("MODELNAME");
  std::string GridFile=eBlMODEL.ListStringValues.at("GridFile");
  std::string HisPrefix=eBlMODEL.ListStringValues.at("HisPrefix");
  std::string BoundFile="unset";
  TripleModelDesc eTriple{eModelName, GridFile, BoundFile, HisPrefix, {}};
  GridArray GrdArr=RETRIEVE_GRID_ARRAY(eTriple);
  //
  // Doing the conversion
  //
  SingleBlock eBlPROC=eFull.ListBlock.at("PROC");
  std::string TheMethod=eBlPROC.ListStringValues.at("Method");
  std::string GridFileOut = eBlPROC.ListStringValues.at("GridFileOut");
  double rx0max=eBlPROC.ListDoubleValues.at("rx0max");
  int NeighborLevel=eBlPROC.ListIntValues.at("NeighborLevel");
  MyMatrix<double> RMat1=GetRoughnessFactor(GrdArr.GrdArrRho.DEP, GrdArr);
  std::cerr << " rx0(input)=" << RMat1.maxCoeff() << "\n";
  MyMatrix<double> DEPnew;
  bool DoOper=false;
  if (TheMethod == "LinearProgramming") {
    DEPnew = DoLinearProgrammingSmoothing(GrdArr, rx0max, NeighborLevel);
    DoOper=true;
  }
  if (TheMethod == "MartinhoBatteen") {
    DEPnew = DoMartinhoBatteenSmoothing(GrdArr, rx0max);
    DoOper=true;
  }
  if (!DoOper) {
    std::cerr << "Failed to find a matching function for bathymetry smoothing\n";
    std::cerr << "Please correct\n";
    throw TerminalException{1};
  }
  GrdArr.GrdArrRho.DEP=DEPnew;
  MyMatrix<double> RMat2=GetRoughnessFactor(DEPnew, GrdArr);
  std::cerr << "rx0(output)=" << RMat2.maxCoeff() << "\n";
  //
  // Now writing the grid
  //
  WriteGrid(GridFileOut, GrdArr);
}











#endif
