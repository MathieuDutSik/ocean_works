#include "Model_grids.h"

int main(int argc, char *argv[])
{
  srand_random_set();
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
  try {
    if (argc != 4) {
      std::cerr << "GRID_SetStraightLine [GridFile] [SetVal] [WWMboundary]\n";
      std::cerr << "\n";
      std::cerr << "GridFile a WWM grid file (type determined by suffix\n";
      std::cerr << "SetVal a value to be set\n";
      std::cerr << "WWMboundary is the output file\n";
      return -1;
    }
    std::string GridFile=argv[1];
    std::cerr << "GridFile = " << GridFile << "\n";
    //
    int SetVal=-400;
    sscanf(argv[2], "%d", &SetVal);
    //
    std::string BndFile=argv[3];
    //
    GridArray GrdArr=ReadUnstructuredGrid(GridFile, "unset");
    int nbNode=GrdArr.GrdArrRho.LON.rows();
    MyVector<int> Status=GetBoundaryStatus(GrdArr.INE, nbNode);
    int st0=0, st1=0, stM1=0;
    for (int iNode=0; iNode<nbNode; iNode++) {
      if (Status(iNode) == 0) st0++;
      if (Status(iNode) == 1) st1++;
      if (Status(iNode) == -1) stM1++;
    }
    std::cerr << "st0=" << st0 << " st1=" << st1 << " stM1=" << stM1 << "\n";

    MyVector<int> IOBP(nbNode);
    std::vector<int> LBnd;
    MyMatrix<double> LON = GrdArr.GrdArrRho.LON;
    MyMatrix<double> LAT = GrdArr.GrdArrRho.LAT;
    double minLON=LON.minCoeff();
    double maxLON=LON.maxCoeff();
    double minLAT=LAT.minCoeff();
    double maxLAT=LAT.maxCoeff();
    std::cerr << "LON(min/max)=" << minLON << " / " << maxLON << "\n";
    std::cerr << "LAT(min/max)=" << minLAT << " / " << maxLAT << "\n";
    double eX0 = 18.33;
    double eY0 = 40.29;
    double eX1 = 19.55;
    double eY1 = 41.1;
    double line_slope = (eY1 - eY0) / (eX1 - eX0);
    double line_cst = eY0 - eX0 * line_slope;
    for (int iNode=0; iNode<nbNode; iNode++) {
      int eIOBP;
      if (Status(iNode) == -1)
	eIOBP=1;
      else
	eIOBP=0;
      IOBP(iNode)=eIOBP;
      double eX=LON(iNode,0);
      double eY=LAT(iNode,0);
      double Yvalue=line_cst + eX * line_slope;
      if (eIOBP == 1 && eY < 41.5 && eX > 18 && eY < Yvalue)
	LBnd.push_back(iNode);
    }
    int nbBnd=LBnd.size();
    std::cerr << "nbBnd=" << nbBnd << "\n";
    //    int iNodeSel=-1;
    //    int jNodeSel=-1;
    int nbMatch=0;
    double tolLL=0.0005;
    std::vector<int> SelectedStraight;
    for (int iBnd=0; iBnd<nbBnd; iBnd++)
      for (int jBnd=iBnd+1; jBnd<nbBnd; jBnd++) {
	std::cerr << "iBnd=" << iBnd << " jBnd=" << jBnd << " nbMatch=" << nbMatch << "\n";
	int iNode=LBnd[iBnd];
	int jNode=LBnd[jBnd];
	//            std::cerr << "iNode=" << iNode << "\n";
	double eX=LON(iNode,0);
	double eY=LAT(iNode,0);
	double fX=LON(jNode,0);
	double fY=LAT(jNode,0);
	//      std::cerr << "After eX, eY\n";
	double deltaX = eX - fX;
	double deltaY = eY - fY;
	double dist=fabs(deltaX) + fabs(deltaY);
	//      std::cerr << "Before test\n";
	if (dist > 0.2) {
	  std::vector<int> TheSel;
	  //	std::cerr << "Before loop\n";
	  for (int kBnd=0; kBnd<nbBnd; kBnd++) {
	    int kNode=LBnd[kBnd];
	    //	  std::cerr << "kNode=" << kNode << "\n";
	    double gX=LON(kNode,0);
	    double gY=LAT(kNode,0);
	    //	  std::cerr << "We have gX, gY\n";
	    double det=deltaX * (gY - eY) - deltaY * (gX - eX);
	    if (fabs(det) < tolLL)
	      TheSel.push_back(kNode);
	  }
	  //	std::cerr << "After loop\n";
	  int siz=TheSel.size();
	  if (siz > nbMatch) {
	    nbMatch = siz;
	    SelectedStraight = TheSel;
	  }
	}
	//      std::cerr << "After test\n";
      }
    std::cerr << "nbMatch=" << nbMatch << "\n";
    for (auto & iNode : SelectedStraight) {
      IOBP(iNode) = SetVal;
      double eX = LON(iNode,0);
      double eY = LAT(iNode,0);
      std::cerr << "iNode=" << iNode << " eX=" << eX << " eY=" << eY << "\n";
    }
    GrdArr.IOBP=IOBP;
    if (IsExistingFile(BndFile) == true) {
      std::cerr << "Error, please remove the file\n";
      std::cerr << "BndFile = " << BndFile << "\n";
      std::cerr << "or correct\n";
      throw TerminalException{1};
    }
    WriteWWMboundaryGR3(BndFile, GrdArr);
    int minIOBP=IOBP.minCoeff();
    int maxIOBP=IOBP.maxCoeff();
    size_t len=1 + maxIOBP - minIOBP;
    std::vector<int> ListNbMatch(len, 0);
    for (int iNode=0; iNode<nbNode; iNode++) {
      int eIOBP=IOBP(iNode);
      int pos=eIOBP - minIOBP;
      ListNbMatch[pos]++;
    }
    for (int eIOBP=minIOBP; eIOBP<=maxIOBP; eIOBP++) {
      int pos=eIOBP - minIOBP;
      int eNb=ListNbMatch[pos];
      std::cerr << " eIOBP=" << eIOBP << " eNb=" << eNb << "\n";
    }
    std::cerr << "Normal termination of GRID_SetStraightLine\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Error in GRID_SetStraightLine\n";
    exit(e.eVal);
  }
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "runtime = " << std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count() << "\n";
}
