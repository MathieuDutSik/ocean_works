#include "Model_grids.h"

int main(int argc, char *argv[])
{
  try {
    if (argc != 12) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "ComputeWWMboundary [GridFile] [NorthStat] [WestStat] [SouthStat] [EastStat] [NorthVal] [WestVal] [SouthVal] [EastVal] [tolLL] [WWMboundary]\n";
      std::cerr << "\n";
      std::cerr << "GridFile a WWM grid file (type determined by suffix\n";
      std::cerr << "NorthStat, WestStat, SouthStat and EastStat specifies boundary status\n";
      std::cerr << "tolLL : tolerance for setting up the boundary values\n";
      std::cerr << "WWMboundary is the output file\n";
      return -1;
    }
    std::string GridFile=argv[1];
    std::cerr << "GridFile = " << GridFile << "\n";
    //
    int NorthStat=-400, WestStat=-400, SouthStat=-400, EastStat=-400;
    int NorthVal=-400, WestVal=-400, SouthVal=-400, EastVal=-400;
    sscanf(argv[2], "%d", &NorthStat);
    sscanf(argv[3], "%d", &WestStat);
    sscanf(argv[4], "%d", &SouthStat);
    sscanf(argv[5], "%d", &EastStat);
    sscanf(argv[6], "%d", &NorthVal);
    sscanf(argv[7], "%d", &WestVal);
    sscanf(argv[8], "%d", &SouthVal);
    sscanf(argv[9], "%d", &EastVal);
    std::cerr << "NorthStat = " << NorthStat << " argv=" << argv[2] << " NorthVal=" << NorthVal << "\n";
    std::cerr << "WestStat  = " << WestStat  << " argv=" << argv[3] << "  WestVal=" << WestVal  << "\n";
    std::cerr << "SouthStat = " << SouthStat << " argv=" << argv[4] << " SouthVal=" << SouthVal << "\n";
    std::cerr << "EastStat  = " << EastStat  << " argv=" << argv[5] << "  EastVal=" << EastVal  << "\n";
    //
    double tolLL;
    sscanf(argv[10], "%lf", &tolLL);
    //
    std::string BndFile=argv[11];
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
    for (int iNode=0; iNode<nbNode; iNode++) {
      int eIOBP;
      if (Status(iNode) == -1)
	eIOBP=1;
      else
	eIOBP=0;
      IOBP(iNode)=eIOBP;
    }
    //
    std::vector<double> ListLON(nbNode);
    std::vector<double> ListLAT(nbNode);
    for (int iNode=0; iNode<nbNode; iNode++) {
      ListLON[iNode]=GrdArr.GrdArrRho.LON(iNode,0);
      ListLAT[iNode]=GrdArr.GrdArrRho.LAT(iNode,0);
    }
    double MinLon=VectorMin(ListLON);
    double MaxLon=VectorMax(ListLON);
    double MinLat=VectorMin(ListLAT);
    double MaxLat=VectorMax(ListLAT);
    int nbSouth=0, nbNorth=0, nbWest=0, nbEast=0;
    for (int iNode=0; iNode<nbNode; iNode++) {
      double eLon=ListLON[iNode];
      double eLat=ListLAT[iNode];
      if (IOBP(iNode) == 1) {
	if (eLon < MinLon + tolLL) {
	  nbWest++;
	  if (WestStat == 1)
	    IOBP(iNode)=WestVal;
	}
	if (eLon > MaxLon - tolLL) {
	  nbEast++;
	  if (EastStat == 1)
	    IOBP(iNode)=EastVal;
	}
	if (eLat < MinLat + tolLL) {
	  nbSouth++;
	  if (SouthStat == 1)
	    IOBP(iNode)=SouthVal;
	}
	if (eLat > MaxLat - tolLL) {
	  nbNorth++;
	  if (NorthStat == 1)
	    IOBP(iNode)=NorthVal;
	}
      }
    }
    std::cerr << "nbSouth=" << nbSouth << " nbNorth=" << nbNorth << " nbEast=" << nbEast << " nbWest=" << nbWest << "\n";
    /*  I do not know why I wrote that segment of code
	for (int iNode=0; iNode<nbNode; iNode++) {
	if (IOBP(iNode) == 1)
	IOBP(iNode)=0;
	}*/
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
    std::cerr << "Normal termination of ComputeWWMboundary\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
