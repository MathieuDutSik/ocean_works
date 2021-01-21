#include "Model_grids.h"
int main(int argc, char *argv[])
{
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
  try {
    if (argc != 3) {
      std::cerr << "ShowGridInformation is used as\n";
      std::cerr << "ShowGridInformation [GridFileIN] [BoundFileIN]\n";
      std::cerr << "with GridFileIN    the input grid\n";
      std::cerr << "with BoundFileIN   the input boundary (put unset if not available)\n";
      return -1;
    }
    std::string GridFileIN  = argv[1];
    std::string BoundFileIN = argv[2];
    std::cerr << " GridFileIN = " << GridFileIN << "\n";
    std::cerr << "BoundFileIN = " << BoundFileIN << "\n";
    GridArray GrdArr=ReadUnstructuredGrid(GridFileIN, BoundFileIN);
    int mnp=GrdArr.GrdArrRho.LON.rows();
    int mne=GrdArr.INE.rows();
    std::cerr << "mnp=" << mnp << " mne=" << mne << "\n";
    //
    // longitude / latitude
    //
    double minLon=GrdArr.GrdArrRho.LON.minCoeff();
    double maxLon=GrdArr.GrdArrRho.LON.maxCoeff();
    double minLat=GrdArr.GrdArrRho.LAT.minCoeff();
    double maxLat=GrdArr.GrdArrRho.LAT.maxCoeff();
    std::cerr << "LON(min/max)=" << minLon << " / " << maxLon << " LAT(min/max)=" << minLat << " / " << maxLat << "\n";
    //
    // Vertex degrees
    //
    GraphSparseImmutable GR = GetUnstructuredVertexAdjInfo(GrdArr.INE, mnp);
    std::vector<int> ListDeg(mnp);
    for (int iVert=0; iVert<mnp; iVert++)
      ListDeg[iVert] = int(GR.Adjacency(iVert).size());
    std::cerr << "Degree of the vertices\n";
    ShowAttainmentVector(std::cerr, ListDeg);
    //
    // Boundary values.
    //
    int sizIOBP=GrdArr.IOBP.size();
    if (sizIOBP > 0) {
      std::vector<int> IOBP_vect=StdVectorFromVector(GrdArr.IOBP);
      std::cerr << "The IOBP values are shown\n";
      ShowAttainmentVector(std::cerr, IOBP_vect);
    }
    else {
      std::cerr << "No boundary avalaible\n";
    }
    //
    // Boundaries
    //
    std::vector<std::vector<int>> ListListBnd = GetListBoundaryCycles(GrdArr.INE, mnp);
    int nbBnd=ListListBnd.size();
    std::cerr << "|ListListBnd|=" << nbBnd << "\n";
    std::vector<int> ListLenBnd;
    std::vector<int> ListStatus(mnp,0);
    std::vector<double> ListLenIsland;
    for (auto & eListBnd : ListListBnd) {
      ListLenBnd.push_back(int(eListBnd.size()));
      for (auto & eVal : eListBnd)
	ListStatus[eVal] = 1;
      int len=eListBnd.size();
      double sumLenKM=0;
      for (int i=0; i<len; i++) {
	int j=NextIdx(len,i);
	int eNode1=eListBnd[i];
	int eNode2=eListBnd[j];
	double eLon1=GrdArr.GrdArrRho.LON(eNode1,0);
	double eLat1=GrdArr.GrdArrRho.LAT(eNode1,0);
	double eLon2=GrdArr.GrdArrRho.LON(eNode2,0);
	double eLat2=GrdArr.GrdArrRho.LAT(eNode2,0);
	double distKM = GeodesicDistanceKM(eLon1, eLat1, eLon2, eLat2);
	sumLenKM += distKM;
      }
      ListLenIsland.push_back(sumLenKM);
    }
    //
    // Determination of the main boundary
    //
    int MaxLenBnd = VectorMax(ListLenBnd);
    double TotSumIslandKM=0;
    int iBndMax=-1;
    for (int iBnd=0; iBnd<nbBnd; iBnd++) {
      if (ListLenBnd[iBnd] < MaxLenBnd) {
	TotSumIslandKM += ListLenIsland[iBnd];
      }
      else {
	iBndMax = iBnd;
      }
    }
    double avgLenIslandKM = TotSumIslandKM / double(nbBnd-1);
    std::cerr << "Number of boundary points = " << VectorSum(ListStatus) << "\n";
    std::cerr << "Minimum length island=" << VectorMin(ListLenIsland) << " avg length island=" << avgLenIslandKM << "\n";
    //
    if (sizIOBP > 0) {
      int lenLandBound=ListListBnd[iBndMax].size();
      double sumLandBndKM=0;
      for (int iLandBound=0; iLandBound<lenLandBound; iLandBound++) {
	int jLandBound=NextIdx(lenLandBound, iLandBound);
	int i1=ListListBnd[iBndMax][iLandBound];
	int i2=ListListBnd[iBndMax][jLandBound];
	if (GrdArr.IOBP(i1) != 2 || GrdArr.IOBP(i2) != 2) {
	  double eLon1=GrdArr.GrdArrRho.LON(i1,0);
	  double eLat1=GrdArr.GrdArrRho.LAT(i1,0);
	  double eLon2=GrdArr.GrdArrRho.LON(i2,0);
	  double eLat2=GrdArr.GrdArrRho.LAT(i2,0);
	  double distKM = GeodesicDistanceKM(eLon1, eLat1, eLon2, eLat2);
	  sumLandBndKM += distKM;
	}
      }
      std::cerr << "lenLandBound=" << lenLandBound << " Length of land boundary=" << sumLandBndKM << "\n";
    }
    //
    std::cerr << "List of boundary length\n";
    std::cerr << "MaxLen=" << VectorMax(ListLenBnd) << " MinLen=" << VectorMin(ListLenBnd) << "\n";
    //
    std::cerr << "Normal termination of ShowGridInformation\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Error in ShowGridInformation\n";
    exit(e.eVal);
  }
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "runtime = " << std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count() << "\n";
}
