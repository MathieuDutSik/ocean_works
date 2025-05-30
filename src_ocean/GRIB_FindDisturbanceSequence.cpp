// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Model_grids.h"
#include "NamelistExampleOcean.h"
int main(int argc, char *argv[]) {
  srand_random_set();
  HumanTime time1;
  try {
    FullNamelist eFull = NAMELIST_ComparisonSequentialRuns();
    if (argc != 2) {
      std::cerr << "GRIB_PrintDisturbanceSequence [file.nml]\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    //
    SingleBlock BlPROC = eFull.get_block("PROC");
    std::string HisPrefix = BlPROC.get_string("HisPrefix");
    std::string ModelName = BlPROC.get_string("ModelName");
    std::string shortName = BlPROC.get_string("shortName");
    //
    std::cerr << "ModelName=" << ModelName << "\n";
    bool RetAllStates = RetrieveAllStates(ModelName);
    if (!RetAllStates) {
      std::cerr << "We need to have \"retrieveallstates\" in the ModelName\n";
      throw TerminalException{1};
    }
    std::string GridFile = "unset"; // for grib there is no separate grid
    TripleModelDesc eTriple{ModelName, GridFile, "unset", HisPrefix, {}};
    GridArray GrdArr = RETRIEVE_GRID_ARRAY(eTriple);
    ArrayHistory eArr = ReadArrayHistory(eTriple);
    TotalArrGetData TotalArr{GrdArr, eArr};
    //
    int GEOSELECTION = BlPROC.get_int("GEOSELECTION");
    double MinLON = BlPROC.get_double("MinLON");
    double MaxLON = BlPROC.get_double("MaxLON");
    double MinLAT = BlPROC.get_double("MinLAT");
    double MaxLAT = BlPROC.get_double("MaxLAT");
    std::vector<double> LONPOLY = BlPROC.get_list_double("LONPOLY");
    std::vector<double> LATPOLY = BlPROC.get_list_double("LATPOLY");
    std::vector<std::pair<int, int>> ListPIdx;
    int nbRow = GrdArr.GrdArrRho.LON.rows();
    int nbCol = GrdArr.GrdArrRho.LON.cols();
    for (int iRow = 0; iRow < nbRow; iRow++) {
      for (int iCol = 0; iCol < nbCol; iCol++) {
        double eLon = GrdArr.GrdArrRho.LON(iRow, iCol);
        double eLat = GrdArr.GrdArrRho.LAT(iRow, iCol);
        bool eStatusGeo = true;
        if (GEOSELECTION == 1) {
          if (!(eLon >= MinLON && eLon <= MaxLON && eLat >= MinLAT &&
                eLat <= MaxLAT))
            eStatusGeo = false;
        }
        if (GEOSELECTION == 2) {
          bool IsInside = IsPointInside(eLon, eLat, LONPOLY, LATPOLY);
          if (!IsInside)
            eStatusGeo = false;
        }
        if (eStatusGeo)
          ListPIdx.push_back({iRow, iCol});
      }
    }
    int nbPIdx = ListPIdx.size();
    //
    int nbTime = eArr.FullOrganizedInfo[shortName].size();
    std::vector<double> ListSpreading(nbTime);
    for (int iTime = 0; iTime < nbTime; iTime++) {
      std::cerr << "iTime=" << iTime << " / " << nbTime << "\n";
      std::pair<double, std::vector<GRIB_MessageInfo>> ePair =
          eArr.FullOrganizedInfo[shortName][iTime];
      std::vector<MyVector<double>> ListVector;
      for (auto const &eMesg : ePair.second) {
        MyMatrix<double> eMat = GRIB_ReadFromMessageInfo(eMesg);
        MyVector<double> V(nbPIdx);
        for (int iP = 0; iP < nbPIdx; iP++) {
          int iRow = ListPIdx[iP].first;
          int iCol = ListPIdx[iP].second;
          V(iP) = eMat(iRow, iCol);
        }
        ListVector.push_back(V);
      }
      MyVector<double> Vsum = ZeroVector<double>(nbPIdx);
      for (auto &eV : ListVector)
        Vsum += eV;
      MyVector<double> Vavg = Vsum / nbPIdx;
      double SumError = 0;
      for (auto &eV : ListVector) {
        MyVector<double> Vdiff = eV - Vavg;
        Vdiff = Vdiff.cwiseAbs();
        double eSum = Vdiff.sum();
        SumError += eSum;
      }
      double AvgError = SumError / (nbPIdx * ListVector.size());
      ListSpreading[iTime] = AvgError;
    }
    double TheAvg = VectorSum(ListSpreading) / nbTime;
    std::cerr << "TheAvg = " << TheAvg << "\n";
    std::cerr << "Normal termination of GRIB_FindDisturbanceSequence\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRIB_FindDisturbanceSequence\n";
    exit(e.eVal);
  }
  runtime(time1);
}
