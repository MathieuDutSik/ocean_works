// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_OCEAN_TRANSECT_H_
#define SRC_OCEAN_TRANSECT_H_

#include "Basic_plot.h"
#include "CommonFuncModel.h"
#include "Interpolation.h"
#include "Kernel_Transect.h"
#include "Model_grids.h"
#include "Model_interpolation.h"
#include "NamelistExampleOcean.h"
#include "SphericalGeom.h"
#include <string>
#include <vector>

TransectInformation_3D
RetrievePointTransectRecord(TotalArrGetData const &TotalArr,
                            MyMatrix<double> const &ListXY,
                            double const &VertResolM) {
  Eigen::Tensor<double, 3> VertCoord =
      RetrieveStandardVerticalCoordinate(TotalArr);
  SingleArrayInterpolationGen eInterp = {
      ComputeArrayInterpolation_ListXY(TotalArr.GrdArr, ListXY), {}};

  std::vector<PairLL> ListPairLL;
  int nbPoint = ListXY.cols();
  for (int iPt = 0; iPt < nbPoint; iPt++) {
    double eLon = ListXY(0, iPt);
    double eLat = ListXY(1, iPt);
    ListPairLL.push_back({eLon, eLat});
  }
  TransectInformation eTrans;
  eTrans.ListPairLL = ListPairLL;
  eTrans.ListRec = {eInterp};
  return GetTransectInformation_3D(eTrans, TotalArr.GrdArr, VertCoord,
                                   VertResolM);
}

std::vector<PointOutTrans> ReadStationCoordinate(SingleBlock const &eBlPLOT) {
  std::string type = eBlPLOT.get_string("TypeListPoint");
  if (type == "fileRovinj") {
    std::string eFile = eBlPLOT.get_string("ListPointFile");
    return ReadStationCoordinate_File(eFile);
  }
  if (type == "namelist") {
    std::vector<double> ListPointLon =
        eBlPLOT.get_list_double("ListPointLon");
    std::vector<double> ListPointLat =
        eBlPLOT.get_list_double("ListPointLat");
    std::vector<double> ListPointDep =
        eBlPLOT.get_list_double("ListPointDepth");
    std::vector<std::string> ListPointName =
        eBlPLOT.get_list_string("ListPointName");
    int nbEnt = ListPointLon.size();
    std::vector<PointOutTrans> ListPoint;
    for (int iEnt = 0; iEnt < nbEnt; iEnt++) {
      double lon = ListPointLon[iEnt];
      double lat = ListPointLat[iEnt];
      double dep = ListPointDep[iEnt];
      std::string name = ListPointName[iEnt];
      ListPoint.push_back({lon, lat, dep, name});
    }
    return ListPoint;
  }
  if (type == "empty")
    return {};
  std::cerr << "Possible values for TypeListPoint are \n";
  std::cerr << "empty, namelist and fileRovinj\n";
  std::cerr << "TypeListPoint = " << type << "\n";
  throw TerminalException{1};
}

void SetDefaultDrawLinesArr(DrawLinesArr &eDrawArr) {
  eDrawArr.IsTimeSeries = false;
  eDrawArr.PairComparison = false;
  eDrawArr.DoExplicitLabel = false;
  // Maybe DrawHorizVertLines should be put as input parameter.
  eDrawArr.DrawHorizVertLines = false;
}

void WriteDrawArrFile(std::string const& FileName, DrawLinesArr const& eDrawArr) {
  int nbArr = eDrawArr.ListListVect.size();
  int nbEntry = eDrawArr.ListListVect[0].size();
  //
  std::string FileOut = FileName + ".txt";
  std::ofstream os(FileOut);
  if (eDrawArr.IsTimeSeries)
    os << "ListTime";
  else
    os << "ListX";
  for (int iArr=0; iArr<nbArr; iArr++)
    os << "; " << eDrawArr.ListName_plot[iArr];
  os << "\n";
  //
  for (int iEntry=0; iEntry<nbEntry; iEntry++) {
    if (eDrawArr.IsTimeSeries && !eDrawArr.DoExplicitLabel) {
      std::string strPres = DATE_ConvertMjd2mystringPres(eDrawArr.ListX(iEntry));
      os << strPres;
    } else {
      os << eDrawArr.ListX(iEntry);
    }
    for (int iArr=0; iArr<nbArr; iArr++)
      os << "; " << eDrawArr.ListListVect[iArr](iEntry);
    os << "\n";
  }
}


void TRANSECT_Plot(FullNamelist const &eFull) {
  SingleBlock eBlPLOT = eFull.get_block("PLOT");
  //
  // Reading grid arrays and the like
  //
  PermanentInfoDrawing ePerm = GET_PERMANENT_INFO(eFull);
  NCLcaller<GeneralType> eCall(ePerm.NPROC);
  //
  SingleBlock eBlPROC = eFull.get_block("PROC");
  std::vector<std::string> ListModelName =
      eBlPROC.get_list_string("ListMODELNAME");
  std::vector<std::string> ListGridFile =
      eBlPROC.get_list_string("ListGridFile");
  std::vector<std::string> ListHisPrefix =
      eBlPROC.get_list_string("ListHisPrefix");
  std::vector<std::string> ListRunName =
      eBlPROC.get_list_string("ListRunName");
  bool PrintTextFiles = eBlPROC.get_bool("PrintTextFiles");
  int nbGrid = ListGridFile.size();
  size_t nbGrid_t = ListGridFile.size();
  std::cerr << "nbGrid=" << nbGrid << "\n";
  if (nbGrid_t != ListHisPrefix.size() || nbGrid_t != ListModelName.size()) {
    std::cerr << "Error for the length of\n";
    std::cerr << "ListMODELNAME, ListGridFile, ListHisPrefix\n";
    std::cerr << "which should be of the same length\n";
    throw TerminalException{1};
  }
  std::vector<GridArray> ListGrdArr(nbGrid);
  std::vector<ArrayHistory> ListArrayHistory(nbGrid);
  std::vector<TotalArrGetData> ListTotalArr(nbGrid);

  for (int iGrid = 0; iGrid < nbGrid; iGrid++) {
    std::cerr << "iGrid=" << iGrid << " / " << nbGrid << "\n";
    std::string eModelName = ListModelName[iGrid];
    std::string GridFile = ListGridFile[iGrid];
    std::string HisPrefix = ListHisPrefix[iGrid];
    TripleModelDesc eTriple{eModelName, GridFile, "unset", HisPrefix, {}};
    GridArray GrdArr = RETRIEVE_GRID_ARRAY(eTriple);
    ListGrdArr[iGrid] = GrdArr;
    std::cerr << "Before call to ReadArrayHistory\n";
    ArrayHistory eArr = ReadArrayHistory(eTriple);
    std::cerr << " After call to ReadArrayHistory\n";
    ListArrayHistory[iGrid] = eArr;
    TotalArrGetData TotalArr = RetrieveTotalArr(eTriple);
    ListTotalArr[iGrid] = TotalArr;
  }
  //
  // Reading vertical transect point information
  //
  bool VariableRange = eBlPLOT.get_bool("VariableRange");
  std::string VariableRangeRounding = eBlPLOT.get_string("VariableRangeRounding");
  bool DoTitle = eBlPLOT.get_bool("DoTitle");
  std::vector<PointOutTrans> ListPointOut = ReadStationCoordinate(eBlPLOT);
  double VertResolM = eBlPLOT.get_double("VertResolM");
  int nbPointOut = ListPointOut.size();
  std::cerr << "nbPointOut = " << nbPointOut << "\n";
  std::vector<TransectInformation_3D> ListTrans3D(nbGrid);
  std::vector<int> ListNbVert(nbPointOut);
  if (nbPointOut > 0) {
    MyMatrix<double> ListXY(2, nbPointOut);
    for (int iPt = 0; iPt < nbPointOut; iPt++) {
      ListXY(0, iPt) = ListPointOut[iPt].lon;
      ListXY(1, iPt) = ListPointOut[iPt].lat;
    }
    for (int iGrid = 0; iGrid < nbGrid; iGrid++)
      ListTrans3D[iGrid] =
          RetrievePointTransectRecord(ListTotalArr[iGrid], ListXY, VertResolM);
    for (int iPt = 0; iPt < nbPointOut; iPt++) {
      std::vector<int> ListNB(nbGrid);
      for (int iGrid = 0; iGrid < nbGrid; iGrid++) {
        int nbVert = ListTrans3D[iGrid].ListVertPos.size();
        int eFirstWet = ListTrans3D[iGrid].FirstWetIndex[iPt];
        int eNB = nbVert - eFirstWet;
        ListNB[iGrid] = eNB;
      }
      ListNbVert[iPt] = VectorMax(ListNB);
    }
  }
  //
  // Reading plotting information.
  //
  std::vector<TransectInformation> ListTransect =
      RetrieveListTransect(eBlPLOT, ListGrdArr);
  int nbTrans = ListTransect.size();
  std::cerr << "nbTrans = " << nbTrans << "\n";
  SingleBlock eBlockVAR = eFull.get_block("VARS");
  std::vector<std::string> ListVarOut = eBlockVAR.ExtractMatchingBool();
  VarQuery eQuery;
  std::vector<std::string> ListNatureQuery =
      eBlPROC.get_list_string("ListNatureQuery");
  eQuery.TimeFrameDay = eBlPROC.get_double("TimeFrameDay");
  std::vector<VarQuery> ListQuery =
      GetIntervalGen_Query(eBlPROC, ListArrayHistory);
  int nbTime = ListQuery.size();
  std::cerr << "nbTime=" << nbTime << "\n";
  PlotBound ePlotBound = ReadPlotBound(eFull);
  for (int iTime = 0; iTime < nbTime; iTime++) {
    VarQuery eQuery = ListQuery[iTime];
    std::string strPres = DATE_ConvertMjd2mystringPres(eQuery.eTimeDay);
    std::string strFile = DATE_ConvertMjd2mystringFile(eQuery.eTimeDay);
    std::cerr << "iTime=" << iTime << "/" << nbTime << " date=" << strPres
              << "\n";
    for (auto &eVarName : ListVarOut) {
      for (auto &eNatureQuery : ListNatureQuery) {
        eQuery.NatureQuery = eNatureQuery;
        std::vector<RecVar> ListRecVar(nbGrid);
        for (int iGrid = 0; iGrid < nbGrid; iGrid++)
          ListRecVar[iGrid] = ModelSpecificVarSpecificTimeGeneral(
              ListTotalArr[iGrid], eVarName, eQuery, ePlotBound);
        for (int iTrans = 0; iTrans < nbTrans; iTrans++) {
          std::string strSp, strUnder;
          if (nbTrans == 1) {
            strSp = "";
            strUnder = "";
          } else {
            strSp = " " + std::to_string(iTrans + 1);
            strUnder = "_" + std::to_string(iTrans + 1);
          }
          int nbPoint = ListTransect[iTrans].ListPairLL.size();
          DrawLinesArr eDrawArr;
          SetDefaultDrawLinesArr(eDrawArr);
          eDrawArr.DoTitle = DoTitle;
          eDrawArr.TitleStr =
              "Transect" + strSp + " of " + eVarName + " at " + strPres;
          std::string fVarName =
              "transect" + strUnder + "_" + eVarName + "_" + strFile;
          eDrawArr.VarName = fVarName;
          eDrawArr.ListName_plot = ListRunName;
          eDrawArr.YAxisString =
              ListRecVar[0].RecS.VarName2 + "(" + ListRecVar[0].RecS.Unit + ")";
          //
          // Determination of dimension variable
          //
          eDrawArr.ListX = ListTransect[iTrans].ListDimVar;
          //
          // Interpolation to the transect
          //
          std::vector<MyVector<double>> ListListVect(nbGrid);
          std::vector<double> ListMax(nbGrid);
          std::vector<double> ListMin(nbGrid);
          for (int iGrid = 0; iGrid < nbGrid; iGrid++) {
            RecVar eRecVar = INTERPOL_SingleRecVarInterpolation(
                ListTransect[iTrans].ListRec[iGrid], ListRecVar[iGrid]);
            MyVector<double> eVect(nbPoint);
            for (int iPoint = 0; iPoint < nbPoint; iPoint++) {
              eVect(iPoint) = eRecVar.F(iPoint, 0);
              // std::cerr << "iPoint=" << iPoint << " val=" <<
              // eRecVar.F(iPoint,0) << "\n";
            }
            ListListVect[iGrid] = eVect;
            ListMax[iGrid] = eVect.maxCoeff();
            ListMin[iGrid] = eVect.minCoeff();
          }
          //
          // Determinantion of min/max ranges
          //
          double TheMax, TheMin;
          if (!VariableRange) {
            TheMax = ListRecVar[0].RecS.maxval;
            TheMin = ListRecVar[0].RecS.minval;
          } else {
            TheMax = VectorMax(ListMax);
            TheMin = VectorMin(ListMin);
          }
          TheMax = ApplyRounding(TheMax, VariableRangeRounding);
          TheMin = ApplyRounding(TheMin, VariableRangeRounding);
          eDrawArr.TheMax = TheMax;
          eDrawArr.TheMin = TheMin;

          eDrawArr.ListListVect = ListListVect;
          std::string FileName =
              ePerm.eDir + "Transect_" + std::to_string(iTrans + 1) + "_" +
              eVarName + "_" + StringNumber(iTime, 4) + "_at_" + strFile;
          LINES_PLOT(FileName, eDrawArr, eCall, ePerm);
          if (PrintTextFiles) {
            WriteDrawArrFile(FileName, eDrawArr);
          }
        }
        if (nbPointOut > 0) {
          std::vector<MyMatrix<double>> ListData(nbGrid);
          for (int iGrid = 0; iGrid < nbGrid; iGrid++)
            ListData[iGrid] = TransectInterpolation_3D(ListTrans3D[iGrid],
                                                       ListRecVar[iGrid].Tens3);
          for (int iPt = 0; iPt < nbPointOut; iPt++) {
            double eLon = ListPointOut[iPt].lon;
            double eLat = ListPointOut[iPt].lat;
            std::string strLon = ConvertLLDoubleToDegMinSec(eLon);
            std::string strLat = ConvertLLDoubleToDegMinSec(eLat);
            std::cerr << "iTime=" << iTime << "/" << nbTime << " iPt=" << iPt
                      << "/" << nbPointOut << " name=" << ListPointOut[iPt].name
                      << " lon=" << ListPointOut[iPt].lon
                      << " lat=" << ListPointOut[iPt].lat << "\n";
            std::cerr << "    strLon=" << strLon << " strLat=" << strLat
                      << "\n";
            std::vector<MyVector<double>> ListListVect(nbGrid);
            std::vector<double> ListMax(nbGrid);
            std::vector<double> ListMin(nbGrid);
            int nbVert = ListNbVert[iPt];
            std::cerr << "nbVert=" << nbVert << "\n";
            for (int iGrid = 0; iGrid < nbGrid; iGrid++) {
              std::cerr << "  iGrid=" << iGrid << " / " << nbGrid << "\n";
              int nbData = ListData[iGrid].rows();
              MyVector<double> eVect(nbVert);
              int eFirstWet = ListTrans3D[iGrid].FirstWetIndex[iPt];
              std::cerr << "Before eVect assignment eFirstWet=" << eFirstWet
                        << "\n";
              std::cerr << "ListData[iGrid](rows/cols)="
                        << ListData[iGrid].rows() << " / "
                        << ListData[iGrid].cols() << "\n";
              for (int i = 0; i < nbVert; i++) {
                int j = nbData - 1 - i;
                double eVal = ListData[iGrid](j, iPt);
                int eMSK = ListTrans3D[iGrid].MSK(j, iPt);
                std::cerr << "    i=" << i << " eVal=" << eVal
                          << " eMSK=" << eMSK << "\n";
                eVect(i) = eVal;
              }
              std::cerr << "After eVect assignment\n";
              ListListVect[iGrid] = eVect;
              std::cerr << "After ListListVect assignation\n";
              double eMax = 0, eMin = 0;
              if (nbVert > 0) {
                eMax = eVect.maxCoeff();
                eMin = eVect.minCoeff();
              }
              ListMax[iGrid] = eMax;
              std::cerr << "After ListMax assignation\n";
              ListMin[iGrid] = eMin;
              std::cerr << "After ListMin assignation\n";
            }
            std::cerr << "ListListVect / ListMax / ListMin built\n";
            MyVector<double> ListX(nbVert);
            int nbData0 = ListData[0].rows();
            for (int i = 0; i < nbVert; i++)
              ListX(i) = -ListTrans3D[0].ListVertPos(nbData0 - 1 - i);
            std::cerr << "ListX constructed\n";
            std::string strNb = " " + std::to_string(iPt + 1);
            //
            std::string name = ListPointOut[iPt].name;
            DrawLinesArr eDrawArr;
            SetDefaultDrawLinesArr(eDrawArr);
            eDrawArr.DoTitle = true;
            eDrawArr.TitleStr =
                "Trans. " + name + " of " + eVarName + " at " + strPres;
            std::string fVarName =
              "Vertical_transect_" + name + "_" + eVarName
              + "_" + StringNumber(iTime,4) + "_" + strFile;
            eDrawArr.VarName = fVarName;
            eDrawArr.ListName_plot = ListRunName;
            eDrawArr.YAxisString = ListRecVar[0].RecS.VarName2 + "(" +
                                   ListRecVar[0].RecS.Unit + ")";
            eDrawArr.ListX = ListX;
            std::cerr << "eDrawArr done\n";
            //
            double TheMax, TheMin;
            double vectMax = VectorMax(ListMax), vectMin = VectorMin(ListMax);
            std::cerr << "vectMax=" << vectMax << " vectMin=" << vectMin
                      << "\n";
            if (!VariableRange) {
              TheMax = ListRecVar[0].RecS.maxval;
              TheMin = ListRecVar[0].RecS.minval;
            } else {
              TheMax = VectorMax(ListMax);
              TheMin = VectorMin(ListMin);
            }
            std::cerr << "Before rounding TheMin=" << TheMin << " TheMax=" << TheMax << "\n";
            TheMax = ApplyRounding(TheMax, VariableRangeRounding);
            TheMin = ApplyRounding(TheMin, VariableRangeRounding);
            std::cerr << "Final TheMin=" << TheMin << " TheMax=" << TheMax << "\n";
            eDrawArr.TheMax = TheMax;
            eDrawArr.TheMin = TheMin;
            //
            if (nbVert > 0 && TheMax < 10000) {
              eDrawArr.ListListVect = ListListVect;
              std::string FileName = ePerm.eDir + "VertTransect_" + name + "_" +
                                     eVarName + "_" + StringNumber(iTime, 4) +
                                     "_at_" + strFile;
              LINES_PLOT(FileName, eDrawArr, eCall, ePerm);
              if (PrintTextFiles) {
                WriteDrawArrFile(FileName, eDrawArr);
              }
            }
          }
        }
      }
    }
  }
}

// clang-format off
#endif  // SRC_OCEAN_TRANSECT_H_
// clang-format on
