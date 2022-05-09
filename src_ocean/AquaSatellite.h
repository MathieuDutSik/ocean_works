// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_OCEAN_AQUASATELLITE_H_
#define SRC_OCEAN_AQUASATELLITE_H_

#include "Basic_plot.h"
#include "Data_Access.h"
#include "NCL_Kernel.h"
#include "Timings.h"
#include <map>
#include <string>
#include <utility>
#include <vector>

FullNamelist NAMELIST_GetStandardAQUA() {
  std::map<std::string, SingleBlock> ListBlock;
  //  PROC
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::vector<double>> ListListDoubleValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  ListStringValues1["AquaLink"] =
      "https://oceandata.sci.gsfc.nasa.gov/MODIS-Aqua/L2/";
  ListStringValues1["PicPrefix"] = "/irrelevant";
  ListStringValues1["StorageFile"] = "unset StorageFile";
  ListBoolValues1["SendingEmail"] = false;
  ListListStringValues1["ListEmailAddress"] = {};
  ListStringValues1["__NaturePlot"] = "SINGLE";
  ListBoolValues1["KeepNC_NCL"] = false;
  ListBoolValues1["InPlaceRun"] = false;
  ListBoolValues1["PrintDebugInfo"] = false;
  ListBoolValues1["OnlyCreateFiles"] = false;
  ListBoolValues1["FirstCleanDirectory"] = true;
  ListStringValues1["Pcolor_method"] = "ncl";
  ListStringValues1["Quiver_method"] = "ncl";
  ListStringValues1["Lines_method"] = "ncl";
  ListStringValues1["Scatter_method"] = "ncl";
  ListIntValues1["NPROC"] = 2;
  ListStringValues1["Extension"] = "png";
  ListBoolValues1["StoreDownloadedData"] = false;
  ListStringValues1["StorePrefix"] = "/unset";
  //  ListBoolValues1["UseFixedList"]=false;
  ListStringValues1["MethodOperation"] = "unset";
  ListStringValues1["GetfileLinkPrefix"] = "";
  ListListStringValues1["ListNakedFile"] = {};
  ListIntValues1["EarliestDay_shift"] = -4;
  ListStringValues1["BEGTC"] = "19900101.000000";
  ListStringValues1["ENDTC"] = "19900201.000000";
  SingleBlock BlockPROC;
  BlockPROC.ListIntValues = ListIntValues1;
  BlockPROC.ListBoolValues = ListBoolValues1;
  BlockPROC.ListDoubleValues = ListDoubleValues1;
  BlockPROC.ListListDoubleValues = ListListDoubleValues1;
  BlockPROC.ListStringValues = ListStringValues1;
  BlockPROC.ListListStringValues = ListListStringValues1;
  ListBlock["PROC"] = BlockPROC;
  //  AQUA
  std::map<std::string, int> ListIntValues2;
  std::map<std::string, bool> ListBoolValues2;
  std::map<std::string, double> ListDoubleValues2;
  std::map<std::string, std::vector<double>> ListListDoubleValues2;
  std::map<std::string, std::string> ListStringValues2;
  std::map<std::string, std::vector<std::string>> ListListStringValues2;
  ListListStringValues2["ListVarName"] = {};
  ListDoubleValues2["MinLon"] = 12;
  ListDoubleValues2["MaxLon"] = 20;
  ListDoubleValues2["MinLat"] = 39;
  ListDoubleValues2["MaxLat"] = 46;
  ListBoolValues2["PlotSatelliteCover"] = false;
  ListStringValues2["CoastlineResolution_Lines"] = "LowRes";
  ListIntValues2["GridSubsample_lines_row"] = 10;
  ListIntValues2["GridSubsample_lines_col"] = 10;
  ListBoolValues2["GlobalGrid_lines"] = true;
  ListBoolValues2["PlotFullGeographicExtent"] = false;
  ListBoolValues2["SetMaximumRangeDynamic"] = true;
  SingleBlock BlockAQUA;
  BlockAQUA.ListIntValues = ListIntValues2;
  BlockAQUA.ListBoolValues = ListBoolValues2;
  BlockAQUA.ListDoubleValues = ListDoubleValues2;
  BlockAQUA.ListListDoubleValues = ListListDoubleValues2;
  BlockAQUA.ListStringValues = ListStringValues2;
  BlockAQUA.ListListStringValues = ListListStringValues2;
  ListBlock["AQUA"] = BlockAQUA;
  //
  return {ListBlock, "undefined"};
}

void AquaDownloading(FullNamelist const &eFull) {
  SingleBlock BlPROC = eFull.ListBlock.at("PROC");
  SingleBlock BlAQUA = eFull.ListBlock.at("AQUA");
  double MinLon = BlAQUA.ListDoubleValues.at("MinLon");
  double MaxLon = BlAQUA.ListDoubleValues.at("MaxLon");
  double MinLat = BlAQUA.ListDoubleValues.at("MinLat");
  double MaxLat = BlAQUA.ListDoubleValues.at("MaxLat");
  bool PlotSatelliteCover = BlAQUA.ListBoolValues.at("PlotSatelliteCover");
  std::vector<std::string> ListVarName =
      BlAQUA.ListListStringValues.at("ListVarName");
  std::string CoastlineResolution_Lines =
      BlAQUA.ListStringValues.at("CoastlineResolution_Lines");
  int GridSubsample_lines_row =
      BlAQUA.ListIntValues.at("GridSubsample_lines_row");
  int GridSubsample_lines_col =
      BlAQUA.ListIntValues.at("GridSubsample_lines_col");
  bool GlobalGrid_lines = BlAQUA.ListBoolValues.at("GlobalGrid_lines");
  bool PlotFullGeographicExtent =
      BlAQUA.ListBoolValues.at("PlotFullGeographicExtent");
  bool SetMaximumRangeDynamic =
      BlAQUA.ListBoolValues.at("SetMaximumRangeDynamic");
  std::vector<std::string> ListEmailAddress =
      BlPROC.ListListStringValues.at("ListEmailAddress");
  bool SendingEmail_L = BlPROC.ListBoolValues.at("SendingEmail");
  //
  std::string rndString = random_string(20);
  std::string TheDir = "/tmp/AQUA_DOWN_" + rndString + "/";
  CreateDirectory(TheDir);
  //
  PermanentInfoDrawing ePerm = GET_PERMANENT_INFO(eFull);
  NCLcaller<GeneralType> eCall(ePerm.NPROC);
  //
  QuadArray eQuad{MinLon, MaxLon, MinLat, MaxLat};
  auto IsGeographicallyCorrect = [&](DATAaqua const &eDATA) -> bool {
    std::vector<double> LONcornerDom{MinLon, MinLon, MaxLon, MaxLon};
    std::vector<double> LATcornerDom{MinLat, MaxLat, MaxLat, MinLat};
    std::vector<double> LONcornerSat(4), LATcornerSat(4);
    int nbRow = eDATA.LON.rows();
    int nbCol = eDATA.LON.cols();
    for (int i = 0; i < 4; i++) {
      int iRow, iCol;
      if (i == 0) {
        iRow = 0;
        iCol = 0;
      }
      if (i == 1) {
        iRow = nbRow - 1;
        iCol = 0;
      }
      if (i == 2) {
        iRow = nbRow - 1;
        iCol = nbCol - 1;
      }
      if (i == 3) {
        iRow = 0;
        iCol = nbCol - 1;
      }
      LONcornerSat[i] = eDATA.LON(iRow, iCol);
      LATcornerSat[i] = eDATA.LAT(iRow, iCol);
    }
    return IsContainedInDomain(LONcornerDom, LATcornerDom, LONcornerSat,
                               LATcornerSat);
  };
  auto GetCoreString = [&](std::string const &DataFile) -> std::string {
    std::vector<std::string> LStrB = STRING_Split(DataFile, "/");
    std::string NakedFileName = LStrB[LStrB.size() - 1];
    std::vector<std::string> LStrC = STRING_Split(NakedFileName, ".");
    std::string CoreString = LStrC[0];
    return CoreString;
  };
  std::string StorageFile = BlPROC.ListStringValues.at("StorageFile");
  bool StoreDownloadedData = BlPROC.ListBoolValues.at("StoreDownloadedData");
  std::string StorePrefix = BlPROC.ListStringValues.at("StorePrefix");
  std::vector<std::string> ListStringLink = ReadFullFile(StorageFile);
  auto RealFuncInsertLink = [&](std::string const &eLink) -> void {
    std::ofstream os(StorageFile, std::ofstream::app);
    os << eLink << "\n";
    ListStringLink.push_back(eLink);
    std::cerr << "Now |ListStringLink|=" << ListStringLink.size() << "\n";
  };
  auto GetLastSlashEntry = [&](std::string const &eLink) -> std::string {
    std::vector<std::string> LStrB = STRING_Split(eLink, "/");
    std::string NakedFileName = LStrB[LStrB.size() - 1];
    //  std::string DataFile=TheDir + "/" + NakedFileName;
    return NakedFileName;
  };
  auto DownLink = [&](std::string const &eLink) -> std::string {
    std::vector<std::string> LStrB = STRING_Split(eLink, "/");
    std::string NakedFileName = LStrB[LStrB.size() - 1];
    std::string tmpStdout = TheDir + "tmp_" + NakedFileName + "_stdout";
    //  std::string tmpStderr=TheDir + "tmp_" + NakedFileName + "_stdout";
    unsigned int TimeSleep = 10;
    while (true) {
      std::string eComm =
          "(cd " + TheDir + " && wget " + eLink + " > " + tmpStdout + ")";
      int iret = system(eComm.c_str());
      if (iret != 0) {
        std::cerr << "iret=" << iret << "\n";
        std::cerr << "Failed in execution of command\n";
        std::cerr << "eComm=" << eComm << "\n";
        sleep(TimeSleep);
        TimeSleep *= 2;
      } else {
        break;
      }
    }
    std::string DataFile = TheDir + NakedFileName;
    return DataFile;
  };
  auto CreateStorageDirectory = [&](std::string const &DataFile) -> void {
    std::string CoreString = GetCoreString(DataFile);
    std::vector<int> LTime = ConvertAquaNameToTime(CoreString);
    int eYear = LTime[0];
    int eMonth = LTime[1];
    int eDay = LTime[2];
    std::string StorDir = StorePrefix + StringNumber(eYear, 4) + "/" +
                          StringNumber(eMonth, 2) + "/" +
                          StringNumber(eDay, 2) + "/";
    CreateDirectory(StorDir);
  };
  auto GetStorageFileName = [&](std::string const &DataFile) -> std::string {
    std::string CoreString = GetCoreString(DataFile);
    std::vector<int> LTime = ConvertAquaNameToTime(CoreString);
    int eYear = LTime[0];
    int eMonth = LTime[1];
    int eDay = LTime[2];
    return StorePrefix + StringNumber(eYear, 4) + "/" +
           StringNumber(eMonth, 2) + "/" + StringNumber(eDay, 2) + "/" +
           CoreString + ".nc";
  };
  auto GetNameUselessStorage = [&](std::string const &DataFile) -> std::string {
    std::string CoreString = GetCoreString(DataFile);
    std::vector<int> LTime = ConvertAquaNameToTime(CoreString);
    int eYear = LTime[0];
    int eMonth = LTime[1];
    int eDay = LTime[2];
    CreateStorageDirectory(DataFile);
    return StorePrefix + StringNumber(eYear, 4) + "/" +
           StringNumber(eMonth, 2) + "/" + StringNumber(eDay, 2) +
           "/ListUselessFiles";
  };
  auto GetListUselessFiles =
      [&](std::string const &DataFile) -> std::vector<std::string> {
    std::string NameStorUseless = GetNameUselessStorage(DataFile);
    if (IsExistingFile(NameStorUseless))
      return ReadFullFile(NameStorUseless);
    return {};
  };
  auto IsKnownUselessFile = [&](std::string const &DataFile) -> bool {
    std::string LastFileName = GetLastSlashEntry(DataFile);
    std::vector<std::string> ListUseless = GetListUselessFiles(LastFileName);
    if (PositionVect(ListUseless, LastFileName) == -1)
      return false;
    return true;
  };
  auto InsertInUselessFiles = [&](std::string const &DataFile) -> void {
    std::string LastFileName = GetLastSlashEntry(DataFile);
    std::string NameStorUseless = GetNameUselessStorage(DataFile);
    std::ofstream os(NameStorUseless, std::ofstream::app);
    os << LastFileName << "\n";
  };
  auto SetFile = [&](std::string const &eLink) -> std::pair<bool, std::string> {
    if (StoreDownloadedData) {
      std::string DataFile = TheDir + GetLastSlashEntry(eLink);
      std::string StorageFile = GetStorageFileName(DataFile);
      if (IsExistingFile(StorageFile)) {
        CopyOperation(StorageFile, DataFile);
        return {true, DataFile};
      } else {
        if (IsKnownUselessFile(DataFile)) {
          return {false, DataFile};
        } else {
          return {true, DownLink(eLink)};
        }
      }
    } else {
      return {true, DownLink(eLink)};
    }
  };
  auto StorageOperation = [&](std::string const &DataFile,
                              bool const &testIncl) -> void {
    if (testIncl) {
      std::string StorageFile = GetStorageFileName(DataFile);
      if (!IsExistingFile(StorageFile)) {
        CreateStorageDirectory(DataFile);
        CopyOperation(DataFile, StorageFile);
      }
    } else {
      InsertInUselessFiles(DataFile);
    }
  };
  //
  //  The drawing operations (main stuff actually)
  //
  auto FullyTreatFile = [&](std::string const &DataFile) -> void {
    std::cerr << "-----------------------------------------------------\n";
    std::string CoreString = GetCoreString(DataFile);
    std::cerr << "CoreString=" << CoreString << "\n";

    DATAaqua eData = ReadAquaSatellite(DataFile, "chlor_a");
    if (eData.LON.minCoeff() < -200 || eData.LON.maxCoeff() > 200 ||
        eData.LAT.minCoeff() < -100 || eData.LAT.maxCoeff() > 100) {
      std::cerr << "DataFile=" << DataFile << "\n";
      std::cerr << "|eData.LON| = " << StringSizeMatrix(eData.LON) << " "
                << MinMaxMatrix(eData.LON) << "\n";
      std::cerr << "|eData.LAT| = " << StringSizeMatrix(eData.LAT) << " "
                << MinMaxMatrix(eData.LAT) << "\n";
      std::cerr << "eData(LON/LAT) is not in the right size\n";
      return;
    }
    //
    double midTime = (eData.minTime + eData.maxTime) / 2;
    std::string strTimePresFirst = DATE_ConvertMjd2mystringPres(eData.minTime);
    std::string strTimePresLast = DATE_ConvertMjd2mystringPres(eData.maxTime);
    std::string strTimePresMid = DATE_ConvertMjd2mystringPres(midTime);
    std::string strTimeFileFirst = DATE_ConvertMjd2mystringFile(eData.minTime);
    std::string strTimeFileLast = DATE_ConvertMjd2mystringFile(eData.maxTime);
    std::string strTimeFileMid = DATE_ConvertMjd2mystringFile(midTime);
    GridArray GrdArr = CURVILINEAR_GRID_ARRAY(eData.LON, eData.LAT);
    MyMatrix<double> LONred = MatrixSubsample(
        eData.LON, GridSubsample_lines_row, GridSubsample_lines_col);
    MyMatrix<double> LATred = MatrixSubsample(
        eData.LAT, GridSubsample_lines_row, GridSubsample_lines_col);
    GridArray GrdArrRed = CURVILINEAR_GRID_ARRAY(LONred, LATred);
    if (PlotSatelliteCover) {
      std::string TitleStr = "Satellite cover around " + strTimePresMid;
      std::string FileName = ePerm.eDir + "Lines_" + CoreString + "_-_" +
                             strTimeFileFirst + "_" + strTimeFileLast;
      std::cerr << "FileName = " << FileName << "\n";
      //
      std::vector<PairLL> ListPairLL = GetGridBoundaryTot(eData.LON, eData.LAT);
      SeqLineSegment eSeq{ListPairLL, true};
      std::vector<SeqLineSegment> ListLineSegment{eSeq};
      //
      RecVar eRecVarRed = GetTrivialArrayPlot(GrdArrRed);
      //
      QuadArray eQuadGA;
      if (GlobalGrid_lines) {
        eQuadGA = {static_cast<double>(-180), static_cast<double>(180),
                   static_cast<double>(-85), static_cast<double>(85)};
      } else {
        eQuadGA = GetQuadArray(GrdArr);
      }
      DrawArr eDrw = BasicArrayDraw(eQuadGA);
      eDrw.DoTitle = true;
      eDrw.TitleStr = TitleStr;
      eDrw.ListLineSegment = ListLineSegment;
      eDrw.GridResolution = CoastlineResolution_Lines;
      eDrw.VarNameUF = "Lines_" + CoreString;
      //
      std::cerr << "  Before PLOT_PCOLOR lines FileName=" << FileName << "\n";
      PLOT_PCOLOR(FileName, GrdArrRed, eDrw, eRecVarRed, eCall, ePerm);
      std::cerr << "  After PLOT_PCOLOR lines\n";
    }
    TotalArrGetData TotalArr;
    TotalArr.GrdArr.ModelName = "TRIVIAL";
    bool testIncl = IsGeographicallyCorrect(eData);
    StorageOperation(DataFile, testIncl);
    std::cerr << "testIncl=" << testIncl << "\n";
    std::vector<std::string> ListPictureFile;
    if (testIncl) {
      for (auto &eVarName : ListVarName) {
        std::cerr << "eVarName = " << eVarName << "\n";
        RecVar eRecVar =
            ModelSpecificVarSpecificTime(TotalArr, eVarName, midTime);
        eRecVar.RecS.strAll = "aquadownloading_modelspecificvar";
        std::string CFname = eRecVar.RecS.CFshortName;
        DATAaqua fData = ReadAquaSatellite(DataFile, CFname);
        std::cerr << "  After ReadAquaSatellite\n";
        std::string eVarName2 = eRecVar.RecS.VarName2;
        //
        std::string TitleStr = eVarName2 + " around " + strTimePresMid;
        std::string FileName = ePerm.eDir + "Pcolor_" + eVarName + "_" +
                               strTimeFileFirst + "_" + strTimeFileLast;
        std::cerr << "FileName = " << FileName << "\n";
        //
        QuadArray eQuadPC;
        if (PlotFullGeographicExtent)
          eQuadPC = GetQuadArray(GrdArr);
        else
          eQuadPC = eQuad;
        DrawArr eDrw = BasicArrayDraw(eQuadPC);
        eDrw.DoTitle = true;
        eDrw.TitleStr = TitleStr;
        eDrw.GridResolution = "MediumRes";
        eDrw.VarNameUF = "PCOLOR_" + CoreString + "_" + eVarName;
        eRecVar.F = fData.VAR;
        int nbVal = 1;
        if (SetMaximumRangeDynamic) {
          int nbRow = fData.MSK.rows();
          int nbCol = fData.MSK.cols();
          std::vector<double> ListVal;
          for (int iRow = 0; iRow < nbRow; iRow++)
            for (int iCol = 0; iCol < nbCol; iCol++) {
              if (fData.MSK(iRow, iCol) == 1) {
                double eLon = fData.LON(iRow, iCol);
                double eLat = fData.LON(iRow, iCol);
                if (IsPointInQuadArray(eQuadPC, eLon, eLat))
                  ListVal.push_back(fData.VAR(iRow, iCol));
              }
            }
          nbVal = ListVal.size();
          std::cerr << "nbRow=" << nbRow << " nbCol=" << nbCol
                    << " nbVal=" << nbVal << "\n";
          if (nbVal > 0) {
            double TheMax = VectorMax(ListVal);
            double TheMin = VectorMin(ListVal);
            std::cerr << "TheMin=" << TheMin << " TheMax=" << TheMax << "\n";
            eRecVar.RecS.minval = TheMin;
            eRecVar.RecS.maxval = TheMax;
          }
        }
        double deltaMM = eRecVar.RecS.maxval - eRecVar.RecS.minval;
        std::cerr << "STRALL 2 : eVarName=" << eVarName
                  << " eRecVar.strAll=" << eRecVar.RecS.strAll << "\n";
        std::cerr << "deltaMM=" << deltaMM << " nbVal=" << nbVal << "\n";
        //
        if (deltaMM > 0 && nbVal > 0) {
          std::cerr << "  Before PLOT_PCOLOR FileName=" << FileName << "\n";
          PLOT_PCOLOR(FileName, GrdArr, eDrw, eRecVar, eCall, ePerm);
          std::cerr << "  After PLOT_PCOLOR\n";
          ListPictureFile.push_back(FileName + "." + ePerm.Extension);
        }
      }
      if (SendingEmail_L && ListPictureFile.size()) {
        std::string Subject = "AQUA Biology at " + strTimePresMid;
        std::string Text = "Dear receivers,\n";
        Text += "\n";
        Text += "please find in attachment picture downloaded from the Aqua "
                "satellite\n";
        Text += "\n";
        Text += "  Mathieu\n";
        std::string FromEmail = "mathieu.dutour@gmail.com";
        for (auto &DestEmail : ListEmailAddress) {
          SendMailOper eSendOper{DestEmail, FromEmail, Subject, Text,
                                 ListPictureFile};
          GeneralType eGen(eSendOper);
          std::cerr << "Sending Email with DestEmail = " << DestEmail << "\n";
          eCall.SubmitJob(eGen);
        }
      }
    }
    std::cerr << "-----------------------------------------------------\n";
  };
  //
  //  The functions for downloading
  //
  auto DownLinkAndTreat = [&](std::string const &eLink) -> void {
    std::pair<bool, std::string> ePair = SetFile(eLink);
    if (ePair.first) {
      std::string DataFile = ePair.second;
      FullyTreatFile(DataFile);
      RemoveFileIfExist(DataFile);
    }
  };
  auto DownloadLinkAndTreat_if_OC = [&](std::string const &eLink) -> void {
    std::vector<std::string> LStr = STRING_Split(eLink, "_");
    std::string strLast = LStr[LStr.size() - 1];
    if (strLast == "OC.nc")
      DownLinkAndTreat(eLink);
  };
  auto DownloadLinkAndTreat_if_OC_and_undone =
      [&](std::string const &eLink) -> void {
    if (PositionVect(ListStringLink, eLink) == -1) {
      RealFuncInsertLink(eLink);
      std::cerr << "Treating eLink=" << eLink << "\n";
      DownloadLinkAndTreat_if_OC(eLink);
    } else {
      std::cerr << "eLink = " << eLink << " already treated\n";
    }
  };
  //
  //  The different possibilities for downloading.
  //
  std::string MethodOperation = BlPROC.ListStringValues.at("MethodOperation");
  bool DidSomething = false;
  if (MethodOperation == "UseFixedList") {
    DidSomething = true;
    std::vector<std::string> ListNakedFile =
        BlPROC.ListListStringValues.at("ListNakedFile");
    std::string GetfileLinkPrefix =
        BlPROC.ListStringValues.at("GetfileLinkPrefix");
    std::string FirstChar = GetfileLinkPrefix.substr(0, 1);
    if (FirstChar != "/" && FirstChar != "h") {
      std::cerr << "First character should be / or h, i.e. full path to file "
                   "or https\n";
      throw TerminalException{1};
    }
    if (FirstChar == "h") {
      for (auto &eLinkShort : ListNakedFile) {
        std::string eFullLink = GetfileLinkPrefix + eLinkShort;
        DownLinkAndTreat(eFullLink);
      }
    }
    if (FirstChar == "/") {
      for (auto &eLinkShort : ListNakedFile) {
        std::string eFullLink = GetfileLinkPrefix + eLinkShort;
        FullyTreatFile(eFullLink);
      }
    }
  }
  if (MethodOperation == "ContinuousOperation") {
    DidSomething = true;
    while (1) {
      std::vector<std::string> ListLink = GetAquaListLinkDownload(eFull);
      for (auto &eLink : ListLink)
        DownloadLinkAndTreat_if_OC_and_undone(eLink);
      sleep(600);
    }
  }
  if (MethodOperation == "FixedPeriod") {
    DidSomething = true;
    std::string strBEGTC = BlPROC.ListStringValues.at("BEGTC");
    std::string strENDTC = BlPROC.ListStringValues.at("ENDTC");
    double MjdFirst = CT2MJD(strBEGTC);
    double MjdLast = CT2MJD(strENDTC);
    std::vector<std::string> ListLink =
        GetAquaListLinkDownload_Specified(eFull, MjdFirst, MjdLast);
    for (auto &eLink : ListLink)
      DownloadLinkAndTreat_if_OC_and_undone(eLink);
  }
  if (!DidSomething) {
    std::cerr << "We did not select a correct method of operation\n";
    std::cerr << "Allowed values for MethodOperation: FixedPeriod, "
                 "ContinuousOperation, UseFixedList\n";
    throw TerminalException{1};
  }
}

// clang-format off
#endif  // SRC_OCEAN_AQUASATELLITE_H_
// clang-format on
