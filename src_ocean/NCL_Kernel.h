// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_OCEAN_NCL_KERNEL_H_
#define SRC_OCEAN_NCL_KERNEL_H_

#include "BasicThreadInclude.h"
#include "Basic_file.h"
#include "EmailBashSystem.h"
#include "Model_interpolation.h"
#include "NamelistExampleOcean.h"
#include <vector>
#include <map>
#include <string>
#include <utility>

std::string NCL_bool(bool eBool) {
  if (eBool)
    return "True";
  return "False";
}

struct AnnotationRec {
  bool DrawAnnotation = false;
  double AnnotationLon;
  double AnnotationLat;
  std::string AnnotationText;
};

struct SeqLineSegment {
  std::vector<PairLL> ListPairLL;
  bool IsClosed;
};

double AveragePairLength(SeqLineSegment const &eSeq) {
  double eSum = 0;
  int nb = eSeq.ListPairLL.size();
  if (!eSeq.IsClosed)
    nb--;
  for (int i = 0; i < nb; i++) {
    int j = (i + 1) % nb;
    double distKM =
        GeodesicDistanceKM_pair(eSeq.ListPairLL[i], eSeq.ListPairLL[j]);
    eSum += distKM;
  }
  double retVal = eSum / double(nb);
  return retVal;
}

struct SingleMarker {
  PairLL Coord;
  int color;
  double thickness;
};

struct DrawArr {
  bool DoTitleString;
  bool DoTitle;
  std::string TitleStr;
  std::string VarNameUF;
  bool DrawRiver;
  bool DrawContourBathy;
  bool PrintMMA;
  std::string LandPortr;
  std::string optStatStr;
  // colormaps
  bool DoColorBar;
  std::string ColorMap;
  std::string cnFillMode;
  bool cnFillOn;
  bool cnLinesOn;
  bool cnLineLabelsOn;
  bool cnSmoothingOn;
  int nbLevelSpa;
  int nbLabelStride;
  // Frame
  QuadArray eQuadFrame;
  // annotations
  AnnotationRec TheAnnot;
  std::vector<std::string> ListInsertLines;
  double vcRefLengthF;
  bool FillLand;
  bool UseNativeGrid;
  std::string GridResolution;
  std::vector<SeqLineSegment> ListLineSegment;
  std::vector<SingleMarker> ListMarker;
};

struct QuadDrawInfo {
  std::string eFrameName;
  int iFrame;
  QuadArray eQuad;
};

struct InterpolToUVpoints {
  GridArray GrdArr;
  SingleArrayInterpolation InterpArr;
};

struct BChoices {
  bool KeepNC_NCL;
  bool InPlaceRun;
  bool PrintDebugInfo;
  bool OnlyCreateFiles;
};

struct TripleNCL {
  std::string eProg;
  std::string TargetFile;
  std::string eFileNC;
  std::string eFileNCL;
  BChoices eBChoice;
};

void ApplyPdfcrop(std::string const &PdfFile) {
  std::string TmpFile = "TMP_" + random_string(20) + ".pdf";
  std::string TheDir = FILE_GetDirectoryOfFileName(PdfFile);
  //  std::cerr << "TheDir=" << TheDir << "\n";
  std::string FileNaked = FILE_GetNakedFilename(PdfFile);
  //  std::cerr << "FileNaked=" << FileNaked << "\n";
  std::string eCommand1 =
      "(cd " + TheDir + " && pdfcrop " + FileNaked + " " + TmpFile + ")";
  std::cerr << "eCommand1=" << eCommand1 << "\n";
  int iret1 = system(eCommand1.c_str());
  //  std::cerr << "iret1=" << iret1 << "\n";
  if (iret1 != 0) {
    std::cerr << "Error at pdfcrop operation\n";
    throw TerminalException{1};
  }
  //
  std::string eCommand2 = "rm " + PdfFile;
  std::cerr << "eCommand2=" << eCommand2 << "\n";
  int iret2 = system(eCommand2.c_str());
  //  std::cerr << "iret2=" << iret2 << "\n";
  if (iret2 != 0) {
    std::cerr << "Error at rm operation\n";
    throw TerminalException{1};
  }
  //
  std::string TmpFileFull = TheDir + TmpFile;
  //  std::cerr << "TmpFileFull=" << TmpFileFull << "\n";
  std::string eCommand3 = "mv " + TmpFileFull + " " + PdfFile;
  std::cerr << "eCommand3=" << eCommand3 << "\n";
  int iret3 = system(eCommand3.c_str());
  //  std::cerr << "iret3=" << iret3 << "\n";
  if (iret3 != 0) {
    std::cerr << "Error at mv operation\n";
    throw TerminalException{1};
  }
}

void ApplyPictureTrim(std::string const &PngFile) {
  std::string eCommand1 = "convert -trim " + PngFile + " " + PngFile;
  std::cerr << "eCommand1=" << eCommand1 << "\n";
  int iret1 = system(eCommand1.c_str());
  if (iret1 == 1 || iret1 == 256) {
    std::cerr << "iret1=" << iret1 << "\n";
    std::cerr << "The operation convert -trim failed\n";
    throw TerminalException{1};
  }
}

std::string GetNakedComm(std::string const &eCommIn) {
  if (eCommIn == "ncl_crop")
    return "ncl";
  return eCommIn;
}

void CroppingOperation(std::string const &eTargetFile) {
  std::string eExtension = FILE_GetExtension(eTargetFile);
  std::cerr << "CroppingOperation eExtension=" << eExtension << "\n";
  if (eExtension == "png") {
    ApplyPictureTrim(eTargetFile);
    return;
  }
  if (eExtension == "pdf") {
    ApplyPdfcrop(eTargetFile);
    return;
  }
  std::cerr << "Right now only png and pdf are supported for cropping\n";
  std::cerr << "eTargetFile=" << eTargetFile << " and eExtension=" << eExtension
            << "\n";
  throw TerminalException{1};
}

void CALL_PICTURE_PROG(TripleNCL const &eTripl) {
  if (eTripl.eBChoice.OnlyCreateFiles) {
    if (eTripl.eBChoice.PrintDebugInfo) {
      std::cerr << "-----------------------------------------------------------"
                   "--------\n";
      std::cerr << "   Begin of CALL_PICTURE_PROG\n";
      std::cerr
          << "   Only creating NC, PY, NCL files, so no operations needed\n";
    }
    return;
  }
  //
  std::string eComm = GetNakedComm(eTripl.eProg);
  if (eTripl.eBChoice.PrintDebugInfo) {
    std::cerr << "-------------------------------------------------------------"
                 "------\n";
    std::cerr << "   Begin of CALL_PICTURE_PROG\n";
    std::cerr << "   eTripl.eProg=" << eTripl.eProg << " eComm=" << eComm
              << "\n";
  }
  auto GetFinalTargetFile = [](std::string const &eFileIn) -> std::string {
    std::vector<std::string> LVect = STRING_Split(eFileIn, "_storsave");
    if (LVect.size() == 1)
      return eFileIn;
    if (LVect.size() == 2)
      return LVect[0] + LVect[1];
    std::cerr << "Error in the FinalTargetFile\n";
    throw TerminalException{1};
  };
  std::string FinalTargetFile = GetFinalTargetFile(eTripl.TargetFile);
  if (eTripl.eBChoice.PrintDebugInfo)
    std::cerr << "eTripl.TargetFile=" << eTripl.TargetFile
              << " FinalTargetFile=" << FinalTargetFile << "\n";
  auto GetCommand = [&]() -> std::string {
    std::string devNull = " > /dev/null 2> /dev/null";
    if (FinalTargetFile == eTripl.TargetFile) {
      if (eTripl.eBChoice.InPlaceRun) {
        std::string TheDir = FILE_GetDirectoryOfFileName(eTripl.TargetFile);
        return "(cd " + TheDir + " && " + eComm + " " +
               FILE_GetNakedFilename(eTripl.eFileNCL) + devNull + ")";
      }
      return eComm + " " + eTripl.eFileNCL + devNull;
    } else {
      if (eTripl.eBChoice.InPlaceRun) {
        std::string TheDir = FILE_GetDirectoryOfFileName(eTripl.TargetFile);
        return "(cd " + TheDir + " && " + eComm + " " +
               FILE_GetNakedFilename(eTripl.eFileNCL) + devNull + " && mv " +
               eTripl.TargetFile + " " + FinalTargetFile + ")";
      }
      return "(" + eComm + " " + eTripl.eFileNCL + devNull + " && mv " +
             eTripl.TargetFile + " " + FinalTargetFile + ")";
    }
  };
  std::string eCommand = GetCommand();
  if (eTripl.eBChoice.PrintDebugInfo) {
    std::cerr << "eCommand = " << eCommand << "\n";
  }
  int iret = system(eCommand.c_str());
  //  std::cerr << "iret=" << iret << "\n";
  if (iret == -1) {
    printf("Oh dear, something went wrong with ncl! %s\n", strerror(errno));
    throw TerminalException{1};
  }
  if (eTripl.eBChoice.PrintDebugInfo) {
    std::cerr << "eFileNC  = " << eTripl.eFileNC << "\n";
    std::cerr << "eFileNCL = " << eTripl.eFileNCL << "\n";
  }
  if (!IsExistingFile(FinalTargetFile)) {
    std::cerr << "The following TargetFile was not created\n";
    std::cerr << "FinalTargetFile = " << FinalTargetFile << "\n";
    std::cerr << "TargetFile = " << eTripl.TargetFile << "\n";
    std::cerr << "eFileNC    = " << eTripl.eFileNC << "\n";
    std::cerr << "eFileNCL   = " << eTripl.eFileNCL << "\n";
    std::cerr << "Please debug\n";
    throw TerminalException{1};
  }
  if (eTripl.eProg == "ncl_crop")
    CroppingOperation(FinalTargetFile);
  if (!eTripl.eBChoice.KeepNC_NCL) {
    RemoveFileIfExist(eTripl.eFileNC);
    RemoveFileIfExist(eTripl.eFileNCL);
  }
}

struct GeneralType {
  GeneralType() { typeoper = "unset"; }
  GeneralType(TripleNCL const &fTripl) {
    typeoper = "picture";
    eTripl = fTripl;
  }
  GeneralType(SendMailOper const &fSendOper) {
    typeoper = "sendattachment";
    eSendOper = fSendOper;
  }
  GeneralType(BashOper const &fBashOper) {
    typeoper = "bashoperation";
    eBashOper = fBashOper;
  }
  std::string typeoper;
  TripleNCL eTripl;
  SendMailOper eSendOper;
  BashOper eBashOper;
};

void OPERATION_PROG(GeneralType const &eGen) {
  if (eGen.typeoper == "picture")
    return CALL_PICTURE_PROG(eGen.eTripl);
  if (eGen.typeoper == "sendattachment")
    return CALL_SEND_ATTACHMENT(eGen.eSendOper);
  if (eGen.typeoper == "bashoperation")
    return CALL_BASH_OPERATION(eGen.eBashOper);
  std::cerr << "Failed to find matching operation\n";
  std::cerr << "typeoper=" << eGen.typeoper << "\n";
  throw TerminalException{1};
}

// My attempt at a thread-pool.
// One problem of code below is that only one thread can submit jobs to the
// bank of threads.
template <typename Toperation> struct NCLcaller {
  NCLcaller() = delete;

  NCLcaller(int const &nproc) : NPROC(nproc) {
    ListExch = new int[NPROC];
    ListTerm = new int[NPROC];
    ListOperation = new Toperation[NPROC];
    ListCond = std::vector<std::condition_variable>(nproc);
    NbRunningJob = 0;
    InWhile = 0;
    for (int iProc = 0; iProc < NPROC; iProc++) {
      ListExch[iProc] = 0;
      ListTerm[iProc] = 0;
    }
    //
    auto IterationLoop = [&](int iproc, int *Exch, int *Term) {
      int IsFirst = 1;
      while (true) {
        if (IsFirst == 1) {
          InWhile++;
          IsFirst = 0;
        }
        if (*Exch == 1) {
          OPERATION_PROG(ListOperation[iproc]);
          *Exch = 0;
          NbRunningJob--;
          sub_cv.notify_one();
        }
        if (*Term == -1) {
          break;
        }
        if (*Exch == 0) {
          std::unique_lock<std::mutex> lk(inst_mtx);
          ListCond[iproc].wait(lk, [&] { return *Exch == 1 || *Term == -1; });
        }
      }
      InWhile--;
      fin_cv.notify_one();
    };
    for (int iProc = 0; iProc < NPROC; iProc++)
      ListThr.push_back(std::thread(IterationLoop, iProc, &(ListExch[iProc]),
                                    &(ListTerm[iProc])));
    for (int iProc = 0; iProc < NPROC; iProc++)
      ListThr[iProc].detach();
  }
  ~NCLcaller() {
    for (int iProc = 0; iProc < NPROC; iProc++) {
      ListTerm[iProc] = -1;
      ListCond[iProc].notify_one();
    }
    std::unique_lock<std::mutex> lk(fin_mtx);
    fin_cv.wait(lk, [&] { return InWhile == 0; });
    delete[] ListExch;
    delete[] ListTerm;
    delete[] ListOperation;
  }
  //
  void SubmitJob(Toperation const &eJob) {
    std::unique_lock<std::mutex> lk(sub_mtx);
    sub_cv.wait(lk, [&] { return NbRunningJob < NPROC; });
    int iProcFound = -1;
    for (int iProc = 0; iProc < NPROC; iProc++)
      if (ListExch[iProc] == 0)
        iProcFound = iProc;
    if (iProcFound == -1) {
      std::cerr << "Failed to find the processor.\n";
      std::cerr << "Bug is in NCLcaller\n";
      std::cerr << "NbRunningJob=" << NbRunningJob << " NPROC=" << NPROC
                << "\n";
      for (int iProc = 0; iProc < NPROC; iProc++)
        std::cerr << "iProc=" << iProc << " Exch=" << ListExch[iProc]
                  << " Term=" << ListTerm[iProc] << "\n";
      throw TerminalException{1};
    }
    ListOperation[iProcFound] = eJob;
    ListExch[iProcFound] = 1;
    NbRunningJob++;
    ListCond[iProcFound].notify_one();
  }

private:
  std::mutex fin_mtx;
  std::condition_variable fin_cv;
  //
  std::atomic<int> NbRunningJob;
  std::atomic<int> InWhile;
  //
  std::mutex sub_mtx;
  std::condition_variable sub_cv;
  //
  std::mutex inst_mtx;
  int NPROC;
  std::vector<std::thread> ListThr;
  int *ListExch;
  int *ListTerm;
  Toperation *ListOperation;
  std::vector<std::condition_variable> ListCond;
};

struct PermanentInfoDrawing {
  TempDirectory PrefixTemp;
  std::string eDir;
  std::string Extension;
  std::string PicPrefix;
  FullNamelist eFull;
  BChoices eBChoice;
  int NPROC;
  //
  // The basic setting for the drawing of data
  //
  DrawArr eDrawArr;
  std::map<std::string, std::string> ChoiceProgram;
  //
  // For finite element plot, we need data for the current
  //
  std::vector<QuadDrawInfo> ListQuadInfo;
  std::vector<InterpolToUVpoints> ListInterpol;
  //
  // For transect plot of 3D variables
  //
  std::vector<TransectInformation_3D> ListTransect;
};

TempDirectory PLOT_CreatePrefixTemp(FullNamelist const &eFull) {
  std::map<std::string, SingleBlock> ListBlock = eFull.ListBlock;
  SingleBlock eBlPROC = ListBlock.at("PROC");
  bool KeepNC_NCL = eBlPROC.ListBoolValues.at("KeepNC_NCL");
  std::string eRand = random_string(20);
  std::string Nature = "unset_Nature";
  auto iter = eBlPROC.ListStringValues.find("__NaturePlot");
  if (iter != eBlPROC.ListStringValues.end())
    Nature = eBlPROC.ListStringValues.at("__NaturePlot");
  std::string PrefixTemp = "/tmp/PLOT_" + Nature + "_" + eRand + "/";
  if (KeepNC_NCL)
    std::cerr << "PrefixTemp = " << PrefixTemp << "\n";
  return TempDirectory(PrefixTemp);
}

PermanentInfoDrawing GET_PERMANENT_INFO(FullNamelist const &eFull) {
  SingleBlock eBlPROC = eFull.ListBlock.at("PROC");
  std::string PicPrefix = eBlPROC.ListStringValues.at("PicPrefix");
  std::string Extension = eBlPROC.ListStringValues.at("Extension");
  bool KeepNC_NCL = eBlPROC.ListBoolValues.at("KeepNC_NCL");
  bool InPlaceRun = eBlPROC.ListBoolValues.at("InPlaceRun");
  bool PrintDebugInfo = eBlPROC.ListBoolValues.at("PrintDebugInfo");
  bool OnlyCreateFiles = eBlPROC.ListBoolValues.at("OnlyCreateFiles");
  BChoices eBChoice{KeepNC_NCL, InPlaceRun, PrintDebugInfo, OnlyCreateFiles};
  int NPROC = eBlPROC.ListIntValues.at("NPROC");
  TempDirectory PrefixTemp = PLOT_CreatePrefixTemp(eFull);
  std::string Pcolor_method = eBlPROC.ListStringValues.at("Pcolor_method");
  std::string Quiver_method = eBlPROC.ListStringValues.at("Quiver_method");
  std::string Lines_method = eBlPROC.ListStringValues.at("Lines_method");
  std::string Scatter_method = eBlPROC.ListStringValues.at("Scatter_method");
  std::map<std::string, std::string> ChoiceProgram;
  ChoiceProgram["Pcolor_method"] = Pcolor_method;
  ChoiceProgram["Quiver_method"] = Quiver_method;
  ChoiceProgram["Lines_method"] = Lines_method;
  ChoiceProgram["Scatter_method"] = Scatter_method;
  if (eBlPROC.ListBoolValues.at("FirstCleanDirectory"))
    RemoveFileInDirectory(PicPrefix);
  std::string eDir = FILE_GetAbsoluteDirectory(PicPrefix);
  CreateDirectory(eDir);
  //
  PermanentInfoDrawing ePerm;
  ePerm.PrefixTemp = std::move(PrefixTemp);
  ePerm.eDir = eDir;
  ePerm.Extension = Extension;
  ePerm.ChoiceProgram = ChoiceProgram;
  ePerm.PicPrefix = PicPrefix;
  ePerm.eFull = eFull;
  ePerm.eBChoice = eBChoice;
  ePerm.NPROC = NPROC;
  //  A priori it should work all ok for keeping files since
  //  they are all different (this is checked)
  //  if (KeepNC_NCL && NPROC > 1) {
  //    std::cerr << "Cannot have KeepNC_NCL = T and NPROC > 1\n";
  //    throw TerminalException{1};
  //  }
  return ePerm;
}

// clang-format off
#endif  // SRC_OCEAN_NCL_KERNEL_H_
// clang-format on
