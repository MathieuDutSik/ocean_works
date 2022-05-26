// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_OCEAN_BASIC_PLOT_H_
#define SRC_OCEAN_BASIC_PLOT_H_

#include "Basic_Ocean_types.h"
#include "Basic_netcdf.h"
#include "MAT_Matrix.h"
#include "NCL_Kernel.h"
#include "SVGfunctions.h"
#include "Statistics.h"
#include "Triangulations.h"
#include "mjdv2.h"
#include <algorithm>
#include <string>
#include <vector>

std::string FinalFile(bool const &InPlaceRun, std::string const &TargetFile,
                      std::string const &eFile) {
  if (InPlaceRun) {
    std::string eDir = FILE_GetDirectoryOfFileName(TargetFile);
    std::string RetFile = eDir + FILE_GetNakedFilename(eFile);
    //    std::cerr << "FinalFile, RetFile=" << RetFile << "\n";
    if (IsExistingFile(RetFile)) {
      std::cerr << "File RetFile=" << RetFile << " is already existing\n";
      throw TerminalException{1};
    }
    return RetFile;
  }
  if (IsExistingFile(eFile)) {
    std::cerr << "Error in FinalFile : eFile=" << eFile
              << " is already existing\n";
    throw TerminalException{1};
  }
  return eFile;
}

std::string FinalFileInScript(bool const &InPlaceRun,
                              std::string const &eFile) {
  if (InPlaceRun)
    return FILE_GetNakedFilename(eFile);
  return eFile;
}

PlotBound ReadPlotBound(FullNamelist const &eFull) {
  SingleBlock eBlPLOT = eFull.ListBlock.at("PLOT");
  bool VariableRange = eBlPLOT.ListBoolValues.at("VariableRange");
  std::vector<std::string> BoundSingle_var =
      eBlPLOT.ListListStringValues.at("BoundSingle_var");
  std::vector<double> BoundSingle_min =
      eBlPLOT.ListListDoubleValues.at("BoundSingle_min");
  std::vector<double> BoundSingle_max =
      eBlPLOT.ListListDoubleValues.at("BoundSingle_max");
  PlotBound eRecPlot;
  eRecPlot.VariableRange = VariableRange;
  eRecPlot.BoundSingle_var = BoundSingle_var;
  eRecPlot.BoundSingle_min = BoundSingle_min;
  eRecPlot.BoundSingle_max = BoundSingle_max;
  auto search = eBlPLOT.ListListStringValues.find("BoundDiff_var");
  if (search != eBlPLOT.ListListStringValues.end()) {
    std::vector<std::string> BoundDiff_var =
        eBlPLOT.ListListStringValues.at("BoundDiff_var");
    std::vector<double> BoundDiff_min =
        eBlPLOT.ListListDoubleValues.at("BoundDiff_min");
    std::vector<double> BoundDiff_max =
        eBlPLOT.ListListDoubleValues.at("BoundDiff_max");
    eRecPlot.BoundDiff_var = BoundDiff_var;
    eRecPlot.BoundDiff_min = BoundDiff_min;
    eRecPlot.BoundDiff_max = BoundDiff_max;
  }
  return eRecPlot;
}

DrawArr BasicArrayDraw(QuadArray const &eQuad) {
  DrawArr eDrw;
  eDrw.DoColorBar = false;
  eDrw.eQuadFrame = eQuad;
  eDrw.DoTitle = false;
  eDrw.ColorMap = "WhBlGrYeRe";
  eDrw.nbLevelSpa = 50;
  eDrw.cnFillMode = "RasterFill";
  eDrw.cnFillOn = false;
  eDrw.cnLinesOn = false;
  eDrw.cnLineLabelsOn = false;
  eDrw.nbLabelStride = 8;
  eDrw.DrawContourBathy = false;
  eDrw.GridResolution = "HighRes";
  eDrw.DrawRiver = false;
  return eDrw;
}

void PrintDrawArray(std::ostream &os, DrawArr const &eDrawArr) {
  os << "DoColorBar=" << eDrawArr.DoColorBar << "\n";
  os << "DoTitle=" << eDrawArr.DoTitle << "\n";
  os << "ColorMap=" << eDrawArr.ColorMap << "\n";
  os << "nbLevelSpa=" << eDrawArr.nbLevelSpa << "\n";
  os << "cnFillMode=" << eDrawArr.cnFillMode << "\n";
  os << "cnFillOn=" << eDrawArr.cnFillOn << "\n";
  os << "cnLinesOn=" << eDrawArr.cnLinesOn << "\n";
  os << "cnLineLabelsOn=" << eDrawArr.cnLineLabelsOn << "\n";
  os << "nbLabelStride=" << eDrawArr.nbLabelStride << "\n";
  os << "DrawContourBathy=" << eDrawArr.DrawContourBathy << "\n";
  os << "GridResolution=" << eDrawArr.GridResolution << "\n";
  os << "DrawRiver=" << eDrawArr.DrawRiver << "\n";
}

void ADD_RIVER(std::ostream &os, DrawArr const &eDrawArr) {
  if (eDrawArr.DrawRiver) {
    os << "  riv_data = "
          "asciiread(\"JadranRivers_extractor.dat\",(/10583,2/),\"float\")\n";
    os << "  lon=riv_data(:,0)\n";
    os << "  lat=riv_data(:,1)\n";
    os << "  segments=ind(lon.eq.-999)\n";
    os << "  ns=dimsizes(segments)\n";
    os << "  resP = True\n";
    os << "  resP@gsLineThicknessF = 1.5\n";
    os << "  resP@gsLineColor  = \"dodgerblue1\"\n";
    os << "  resP@tfPolyDrawOrder = \"PostDraw\"\n";
    os << "  lines = new(ns(0)-1,graphic)   ; array to hold polylines\n";
    os << "  do i=0,ns-2\n";
    os << "    xp=lon(segments(i)+2 : segments(i+1)-2)\n";
    os << "    yp=lat(segments(i)+2 : segments(i+1)-2)\n";
    os << "    lines(i)=gsn_add_polyline(wks,plot,xp,yp,resP)\n";
    os << "    delete(xp)\n";
    os << "    delete(yp)\n";
    os << "  end do\n";
    os << "  delete(segments)\n";
  }
}

void ADD_LISTLINESEGMENT(std::ostream &os, DrawArr const &eDrawArr) {
  if (eDrawArr.ListLineSegment.size() > 0) {
    os << "  LLS_ListLon = f->LLS_ListLon\n";
    os << "  LLS_ListLat = f->LLS_ListLat\n";
    os << "  LLS_ns = dimsizes(LLS_ListLon)\n";
    os << "  nbLine=LLS_ns(0)/2\n";
    os << "  resLine = True\n";
    os << "  resLine@gsLineThicknessF = 3.0\n";
    os << "  resLine@gsLineColor  = \"dodgerblue1\"\n";
    os << "  resLine@tfPolyDrawOrder = \"PostDraw\"\n";
    os << "  if (nbLine .gt. 0) then\n";
    os << "    linesRect = new(nbLine,graphic)\n";
    os << "    do iLine=0,nbLine-1\n";
    os << "      xp=LLS_ListLon(2*iLine : 2*iLine+1)\n";
    os << "      yp=LLS_ListLat(2*iLine : 2*iLine+1)\n";
    os << "      linesRect(iLine)=gsn_add_polyline(wks,plot,xp,yp,resLine)\n";
    os << "      delete(xp)\n";
    os << "      delete(yp)\n";
    os << "    end do\n";
    os << "  end if\n";
  }
}

void ADD_LISTMARKER(std::ostream &os, DrawArr const &eDrawArr) {
  if (eDrawArr.ListMarker.size() > 0) {
    os << "  LM_ListLon = f->LM_ListLon\n";
    os << "  LM_ListLat = f->LM_ListLat\n";
    os << "  LM_ListThick = f->LM_ListThick\n";
    os << "  LM_ns = dimsizes(LM_ListLon)\n";
    os << "  num_distinct_markers=LM_ns(0)\n";
    os << "  do i = 0, num_distinct_markers-1\n";
    os << "    gsres@gsMarkerColor      = \"dodgerblue1\"\n";
    os << "    gsres@gsMarkerThicknessF = LM_ListThick(i)\n";
    os << "    xp=(/LM_ListLon(i) /)\n";
    os << "    yp=(/LM_ListLat(i) /)\n";
    os << "    gsn_polymarker(wks, map, xp, yp, gsres)\n";
    os << "  end do\n";
  }
}

void ADD_ANNOTATION_TEXT(std::ostream &os, AnnotationRec const &TheAnnot) {
  if (TheAnnot.DrawAnnotation) {
    os << "  label=\"" << TheAnnot.AnnotationText << "\"\n";
    os << "  Xpos=" << TheAnnot.AnnotationLon << "\n";
    os << "  Ypos=" << TheAnnot.AnnotationLat << "\n";
    os << "  txres             = True                         ; Text resources "
          "desired\n";
    os << "  txres@txFont        = \"helvetica\"\n";
    os << "  txres@txFontHeightF=0.02\n";
    os << "  text = gsn_add_text(wks, plot, label, Xpos, Ypos, txres)\n";
  }
}

struct DrawScatterArr {
  std::string VarNameAB_file;
  // Drawing itself
  AnnotationRec TheAnnot;
  bool DoTitle;
  bool AddStatMeasModel;
  std::string NameA_plot;
  std::string NameB_plot;
  std::vector<double> data_rangeA;
  std::vector<double> data_rangeB;
  MyVector<double> eVectA;
  MyVector<double> eVectB;
  int aSize;
  int bSize;
};

//
// We need to deal with situation of computing scatter plot of
// measurement vs model (say Hs from altimeter vs model)
//   or
// Charnock vs U10.
void DEFINE_SCATTER_NC(std::string const &eFileNC,
                       DrawScatterArr const &eDrawScatter) {
  netCDF::NcFile dataFile(eFileNC, netCDF::NcFile::replace,
                          netCDF::NcFile::nc4);
  int nb = eDrawScatter.eVectA.size();
  int aSize = eDrawScatter.aSize;
  int bSize = eDrawScatter.bSize;
  netCDF::NcDim eDimOne = dataFile.addDim("one", 1);
  netCDF::NcDim eDimTwo = dataFile.addDim("two", 2);
  netCDF::NcDim eDimNb = dataFile.addDim("nb", nb);
  netCDF::NcDim eDimASize = dataFile.addDim("aSize", aSize);
  netCDF::NcDim eDimBSize = dataFile.addDim("bSize", bSize);

  std::vector<std::string> ListDimOne{"one"};
  std::vector<std::string> ListDimTwo{"two"};
  std::vector<std::string> ListDimNb{"nb"};
  std::vector<std::string> ListDimABsize{"aSize", "bSize"};

  netCDF::NcVar eVarData_rangeA =
      dataFile.addVar("data_rangeA", "double", ListDimTwo);
  netCDF::NcVar eVarData_rangeB =
      dataFile.addVar("data_rangeB", "double", ListDimTwo);
  netCDF::NcVar eVarListX = dataFile.addVar("ListX", "double", ListDimTwo);
  netCDF::NcVar eVarListY = dataFile.addVar("ListY", "double", ListDimTwo);
  netCDF::NcVar eVarX = dataFile.addVar("X", "double", ListDimNb);
  netCDF::NcVar eVarY = dataFile.addVar("Y", "double", ListDimNb);
  netCDF::NcVar eVarX2D = dataFile.addVar("X2D", "double", ListDimABsize);
  netCDF::NcVar eVarY2D = dataFile.addVar("Y2D", "double", ListDimABsize);
  netCDF::NcVar eVarCanvas = dataFile.addVar("canvas", "double", ListDimABsize);
  //
  std::vector<double> eFieldX(nb), eFieldY(nb);
  double TheMaxX = 0;
  double TheMaxY = 0;
  for (int i = 0; i < nb; i++) {
    double eVal = eDrawScatter.eVectA(i);
    if (eVal > TheMaxX)
      TheMaxX = eVal;
    eFieldX[i] = eVal;
  }
  for (int i = 0; i < nb; i++) {
    double eVal = eDrawScatter.eVectB(i);
    if (eVal > TheMaxY)
      TheMaxY = eVal;
    eFieldY[i] = eVal;
  }
  eVarX.putVar(eFieldX.data());
  eVarY.putVar(eFieldY.data());
  //
  double ePair[2];
  ePair[0] = 0;
  ePair[1] = TheMaxX;
  eVarListX.putVar(ePair);
  ePair[0] = 0;
  ePair[1] = TheMaxY;
  eVarListY.putVar(ePair);
  //
  std::vector<double> data_rangeA = eDrawScatter.data_rangeA;
  std::vector<double> data_rangeB = eDrawScatter.data_rangeB;
  double eFrangeA[2];
  for (int i = 0; i < 2; i++)
    eFrangeA[i] = data_rangeA[i];
  eVarData_rangeA.putVar(eFrangeA);
  double eFrangeB[2];
  for (int i = 0; i < 2; i++)
    eFrangeB[i] = data_rangeB[i];
  eVarData_rangeB.putVar(eFrangeB);
  //
  std::vector<double> X2D(aSize * bSize), Y2D(aSize * bSize);
  double deltaA =
      (data_rangeA[1] - data_rangeA[0]) / static_cast<double>(aSize - 1);
  double deltaB =
      (data_rangeB[1] - data_rangeB[0]) / static_cast<double>(bSize - 1);
  int idx = 0;
  for (int iA = 0; iA < aSize; iA++)
    for (int iB = 0; iB < bSize; iB++) {
      double eA = data_rangeA[0] + iA * deltaA;
      double eB = data_rangeB[0] + iB * deltaB;
      X2D[idx] = eA;
      Y2D[idx] = eB;
      idx++;
    }
  eVarX2D.putVar(X2D.data());
  eVarY2D.putVar(Y2D.data());
  //
  MyMatrix<int> canvasInt(aSize, bSize);
  for (int iA = 0; iA < aSize; iA++)
    for (int iB = 0; iB < bSize; iB++)
      canvasInt(iA, iB) = 0;
  for (int iEnt = 0; iEnt < nb; iEnt++) {
    double eA = eDrawScatter.eVectA(iEnt);
    double eB = eDrawScatter.eVectB(iEnt);
    double iA_d = (eA - data_rangeA[0]) / deltaA;
    double iB_d = (eB - data_rangeB[0]) / deltaB;
    int iA = static_cast<int>(floor(iA_d));
    int iB = static_cast<int>(floor(iB_d));
    if (iA >= 0 && iA < aSize && iB >= 0 && iB < bSize)
      canvasInt(iA, iB)++;
  }
  std::vector<double> canvas(aSize * bSize);
  double MissVal = 0;
  idx = 0;
  for (int iA = 0; iA < aSize; iA++)
    for (int iB = 0; iB < bSize; iB++) {
      double eVal;
      int eValI = canvasInt(iA, iB);
      if (eValI == 0)
        eVal = MissVal;
      else
        eVal = log10(static_cast<double>(eValI));
      //      std::cerr << "iA=" << iA << " iB=" << iB << " eVal=" << eVal <<
      //      "\n";
      canvas[idx] = eVal;
      idx++;
    }
  eVarCanvas.putVar(canvas.data());
  //
}

void PrintMyScriptSubtitle(std::ostream &os) {
  os << "procedure subtitles(wks:graphic,plot:graphic,lstr:string,cstr:string, "
        "\\\n";
  os << "                    rstr:string,tres)\n";
  os << "local txres, font_height, amres\n";
  os << "begin\n";
  os << "  if(tres) then\n";
  os << "    txres = tres     ; Copy resources\n";
  os << "  else\n";
  os << "    txres = True\n";
  os << "  end if\n";
  os << "  ;\n";
  os << "  ; Retrieve font height of left axis string and use to\n";
  os << "  ; calculate size of subtitles.\n";
  os << "  ;\n";
  os << "  if(.not.isatt(txres,\"txFontHeightF\")) then\n";
  os << "    getvalues plot\n";
  os << "      \"tiXAxisFontHeightF\" : font_height\n";
  os << "    end getvalues\n";
  os << "    txres@txFontHeightF = font_height*0.9\n";
  os << "  end if\n";
  os << "  ;\n";
  os << "  ; Set some some annotation resources.\n";
  os << "  ;\n";
  os << "  amres                  = True\n";
  os << "  if(.not.isatt(txres,\"amOrthogonalPosF\")) then\n";
  os << "    amres@amOrthogonalPosF = -0.53   ; Top of plot plus a little "
        "extra\n";
  os << "                                     ; to stay out of the "
        "tickmarks.\n";
  os << "  else\n";
  os << "    amres@amOrthogonalPosF = txres@amOrthogonalPosF\n";
  os << "  end if\n";
  os << "  ;\n";
  os << "  ; Create three strings to put at the top, using a slightly\n";
  os << "  ; smaller font height than the axis titles.\n";
  os << "  ;\n";
  os << "  if(lstr.ne.\"\") then\n";
  os << "    txidl = gsn_create_text(wks, lstr, txres)\n";
  os << "    amres@amJust           = \"BottomLeft\"\n";
  os << "    amres@amParallelPosF   = -0.5   ; Left-justified\n";
  os << "    annoidl = gsn_add_annotation(plot, txidl, amres)\n";
  os << "  end if\n";
  os << "  if(cstr.ne.\"\") then\n";
  os << "    txidc = gsn_create_text(wks, cstr, txres)\n";
  os << "    amres@amJust           = \"BottomCenter\"\n";
  os << "    amres@amParallelPosF   = 0.0   ; Centered\n";
  os << "    annoidc = gsn_add_annotation(plot, txidc, amres)\n";
  os << "  end if\n";
  os << "  if(rstr.ne.\"\") then\n";
  os << "    txidr = gsn_create_text(wks, rstr, txres)\n";
  os << "    amres@amJust           = \"BottomRight\"\n";
  os << "    amres@amParallelPosF   = 0.5   ; Right-justifed\n";
  os << "    annoidr = gsn_add_annotation(plot, txidr, amres)\n";
  os << "  end if\n";
  os << "end\n";
}

struct RealInfoNclExtension {
  std::string eProgReal;
  std::string eExtensionReal;
};

RealInfoNclExtension GetRealInfoNclExtension(std::string const &eExtension) {
  if (eExtension == "pdfcrop")
    return {"ncl_crop", "pdf"};
  if (eExtension == "pngcrop")
    return {"ncl_crop", "png"};
  return {"ncl", eExtension};
}

void PLOT_SCATTER(DrawScatterArr const &eDrawScatter,
                  NCLcaller<GeneralType> &eCall,
                  PermanentInfoDrawing const &ePerm) {
  RealInfoNclExtension eReal = GetRealInfoNclExtension(ePerm.Extension);
  std::string VarNameAB_file = eDrawScatter.VarNameAB_file;
  std::string FileName = ePerm.eDir + VarNameAB_file;
  std::string TargetFile = FileName + "_storsave." + eReal.eExtensionReal;
  std::string TitleStr;
  bool InPlaceRun = ePerm.eBChoice.InPlaceRun;
  std::string eFileNC = FinalFile(InPlaceRun, TargetFile,
                                  ePerm.PrefixTemp.str() + "DataScatter_" +
                                      VarNameAB_file + ".nc");
  std::string eFileNCL = FinalFile(InPlaceRun, TargetFile,
                                   ePerm.PrefixTemp.str() + "ScriptScatter_" +
                                       VarNameAB_file + ".ncl");
  DEFINE_SCATTER_NC(eFileNC, eDrawScatter);
  //
  std::ofstream OUTncl(eFileNCL);
  OUTncl << "load \"$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl\"\n";
  OUTncl << "load \"$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl\"\n";
  PrintMyScriptSubtitle(OUTncl);
  //  OUTncl << "load \"$NCARG_ROOT/myscript/subtitles.ncl\"\n";
  OUTncl << "begin\n";
  OUTncl << "  f = addfile(\"" << FinalFileInScript(InPlaceRun, eFileNC)
         << "\", \"r\")\n";
  OUTncl << "  X = f->X2D\n";
  OUTncl << "  Y = f->Y2D\n";
  OUTncl << "  Vtot = f->canvas\n";
  OUTncl << "  stride=20\n";
  OUTncl << "  spacing=0.025\n";
  OUTncl << "  minValue=min(Vtot)\n";
  OUTncl << "  maxValue=2\n";
  OUTncl << "  wks  = gsn_open_wks (\"" << eReal.eExtensionReal << "\",\""
         << FinalFileInScript(InPlaceRun, FILE_RemoveExtension(TargetFile))
         << "\")\n";
  OUTncl << "  gsn_define_colormap(wks,\"BlAqGrYeOrRevi200\")\n";
  OUTncl << "  vres1 = True               ; plot mods desired\n";
  OUTncl << "  vres1@gsnDraw   = False\n";
  OUTncl << "  vres1@gsnFrame  = False\n";
  OUTncl << "  vres1@gsnMaximize     = True    ; Maximize plot in frame\n";
  OUTncl << "  vres1@gsnPaperOrientation = \"" << ePerm.eDrawArr.LandPortr
         << "\"\n";
  OUTncl << "  vres1@cnFillDrawOrder        = \"PreDraw\"\n";
  OUTncl << "  vres1@cnFillOn             = "
         << NCL_bool(ePerm.eDrawArr.cnFillOn)
         << "   ; turn on color for contours\n";
  OUTncl << "  vres1@cnLinesOn            = "
         << NCL_bool(ePerm.eDrawArr.cnLinesOn)
         << "    ; turn off contour lines\n";
  OUTncl << "  vres1@cnLineLabelsOn       = "
         << NCL_bool(ePerm.eDrawArr.cnLineLabelsOn)
         << "   ; turn off contour line labels\n";
  OUTncl << "  vres1@cnFillMode           = \"AreaFill\"\n";
  OUTncl << "  vres1@gsnSpreadColors      = True               ; use full "
            "color map\n";
  OUTncl << "  vres1@gsnSpreadColorEnd     = -2\n";
  OUTncl << "  LLabel=new(81,string)\n";
  OUTncl << "  LLabel(0)=\"1\"\n";
  OUTncl << "  LLabel(20)=\"3\"\n";
  OUTncl << "  LLabel(40)=\"10\"\n";
  OUTncl << "  LLabel(60)=\"30\"\n";
  OUTncl << "  LLabel(80)=\"100\"\n";
  OUTncl << "  vres1@lbLabelBarOn =  True\n";
  OUTncl << "  vres1@lbLabelStride            = stride\n";
  OUTncl << "  vres1@lbOrientation        = \"Vertical\"     ; Vertical label "
            "bar\n";
  OUTncl << "  vres1@lbLabelStrings = LLabel\n";
  OUTncl << "  vres1@cnLevelSelectionMode = \"ManualLevels\"     ; set manual "
            "contour levels\n";
  OUTncl << "  vres1@cnMinLevelValF       = minValue                ; set min "
            "contour level\n";
  OUTncl << "  vres1@cnMaxLevelValF       = maxValue              ; set max "
            "contour level\n";
  OUTncl << "  vres1@cnLevelSpacingF      = spacing                 ; set "
            "contour spacing\n";
  OUTncl << "  vres1@pmLabelBarOrthogonalPosF = -0.01          ; move label "
            "bar closer\n";
  OUTncl << "  ;  vres1@pmLabelBarDisplayMode = \"Always\"          ; Turn on "
            "a label bar.\n";
  OUTncl << "  vres1@lbPerimOn             = False             ; no box around "
            "it\n";
  OUTncl << "  vres1@lbBoxLinesOn         = False               ; Yes/No "
            "labelbar box lines.\n";
  OUTncl << "  vres1@tiXAxisString  = \"" << eDrawScatter.NameA_plot << "\"\n";
  OUTncl << "  vres1@tiYAxisString    = \"" << eDrawScatter.NameB_plot
         << "\"\n";
  OUTncl << "  vres1@tiXAxisOffsetYF = 0.0\n";
  OUTncl << "  ; vres1@tmYLPrecision = 0\n";
  OUTncl << "  ; First part Hwave\n";
  OUTncl << "  vres1@sfXArray            = X\n";
  OUTncl << "  vres1@sfYArray            = Y\n";
  OUTncl << "  vres1@trGridType          = \"TriangularMesh\"\n";
  OUTncl << "  ;  vres1@trGridType         = \"curvilinear\"\n";
  OUTncl << "  plot = gsn_csm_contour(wks,Vtot,vres1)\n";
  OUTncl << "  resP = True\n";
  OUTncl << "  resP@gsLineThicknessF = 1.5\n";
  OUTncl << "  resP@gsLineColor  = \"black\"\n";
  OUTncl << "  resP@tfPolyDrawOrder = \"PostDraw\"\n";
  if (eDrawScatter.AddStatMeasModel) {
    OUTncl << "  xp=(/-180, 180/)\n";
    OUTncl << "  yp=(/-180, 180/)\n";
    OUTncl << "  line0=gsn_add_polyline(wks,plot,xp,yp,resP)\n";
  }
  OUTncl << "  data_rangeA = f->data_rangeA\n";
  OUTncl << "  data_rangeB = f->data_rangeB\n";
  OUTncl << "  line1=gsn_add_polyline(wks,plot,data_rangeA,data_rangeB,resP)\n";
  ADD_ANNOTATION_TEXT(OUTncl, eDrawScatter.TheAnnot);
  if (eDrawScatter.AddStatMeasModel) {
    T_stat eStat =
        ComputeStatistics_MyVector(eDrawScatter.eVectA, eDrawScatter.eVectB);
    std::string optStatStr = ePerm.eDrawArr.optStatStr;
    T_statString eStatStr =
        ComputeStatisticString_from_Statistics(eStat, optStatStr);
    std::string strWrite = "m=" + eStatStr.strSlope +
                           " c=" + eStatStr.strCorrelation +
                           " s=" + eStatStr.strScatterIndex;
    OUTncl << "  txresB             = True\n";
    OUTncl << "  txresB@txFontHeightF = 0.02\n";
    OUTncl << "  txresB@txFontColor = \"black\"\n";
    OUTncl << "  strLeft=\"\"\n";
    OUTncl << "  strMid=\"" << strWrite << "\"\n";
    OUTncl << "  strRight=\"\"\n";
    OUTncl << "  subtitles(wks, plot, strLeft, strMid, strRight, txresB)\n";
  }
  OUTncl << "  draw(plot)\n";
  OUTncl << "  frame(wks)\n";
  OUTncl << "end\n";
  OUTncl.close();
  //
  TripleNCL eTripl{eReal.eProgReal, TargetFile, eFileNC, eFileNCL,
                   ePerm.eBChoice};
  GeneralType eGen(eTripl);
  eCall.SubmitJob(eGen);
}

void ADD_LISTLINESEGMENT_NC(
    netCDF::NcFile &dataFile,
    std::vector<SeqLineSegment> const &ListLineSegment) {
  int TotalLen = 0;
  for (auto &eSeqLineSegment : ListLineSegment) {
    int len = eSeqLineSegment.ListPairLL.size();
    if (!eSeqLineSegment.IsClosed)
      len--;
    if (len < 0) {
      std::cerr << "We should have len >= 0. len=" << len << "\n";
      throw TerminalException{1};
    }
    TotalLen += len;
  }
  int RelSiz = 2 * TotalLen;
  //  std::cerr << "RelSiz=" << RelSiz << "\n";
  std::vector<double> ListLon(RelSiz), ListLat(RelSiz);
  int idx = 0;
  for (auto &eSeqLineSegment : ListLineSegment) {
    int len = eSeqLineSegment.ListPairLL.size();
    for (int i = 0; i < len - 1; i++) {
      double eLon1 = eSeqLineSegment.ListPairLL[i].eLon;
      double eLat1 = eSeqLineSegment.ListPairLL[i].eLat;
      double eLon2 = eSeqLineSegment.ListPairLL[i + 1].eLon;
      double eLat2 = eSeqLineSegment.ListPairLL[i + 1].eLat;
      ListLon[idx] = eLon1;
      ListLat[idx] = eLat1;
      idx++;
      ListLon[idx] = eLon2;
      ListLat[idx] = eLat2;
      idx++;
    }
    if (eSeqLineSegment.IsClosed) {
      double eLon1 = eSeqLineSegment.ListPairLL[0].eLon;
      double eLat1 = eSeqLineSegment.ListPairLL[0].eLat;
      double eLon2 = eSeqLineSegment.ListPairLL[len - 1].eLon;
      double eLat2 = eSeqLineSegment.ListPairLL[len - 1].eLat;
      ListLon[idx] = eLon1;
      ListLat[idx] = eLat1;
      idx++;
      ListLon[idx] = eLon2;
      ListLat[idx] = eLat2;
      idx++;
    }
  }
  netCDF::NcDim eDimMnp = dataFile.addDim("RelSiz", RelSiz);
  std::vector<std::string> ListDim = {"RelSiz"};
  netCDF::NcVar eVarListLon = dataFile.addVar("LLS_ListLon", "double", ListDim);
  netCDF::NcVar eVarListLat = dataFile.addVar("LLS_ListLat", "double", ListDim);
  eVarListLon.putVar(ListLon.data());
  eVarListLat.putVar(ListLat.data());
}

void ADD_LISTMARKER_NC(netCDF::NcFile &dataFile,
                       std::vector<SingleMarker> const &ListMarker) {
  int nbMarker = ListMarker.size();
  std::vector<double> ListLon(nbMarker), ListLat(nbMarker), ListThick(nbMarker);
  for (int i = 0; i < nbMarker; i++) {
    ListLon[i] = ListMarker[i].Coord.eLon;
    ListLat[i] = ListMarker[i].Coord.eLat;
    ListThick[i] = ListMarker[i].thickness;
  }
  netCDF::NcDim eDimMnp = dataFile.addDim("nbMarker", nbMarker);
  std::vector<std::string> ListDim = {"nbMarker"};
  netCDF::NcVar eVarListLon = dataFile.addVar("LM_ListLon", "double", ListDim);
  netCDF::NcVar eVarListLat = dataFile.addVar("LM_ListLat", "double", ListDim);
  netCDF::NcVar eVarListThick =
      dataFile.addVar("LM_ListThick", "double", ListDim);
  eVarListLon.putVar(ListLon.data());
  eVarListLat.putVar(ListLat.data());
  eVarListThick.putVar(ListThick.data());
}

void DEFINE_QUIVER_NC(std::string const &eFileNC, GridArray const &GrdArr,
                      MyMatrix<double> const &U_rho,
                      MyMatrix<double> const &V_rho,
                      MyMatrix<double> const &F_rho,
                      std::vector<SeqLineSegment> const &ListLineSegment,
                      std::vector<SingleMarker> const &ListMarker) {
  int idx;
  netCDF::NcFile dataFile(eFileNC, netCDF::NcFile::replace,
                          netCDF::NcFile::nc4);
  int eta = GrdArr.GrdArrRho.LON.rows();
  int xi = GrdArr.GrdArrRho.LON.cols();
  int eta_lat = GrdArr.GrdArrRho.LAT.rows();
  int xi_lat = GrdArr.GrdArrRho.LAT.cols();
  int eta_U = U_rho.rows();
  int xi_U = U_rho.cols();
  int eta_V = V_rho.rows();
  int xi_V = V_rho.cols();
  if (eta != eta_lat || eta != eta_U || eta != eta_V) {
    std::cerr << "All dimensions should be the same.\n";
    std::cerr << "Now we have following:\n";
    std::cerr << "eta_lon=" << eta << "\n";
    std::cerr << "eta_lat=" << eta_lat << "\n";
    std::cerr << "eta_U=" << eta_U << "\n";
    std::cerr << "eta_V=" << eta_V << "\n";
    throw TerminalException{1};
  }
  if (xi != xi_lat || xi != xi_U || xi != xi_V) {
    std::cerr << "All dimensions should be the same.\n";
    std::cerr << "Now we have following:\n";
    std::cerr << "xi_lon=" << xi << "\n";
    std::cerr << "xi_lat=" << xi_lat << "\n";
    std::cerr << "xi_U=" << xi_U << "\n";
    std::cerr << "xi_V=" << xi_V << "\n";
    throw TerminalException{1};
  }
  std::string typeName = "double";
  std::string typeNameInt = "int";
  std::string eEta = "eta_rho";
  std::string eXi = "xi_rho";
  std::string Lon = "lon";
  std::string Lat = "lat";
  std::string Uvar = "u";
  std::string Vvar = "v";
  std::string Fvar = "F";
  if (xi > 1) {
    bool ApplyCritValue = true;
    double eCritValue = -10000;
    double dataMiss[1];
    dataMiss[0] = eCritValue;
    std::string MissVal = "_FillValue";
    //
    std::vector<double> valLON(eta * xi), valLAT(eta * xi), valU(eta * xi),
        valV(eta * xi), valF(eta * xi);
    idx = 0;
    //
    netCDF::NcDim eDimEta = dataFile.addDim(eEta, eta);
    netCDF::NcDim eDimXi = dataFile.addDim(eXi, xi);
    bool PositiveOrient = false;
    std::vector<std::string> ListDim;
    if (PositiveOrient) {
      ListDim = {eEta, eXi};
    } else {
      ListDim = {eXi, eEta};
    }
    netCDF::NcVar eVarLON = dataFile.addVar(Lon, typeName, ListDim);
    netCDF::NcVar eVarLAT = dataFile.addVar(Lat, typeName, ListDim);
    //
    netCDF::NcVar eVarU = dataFile.addVar(Uvar, typeName, ListDim);
    eVarU.putAtt(MissVal, netCDF::NcType::nc_DOUBLE, 1, dataMiss);
    //
    netCDF::NcVar eVarV = dataFile.addVar(Vvar, typeName, ListDim);
    eVarV.putAtt(MissVal, netCDF::NcType::nc_DOUBLE, 1, dataMiss);
    //
    netCDF::NcVar eVarF = dataFile.addVar(Fvar, typeName, ListDim);
    eVarF.putAtt(MissVal, netCDF::NcType::nc_DOUBLE, 1, dataMiss);
    //
    // We change the orientation of the grid
    // This is needed for plots of ideal test cases which are not in LON/LAT
    // For those the LON/LAT are simply not read and so the result become
    // inverted plots.
    //
    if (PositiveOrient) {
      for (int i = 0; i < eta; i++)
        for (int j = 0; j < xi; j++) {
          valLON[idx] = GrdArr.GrdArrRho.LON(i, j);
          valLAT[idx] = GrdArr.GrdArrRho.LAT(i, j);
          if (GrdArr.GrdArrRho.MSK(i, j) == 1 || !ApplyCritValue) {
            valU[idx] = U_rho(i, j);
            valV[idx] = V_rho(i, j);
            valF[idx] = F_rho(i, j);
          } else {
            valU[idx] = eCritValue;
            valV[idx] = eCritValue;
            valF[idx] = eCritValue;
          }
          idx++;
        }
    } else {
      for (int j = 0; j < xi; j++)
        for (int i = 0; i < eta; i++) {
          valLON[idx] = GrdArr.GrdArrRho.LON(i, j);
          valLAT[idx] = GrdArr.GrdArrRho.LAT(i, j);
          if (GrdArr.GrdArrRho.MSK(i, j) == 1 || !ApplyCritValue) {
            valU[idx] = U_rho(i, j);
            valV[idx] = V_rho(i, j);
            valF[idx] = F_rho(i, j);
          } else {
            valU[idx] = eCritValue;
            valV[idx] = eCritValue;
            valF[idx] = eCritValue;
          }
          idx++;
        }
    }
    eVarLON.putVar(valLON.data());
    eVarLAT.putVar(valLAT.data());
    eVarU.putVar(valU.data());
    eVarV.putVar(valV.data());
    eVarF.putVar(valF.data());
  } else {
    int mnp = eta;
    int mne = GrdArr.INE.rows();
    std::string eMnp = "mnp";
    std::string eMne = "mne";
    std::string eThree = "three";
    netCDF::NcDim eDimMnp = dataFile.addDim(eMnp, mnp);
    netCDF::NcDim eDimThree = dataFile.addDim(eThree, 3);
    netCDF::NcDim eDimMne = dataFile.addDim(eMne, mne);
    std::vector<std::string> ListDim = {eMnp};
    auto writeMnpVar = [&](std::string const &name,
                           MyVector<double> const &VAR) -> void {
      netCDF::NcVar eVar = dataFile.addVar(name, typeName, ListDim);
      std::vector<double> A(eta);
      for (int i = 0; i < eta; i++)
        A[i] = VAR(i, 0);
      eVar.putVar(A.data());
    };
    writeMnpVar(Lon, GrdArr.GrdArrRho.LON);
    writeMnpVar(Lat, GrdArr.GrdArrRho.LAT);
    writeMnpVar(Uvar, U_rho);
    writeMnpVar(Vvar, V_rho);
    writeMnpVar(Fvar, F_rho);
    std::string Fine = "ele";
    std::vector<std::string> ListDimINE = {eMne, eThree};
    netCDF::NcVar eVarINE = dataFile.addVar(Fine, typeNameInt, ListDimINE);
    //
    std::vector<int> valI(3 * mne);
    idx = 0;
    for (int ie = 0; ie < mne; ie++)
      for (int i = 0; i < 3; i++) {
        int eConn = GrdArr.INE(ie, i);
        valI[idx] = eConn;
        idx++;
      }
    eVarINE.putVar(valI.data());
  }
  if (ListLineSegment.size() > 0)
    ADD_LISTLINESEGMENT_NC(dataFile, ListLineSegment);
  if (ListMarker.size() > 0)
    ADD_LISTMARKER_NC(dataFile, ListMarker);
}

void PLOT_QUIVER(std::string const &FileName, GridArray const &GrdArr,
                 DrawArr const &eDrawArr, RecVar const &eRecVar,
                 NCLcaller<GeneralType> &eCall,
                 PermanentInfoDrawing const &ePerm) {
  RealInfoNclExtension eReal = GetRealInfoNclExtension(ePerm.Extension);
  std::string TargetFile = FileName + "_storsave." + eReal.eExtensionReal;
  RecSymbolic RecS = eRecVar.RecS;
  bool InPlaceRun = ePerm.eBChoice.InPlaceRun;
  std::string eFileNC =
      FinalFile(InPlaceRun, TargetFile,
                ePerm.PrefixTemp.str() + "DataQuiver_" + eDrawArr.VarNameUF +
                    "_" + RecS.strAll + ".nc");
  std::string eFileNCL =
      FinalFile(InPlaceRun, TargetFile,
                ePerm.PrefixTemp.str() + "ScriptQuiver_" + eDrawArr.VarNameUF +
                    "_" + RecS.strAll + ".ncl");
  DEFINE_QUIVER_NC(eFileNC, GrdArr, eRecVar.U, eRecVar.V, eRecVar.F,
                   eDrawArr.ListLineSegment, eDrawArr.ListMarker);
  bool IsSpherical = GrdArr.IsSpherical;
  //  std::cerr << "PLOT_QUIVER, IsSpherical=" << IsSpherical << "\n";
  //
  std::ofstream OUTncl(eFileNCL);
  OUTncl << "load \"$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl\"\n";
  OUTncl << "load \"$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl\"\n";
  OUTncl << "begin\n";
  OUTncl << "  ; \n";
  OUTncl << "  ; Loading fundamental data\n";
  OUTncl << "  ; \n";
  OUTncl << "  f = addfile(\"" << FinalFileInScript(InPlaceRun, eFileNC)
         << "\", \"r\")\n";
  OUTncl << "  lat2d  = f->lat\n";
  OUTncl << "  lon2d  = f->lon\n";
  OUTncl << "  eVarU  = f->u\n";
  OUTncl << "  eVarV  = f->v\n";
  OUTncl << "  speed  = f->F\n";
  OUTncl << "  wks  = gsn_open_wks (\"" << eReal.eExtensionReal << "\",\""
         << FinalFileInScript(InPlaceRun, FILE_RemoveExtension(TargetFile))
         << "\")\n";
  OUTncl << "  vres1=True\n";
  OUTncl << "  vres1@gsnDraw = False\n";
  OUTncl << "  vres1@gsnFrame = False\n";
  OUTncl << "  vres1@gsnPaperOrientation = \"" << eDrawArr.LandPortr << "\"\n";
  OUTncl << "  vres1@gsnScalarContour     = True\n";
  OUTncl << "  vres1@gsnSpreadColors      = True\n";
  OUTncl << "  vres1@gsnSpreadColorEnd    = -3\n";
  OUTncl << "  ;  vres1@gsnSpreadColorStart  = 1\n";
  OUTncl << "  vres1@mpLandFillColor      = \"grey\"\n";
  OUTncl << "  vres1@cnFillDrawOrder      = \"PreDraw\"\n";
  if (eDrawArr.FillLand) {
    OUTncl << "  vres1@cnFillOn             = True\n";
  } else {
    OUTncl << "  vres1@cnFillOn             = False\n";
  }
  OUTncl << "  vres1@cnLinesOn            = " << NCL_bool(eDrawArr.cnLinesOn)
         << "\n";
  OUTncl << "  vres1@cnLineLabelsOn       = "
         << NCL_bool(eDrawArr.cnLineLabelsOn)
         << "   ; turn off contour line labels\n";
  OUTncl << "  vres1@cnLevelSelectionMode = \"ManualLevels\"\n";
  OUTncl << "  vres1@cnMinLevelValF       = " << RecS.minval << "\n";
  OUTncl << "  vres1@cnMaxLevelValF       = " << RecS.maxval << "\n";
  OUTncl << "  vres1@lbLabelStride        = 8\n";
  int nbLevelSpa = eDrawArr.nbLevelSpa;
  double TheLevelSpa =
      (RecS.maxval - RecS.minval) / static_cast<double>(nbLevelSpa);
  OUTncl << "  vres1@cnLevelSpacingF      = " << TheLevelSpa << "\n";
  OUTncl << "  ;\n";
  OUTncl << "  ; The vertical label bar on the right\n";
  OUTncl << "  ;\n";
  OUTncl << "  vres1@lbOrientation        = \"Vertical\"     ; Vertical label "
            "bar\n";
  OUTncl << "  vres1@pmLabelBarOrthogonalPosF = 0.025          ; move label "
            "bar closer\n";
  OUTncl << "  vres1@pmLabelBarDisplayMode = \"Always\"          ; Turn on a "
            "label bar.\n";
  OUTncl << "  vres1@lbPerimOn             = False             ; no box around "
            "it\n";
  OUTncl << "  vres1@lbBoxLinesOn         = False               ; Yes/No "
            "labelbar box lines.\n";
  OUTncl << "  vres1@lbTitleString    = \"" << RecS.VarName1 << " ["
         << RecS.Unit << "]\"\n";
  OUTncl << "  vres1@lbTitleFont      = \"Helvetica\"\n";
  OUTncl << "  vres1@lbTitleFontHeightF = 0.015\n";
  OUTncl << "  vres1@lbTitleDirection     = \"Across\" \n";
  OUTncl << "  vres1@lbTitlePosition = \"Right\"\n";
  OUTncl << "  vres1@lbTitleAngleF = 90\n";
  if (IsSpherical) {
    OUTncl << "  vres1@mpProjection = \"Mercator\"\n";
    OUTncl << "  vres1@mpLimitMode         = \"Corners\"             ; choose "
              "range of map\n";
    OUTncl << "  vres1@mpLeftCornerLatF    = " << eDrawArr.eQuadFrame.MinLat
           << "\n";
    OUTncl << "  vres1@mpLeftCornerLonF    = " << eDrawArr.eQuadFrame.MinLon
           << "\n";
    OUTncl << "  vres1@mpRightCornerLatF   = " << eDrawArr.eQuadFrame.MaxLat
           << "\n";
    OUTncl << "  vres1@mpRightCornerLonF   = " << eDrawArr.eQuadFrame.MaxLon
           << "\n";
    if (eDrawArr.FillLand) {
      OUTncl << "  vres1@mpFillOn      = True\n";
      OUTncl << "  vres1@mpLandFillColor       = \"gray\"            ; set "
                "land to be gray\n";
    } else {
      OUTncl << "  vres1@mpFillOn      = False\n";
    }
    OUTncl << "  vres1@mpDataBaseVersion      = \"" << eDrawArr.GridResolution
           << "\"          ; use high resolution coast. Outcomment and get "
              "rough coastline\n";
  } else {
    OUTncl << "  vres1@sfXArray = f->lon\n";
    OUTncl << "  vres1@sfYArray = f->lat\n";
    if (GrdArr.IsFE == 1)
      OUTncl << "  vres1@trGridType = \"TriangularMesh\"\n";
  }
  OUTncl << "  vres1@pmTickMarkDisplayMode  = \"Always\"           ; turn on "
            "tickmarks\n";
  OUTncl << "  vres1@vcMonoLineArrowColor  = True             ; colored by "
            "their mag\n";
  OUTncl << "  vres1@pmLabelBarWidthF = 0.03\n";
  if (eDrawArr.DoTitle) {
    OUTncl << "  vres1@tiMainString    = \"" << eDrawArr.TitleStr << "\"\n";
    OUTncl << "  vres1@tiMainFont      = \"Helvetica\"\n";
    OUTncl << "  vres1@tiMainFontHeightF=0.015\n";
  }
  OUTncl << "  ; vres1@gsnRightString    = \"Sea surface elevation\"\n";
  OUTncl << "  ; vres1@gsnLeftString     = \"Difference\"\n";
  OUTncl << "  ;  vres1@vcRefMagnitudeF   = 10.\n";
  OUTncl << "  vres1@vcGlyphStyle      = \"CurlyVector\"\n";
  OUTncl << "  ;  vres1@vcGlyphStyle      = \"LineArrow\"\n";
  OUTncl << "  vres1@vcMinDistanceF    = 0.014\n";
  OUTncl << "  vres1@vcRefLengthF      = " << eDrawArr.vcRefLengthF << "\n";
  OUTncl << "  vres1@vcRefAnnoOn       = False\n";
  OUTncl << "  vres1@vcLineArrowHeadMaxSizeF = 0.005\n";
  OUTncl << "  ;  vres1@vcRefAnnoOrthogonalPosF = -0.8\n";
  OUTncl << "  ;  vres1@vcRefAnnoParallelPosF = 1.0\n";
  OUTncl << "  ;  vres1@vcRefMagnitudeF         = 10\n";
  OUTncl << "  ;  vres1@vcRefAnnoString1 = \"10 m/s\"\n";
  OUTncl << "  ;  gsn_define_colormap (wks,\"testcmap\")\n";
  OUTncl << "  ;\n";
  OUTncl << "  ; Specifying the colormap\n";
  OUTncl << "  ;\n";
  OUTncl << "  ;  gsn_define_colormap (wks,\"" << eDrawArr.ColorMap << "\")\n";
  OUTncl << "  gsn_define_colormap (wks,\"BlAqGrYeOrReVi200\")\n";
  OUTncl << "  ;  gsn_define_colormap(wks,\"BlGrYeOrReVi200\")\n";
  OUTncl << "  ;  gsn_define_colormap(wks,\"hotres\")\n";
  OUTncl << "  ;  gsn_define_colormap(wks,\"rainbow\")\n";
  OUTncl << "  ;  gsn_define_colormap(wks,\"ViBlGrWhYeOrRe\")\n";
  OUTncl << "  ;  gsn_define_colormap(wks,\"BlWhRe\")\n";
  OUTncl << "  ;  gsn_define_colormap(wks,\"GrayWhiteGray\")\n";
  OUTncl << "  ;  gsn_define_colormap(wks,\"BlWhRe\")\n";
  OUTncl << "  i = NhlNewColor(wks,0.8,0.8,0.8)      ; add gray to colormap\n";
  OUTncl << "  i = NhlNewColor(wks,0.9,0.9,0.9)      ; add gray to colormap\n";
  OUTncl << "  eVarU@lat2d=lat2d\n";
  OUTncl << "  eVarU@lon2d=lon2d\n";
  OUTncl << "  eVarV@lat2d=lat2d\n";
  OUTncl << "  eVarV@lon2d=lon2d\n";
  OUTncl << "  speed@lat2d=lat2d\n";
  OUTncl << "  speed@lon2d=lon2d\n";
  bool UseFvar = false;
  if (UseFvar) {
    OUTncl << "  fVarU=eVarU/speed\n";
    OUTncl << "  fVarV=eVarV/speed\n";
    OUTncl << "  fVarU@lat2d=lat2d\n";
    OUTncl << "  fVarU@lon2d=lon2d\n";
    OUTncl << "  fVarV@lat2d=lat2d\n";
    OUTncl << "  fVarV@lon2d=lon2d\n";
    OUTncl
        << "  plot = gsn_csm_vector_scalar_map(wks,fVarU,fVarV,speed,vres1)\n";
  } else {
    if (IsSpherical) {
      OUTncl << "  plot = "
                "gsn_csm_vector_scalar_map(wks,eVarU,eVarV,speed,vres1)\n";
    } else {
      OUTncl << "  plot = gsn_csm_vector_scalar(wks,eVarU,eVarV,speed,vres1)\n";
    }
  }
  ADD_LISTLINESEGMENT(OUTncl, eDrawArr);
  ADD_LISTMARKER(OUTncl, eDrawArr);
  ADD_RIVER(OUTncl, eDrawArr);
  ADD_ANNOTATION_TEXT(OUTncl, eDrawArr.TheAnnot);
  for (auto &eLine : eDrawArr.ListInsertLines)
    OUTncl << eLine << "\n";
  OUTncl << "  draw(plot)\n";
  OUTncl << "  frame(wks)\n";
  OUTncl << "end\n";
  OUTncl.close();
  //
  TripleNCL eTripl{eReal.eProgReal, TargetFile, eFileNC, eFileNCL,
                   ePerm.eBChoice};
  GeneralType eGen(eTripl);
  eCall.SubmitJob(eGen);
}

void DEFINE_PCOLOR_NC(std::string const &eFileNC, GridArray const &GrdArr,
                      MyMatrix<double> const &F_rho, bool const &WriteDEP,
                      std::vector<SeqLineSegment> const &ListLineSegment,
                      std::vector<SingleMarker> const &ListMarker) {
  if (WriteDEP && !GrdArr.GrdArrRho.DEP) {
    std::cerr << "You asked for WriteDEP (related to DrawContourBathy)\n";
    std::cerr << "However, the bathymetry is not available\n";
    throw TerminalException{1};
  }
  //  std::cerr << "eFileNC=" << eFileNC << "\n";
  netCDF::NcFile dataFile(eFileNC, netCDF::NcFile::replace,
                          netCDF::NcFile::nc4);
  int eta = GrdArr.GrdArrRho.LON.rows();
  int xi = GrdArr.GrdArrRho.LON.cols();
  if (GrdArr.GrdArrRho.LAT.rows() != eta || GrdArr.GrdArrRho.LAT.cols() != xi) {
    std::cerr << "Dimension errors 1\n";
    std::cerr << "dim(LON)=" << eta << "/" << xi << "\n";
    std::cerr << "dim(LAT)=" << GrdArr.GrdArrRho.LAT.rows() << "/"
              << GrdArr.GrdArrRho.LAT.cols() << "\n";
    throw TerminalException{1};
  }
  if (F_rho.rows() != eta || F_rho.cols() != xi) {
    std::cerr << "Dimension errors 2\n";
    std::cerr << "dim(LON)=" << eta << "/" << xi << "\n";
    std::cerr << "dim(F  )=" << F_rho.rows() << "/" << F_rho.cols() << "\n";
    throw TerminalException{1};
  }
  if (GrdArr.IsFE == 0) {
    bool ApplyCritValue = true;
    double eCritValue = GetStandardMissingValue();
    double dataMiss[1];
    dataMiss[0] = eCritValue;
    std::string MissVal = "_FillValue";
    //
    netCDF::NcDim eDimEta = dataFile.addDim("eta_rho", eta);
    netCDF::NcDim eDimXi = dataFile.addDim("xi_rho", xi);
    std::vector<std::string> ListDim = {"eta_rho", "xi_rho"};
    netCDF::NcVar eVarLON = dataFile.addVar("lon", "double", ListDim);
    netCDF::NcVar eVarLAT = dataFile.addVar("lat", "double", ListDim);
    //
    netCDF::NcVar eVarF = dataFile.addVar("field", "double", ListDim);
    eVarF.putAtt(MissVal, netCDF::NcType::nc_DOUBLE, 1, dataMiss);
    //
    std::vector<double> valLON(eta * xi), valLAT(eta * xi), valF(eta * xi);
    int idx = 0;
    for (int i = 0; i < eta; i++)
      for (int j = 0; j < xi; j++) {
        valLON[idx] = GrdArr.GrdArrRho.LON(i, j);
        valLAT[idx] = GrdArr.GrdArrRho.LAT(i, j);
        double eValF;
        if (GrdArr.GrdArrRho.MSK(i, j) == 1 || !ApplyCritValue) {
          eValF = F_rho(i, j);
        } else {
          eValF = eCritValue;
        }
        valF[idx] = eValF;
        idx++;
      }
    eVarLON.putVar(valLON.data());
    eVarLAT.putVar(valLAT.data());
    eVarF.putVar(valF.data());
    if (WriteDEP) {
      const MyMatrix<double> &DEP = GetDEP(GrdArr.GrdArrRho);
      if (DEP.size() == 0) {
        std::cerr << "For option WriteDEP, we need to have "
                     "GrdArr.GrdArrRho.DEP assigned\n";
        throw TerminalException{1};
      }
      std::vector<double> valD(eta * xi);
      netCDF::NcVar eVarDEP = dataFile.addVar("dep", "double", ListDim);
      idx = 0;
      for (int i = 0; i < eta; i++)
        for (int j = 0; j < xi; j++) {
          double eValD = DEP(i, j);
          valD[idx] = eValD;
          idx++;
        }
      eVarDEP.putVar(valD.data());
    }
  } else {
    int mnp = eta;
    int mne = GrdArr.INE.rows();
    netCDF::NcDim eDimMnp = dataFile.addDim("mnp", mnp);
    netCDF::NcDim eDimThree = dataFile.addDim("three", 3);
    netCDF::NcDim eDimMne = dataFile.addDim("mne", mne);
    std::vector<std::string> ListDim = {"mnp"};
    netCDF::NcVar eVarLON = dataFile.addVar("lon", "double", ListDim);
    netCDF::NcVar eVarLAT = dataFile.addVar("lat", "double", ListDim);
    netCDF::NcVar eVarF = dataFile.addVar("field", "double", ListDim);
    std::vector<std::string> ListDimINE = {"mne", "three"};
    netCDF::NcVar eVarINE = dataFile.addVar("ele", "int", ListDimINE);
    std::vector<double> valLON(mnp), valLAT(mnp), valF(mnp);
    for (int i = 0; i < mnp; i++) {
      valLON[i] = GrdArr.GrdArrRho.LON(i, 0);
      valLAT[i] = GrdArr.GrdArrRho.LAT(i, 0);
      double eValF = F_rho(i, 0);
      valF[i] = eValF;
    }
    eVarLON.putVar(valLON.data());
    eVarLAT.putVar(valLAT.data());
    eVarF.putVar(valF.data());
    //
    std::vector<int> valI(3 * mne);
    int idx = 0;
    for (int ie = 0; ie < mne; ie++)
      for (int i = 0; i < 3; i++) {
        int eConn = GrdArr.INE(ie, i);
        valI[idx] = eConn;
        idx++;
      }
    eVarINE.putVar(valI.data());
    //
    if (WriteDEP) {
      const MyMatrix<double> &DEP = GetDEP(GrdArr.GrdArrRho);
      netCDF::NcVar eVarDEP = dataFile.addVar("dep", "double", ListDim);
      std::vector<double> valD(mnp);
      for (int i = 0; i < mnp; i++)
        valD[i] = DEP(i, 0);
      eVarDEP.putVar(valD.data());
    }
  }
  if (ListLineSegment.size() > 0)
    ADD_LISTLINESEGMENT_NC(dataFile, ListLineSegment);
  if (ListMarker.size() > 0)
    ADD_LISTMARKER_NC(dataFile, ListMarker);
}

void DEFINE_PCOLOR_NC_GRI(std::string const &eFileNC, GridArray const &GrdArr,
                          MyMatrix<double> const &F_rho) {
  netCDF::NcFile dataFile(eFileNC, netCDF::NcFile::replace,
                          netCDF::NcFile::nc4);
  if (GrdArr.IsFE == 0) {
    std::cerr << "DEFINE_PCOLOR_NC_GRI: The code has not been written yet\n";
    std::cerr << "Case GrdArr.IsFE = 0\n";
    throw TerminalException{1};
  } else {
    int mne = GrdArr.INE.rows();
    int TotDim = 4 * mne;
    netCDF::NcDim eDimTot = dataFile.addDim("tdim", TotDim);
    std::vector<std::string> ListDim = {"tdim"};
    netCDF::NcVar eVarLON = dataFile.addVar("lon", "float", ListDim);
    netCDF::NcVar eVarLAT = dataFile.addVar("lat", "float", ListDim);
    netCDF::NcVar eVarF = dataFile.addVar("field", "float", ListDim);
    std::vector<double> valLON(TotDim), valLAT(TotDim), valF(TotDim);
    int idx = 0;
    for (int ie = 0; ie < mne; ie++) {
      int i1 = GrdArr.INE(ie, 0);
      int i2 = GrdArr.INE(ie, 1);
      int i3 = GrdArr.INE(ie, 2);
      valLON[idx] = static_cast<float>(GrdArr.GrdArrRho.LON(i1, 0));
      valLAT[idx] = static_cast<float>(GrdArr.GrdArrRho.LAT(i1, 0));
      valF[idx] = static_cast<float>(F_rho(i1, 0));
      idx++;
      valLON[idx] = static_cast<float>(GrdArr.GrdArrRho.LON(i2, 0));
      valLAT[idx] = static_cast<float>(GrdArr.GrdArrRho.LAT(i2, 0));
      valF[idx] = static_cast<float>(F_rho(i2, 0));
      idx++;
      valLON[idx] = static_cast<float>(GrdArr.GrdArrRho.LON(i3, 0));
      valLAT[idx] = static_cast<float>(GrdArr.GrdArrRho.LAT(i3, 0));
      valF[idx] = static_cast<float>(F_rho(i3, 0));
      idx++;
      valLON[idx] = static_cast<float>(GrdArr.GrdArrRho.LON(i1, 0));
      valLAT[idx] = static_cast<float>(GrdArr.GrdArrRho.LAT(i1, 0));
      valF[idx] = static_cast<float>(F_rho(i1, 0));
      idx++;
    }
    eVarLON.putVar(valLON.data());
    eVarLAT.putVar(valLAT.data());
    eVarF.putVar(valF.data());
  }
}

void DEFINE_MESH_NC(std::string const &eFileNC, GridArray const &GrdArr) {
  netCDF::NcFile dataFile(eFileNC, netCDF::NcFile::replace,
                          netCDF::NcFile::nc4);
  if (GrdArr.IsFE == 1) {
    int mnp = GrdArr.GrdArrRho.LON.rows();
    int mne = GrdArr.INE.rows();
    MyMatrix<int> ListEdges = GetEdgeSet(GrdArr.INE, mnp);
    int nbEdge = ListEdges.rows();
    netCDF::NcDim eDimMnp = dataFile.addDim("mnp", mnp);
    netCDF::NcDim eDimThree = dataFile.addDim("three", 3);
    netCDF::NcDim eDimTwo = dataFile.addDim("two", 2);
    netCDF::NcDim eDimMne = dataFile.addDim("mne", mne);
    netCDF::NcDim eDimEdge = dataFile.addDim("nbEdges", nbEdge);
    netCDF::NcVar eVarLON =
        dataFile.addVar("lon", "double", std::vector<std::string>({"mnp"}));
    netCDF::NcVar eVarLAT =
        dataFile.addVar("lat", "double", std::vector<std::string>({"mnp"}));
    netCDF::NcVar eVarINE = dataFile.addVar(
        "ele", "int", std::vector<std::string>({"mne", "three"}));
    netCDF::NcVar eVarEDGE = dataFile.addVar(
        "edges", "int", std::vector<std::string>({"nbEdges", "two"}));
    std::vector<double> valLON(mnp), valLAT(mnp);
    for (int i = 0; i < mnp; i++) {
      valLON[i] = GrdArr.GrdArrRho.LON(i, 0);
      valLAT[i] = GrdArr.GrdArrRho.LAT(i, 0);
    }
    eVarLON.putVar(valLON.data());
    eVarLAT.putVar(valLAT.data());
    //
    std::vector<int> valI(3 * mne);
    int idx = 0;
    for (int ie = 0; ie < mne; ie++)
      for (int i = 0; i < 3; i++) {
        int eConn = GrdArr.INE(ie, i);
        valI[idx] = eConn;
        idx++;
      }
    eVarINE.putVar(valI.data());
    //
    std::vector<int> valEDGE(3 * nbEdge);
    idx = 0;
    for (int iedge = 0; iedge < nbEdge; iedge++)
      for (int i = 0; i < 2; i++) {
        int eConn = ListEdges(iedge, i);
        valEDGE[idx] = eConn;
        idx++;
      }
    eVarEDGE.putVar(valEDGE.data());
  } else {
    std::cerr
        << "The corresponding code for finite differences need to be written\n";
    throw TerminalException{1};
  }
}

void DEFINE_MESH_SVG(std::string const &eFileSVG, GridArray const &GrdArr) {
  if (GrdArr.IsFE == 1) {
    int mnp = GrdArr.GrdArrRho.LON.rows();
    MyMatrix<int> ListEdges = GetEdgeSet(GrdArr.INE, mnp);
    int nbEdge = ListEdges.rows();
    //
    double eMult = 10000;
    //    double ShiftLon=180;
    //    double ShiftLat=90;
    double ShiftLon = 0;
    double ShiftLat = 0;
    std::vector<SVGgeneral> ListGeneral;
    for (int iEdge = 0; iEdge < nbEdge; iEdge++) {
      int i1 = ListEdges(iEdge, 0);
      int i2 = ListEdges(iEdge, 1);
      /*
      double x1=180 + GrdArr.GrdArrRho.LON(i1);
      double x2=180 + GrdArr.GrdArrRho.LON(i2);
      double y1=90 - GrdArr.GrdArrRho.LAT(i1);
      double y2=90 - GrdArr.GrdArrRho.LAT(i2);*/
      double x1 = ShiftLon + GrdArr.GrdArrRho.LON(i1) * eMult;
      double x2 = ShiftLon + GrdArr.GrdArrRho.LON(i2) * eMult;
      double y1 = ShiftLat + GrdArr.GrdArrRho.LAT(i1) * eMult;
      double y2 = ShiftLat + GrdArr.GrdArrRho.LAT(i2) * eMult;
      coor ePt{x1, y1};
      coor fPt{x2, y2};
      std::vector<int> color{255, 0, 0};
      double Size = 2;
      SVGline eLine{ePt, fPt, {color, Size, ""}};
      ListGeneral.push_back(SVGgeneral(eLine));
    }
    SVGplotDescription eSVGplot;
    eSVGplot.FrameOption = 1;
    eSVGplot.height = 30;
    eSVGplot.width = 30;
    eSVGplot.RoundMethod = 2;
    eSVGplot.ListGeneral = ListGeneral;
    GeneralWriteSVGfile(eFileSVG, eSVGplot);
  } else {
    std::cerr
        << "The corresponding code for finite differences need to be written\n";
    throw TerminalException{1};
  }
}

void PLOT_MESH(DrawArr const &eDrawArr, GridArray const &GrdArr,
               NCLcaller<GeneralType> &eCall,
               PermanentInfoDrawing const &ePerm) {
  RealInfoNclExtension eReal = GetRealInfoNclExtension(ePerm.Extension);
  std::string FileName = ePerm.PicPrefix + "mesh";
  std::string TargetFile = FileName + "_storsave." + eReal.eExtensionReal;
  std::cerr << "Beginning of PLOT_MESH\n";
  int IsFE = GrdArr.IsFE;
  if (IsFE == 0)
    return;
  bool InPlaceRun = ePerm.eBChoice.InPlaceRun;
  std::string eFileNC =
      FinalFile(InPlaceRun, TargetFile, ePerm.PrefixTemp.str() + "mesh.nc");
  std::string eFileNCL =
      FinalFile(InPlaceRun, TargetFile, ePerm.PrefixTemp.str() + "mesh.ncl");
  //  std::cerr << "eFileNC=" << eFileNC << "\n";
  //  std::cerr << "eFileNCL=" << eFileNCL << "\n";
  DEFINE_MESH_NC(eFileNC, GrdArr);
  std::ofstream OUTncl(eFileNCL);
  OUTncl << "load \"$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl\"\n";
  OUTncl << "load \"$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl\"\n";
  OUTncl << "begin\n";
  OUTncl << "  ;\n";
  OUTncl << "  ; Data reading\n";
  OUTncl << "  ;\n";
  OUTncl << "  f = addfile(\"" << FinalFileInScript(InPlaceRun, eFileNC)
         << "\", \"r\")\n";
  OUTncl << "  lat  = f->lat\n";
  OUTncl << "  lon  = f->lon\n";
  OUTncl << "  edges  = f->edges\n";
  OUTncl << "  ms=dimsizes(lat)\n";
  OUTncl << "  mnp=ms(0)\n";
  OUTncl << "  eVarF=new( (/mnp/), float)\n";
  OUTncl << "  eVarF=0\n";
  OUTncl << "  wks  = gsn_open_wks (\"" << eReal.eExtensionReal << "\",\""
         << FinalFileInScript(InPlaceRun, FILE_RemoveExtension(TargetFile))
         << "\")\n";
  OUTncl << "  res2 = True               ; plot mods desired\n";
  OUTncl << "  res2@gsnDraw   = False\n";
  OUTncl << "  res2@gsnFrame  = False\n";
  OUTncl << "  res2@gsnMaximize     = True    ; Maximize plot in frame\n";
  //  OUTncl << "  res2@gsnPaperOrientation  = \"Portrait\"\n";
  OUTncl << "  res2@gsnPaperOrientation = \"" << eDrawArr.LandPortr << "\"\n";
  OUTncl << "  ;\n";
  OUTncl << "  ; General frame information\n";
  OUTncl << "  ;\n";
  OUTncl << "  res2@mpProjection = \"Mercator\"\n";
  OUTncl << "  res2@mpLimitMode         = \"Corners\"             ; choose "
            "range of map\n";
  OUTncl << "  res2@mpLeftCornerLatF    = " << eDrawArr.eQuadFrame.MinLat
         << "\n";
  OUTncl << "  res2@mpLeftCornerLonF    = " << eDrawArr.eQuadFrame.MinLon
         << "\n";
  OUTncl << "  res2@mpRightCornerLatF   = " << eDrawArr.eQuadFrame.MaxLat
         << "\n";
  OUTncl << "  res2@mpRightCornerLonF   = " << eDrawArr.eQuadFrame.MaxLon
         << "\n";
  OUTncl << "  res2@pmTickMarkDisplayMode  = \"Always\"           ; turn on "
            "tickmarks\n";
  OUTncl << "  res2@mpFillOn      = False\n";
  OUTncl << "  res2@sfXArray            = lon\n";
  OUTncl << "  res2@sfYArray            = lat\n";
  if (eDrawArr.UseNativeGrid) {
    OUTncl << "  res2@sfElementNodes      = f->ele\n";
    OUTncl << "  res2@sfFirstNodeIndex    = 0\n";
  }
  OUTncl << "  map = gsn_csm_contour_map(wks,eVarF,res2)\n";
  OUTncl << "  res = True\n";
  OUTncl << "  res@gsLineThicknessF = 1.5\n";
  OUTncl << "  res@gsLineColor  = \"dodgerblue1\"\n";
  OUTncl << "  ns=dimsizes(edges)\n";
  OUTncl << "  lines = new(ns(0),graphic)   ; array to hold polylines\n";
  OUTncl << "  do iedge=0,ns(0)-1\n";
  OUTncl << "    i1=edges(iedge,0)\n";
  OUTncl << "    i2=edges(iedge,1)\n";
  OUTncl << "    xp=(/lon(i1), lon(i2)/)\n";
  OUTncl << "    yp=(/lat(i1), lat(i2)/)\n";
  OUTncl << "    lines(iedge)=gsn_add_polyline(wks,map,xp,yp,res)\n";
  OUTncl << "    delete(xp)\n";
  OUTncl << "    delete(yp)\n";
  OUTncl << "  end do\n";
  OUTncl << "  plot = gsn_csm_contour_map(wks,eVarF,res2)\n";
  OUTncl << "  draw(plot)\n";
  OUTncl << "  frame(wks)\n";
  OUTncl << "end\n";
  OUTncl.close();
  //
  TripleNCL eTripl{eReal.eProgReal, TargetFile, eFileNC, eFileNCL,
                   ePerm.eBChoice};
  GeneralType eGen(eTripl);
  eCall.SubmitJob(eGen);
}

void PLOT_PCOLOR_BASEMAP(std::string const &FileName, GridArray const &GrdArr,
                         DrawArr const &eDrawArr, RecVar const &eRecVar,
                         NCLcaller<GeneralType> &eCall,
                         PermanentInfoDrawing const &ePerm) {
  std::string TargetFile = FileName + "_storsave." + ePerm.Extension;
  RecSymbolic RecS = eRecVar.RecS;
  //  std::cerr << "STRALL RecS.strAll=" << RecS.strAll << "\n";
  bool InPlaceRun = ePerm.eBChoice.InPlaceRun;
  std::string eFileNC =
      FinalFile(InPlaceRun, TargetFile,
                ePerm.PrefixTemp.str() + "DataPcolor_" + eDrawArr.VarNameUF +
                    "_" + RecS.strAll + ".nc");
  std::string eFilePY =
      FinalFile(InPlaceRun, TargetFile,
                ePerm.PrefixTemp.str() + "ScriptPcolor_" + eDrawArr.VarNameUF +
                    "_" + RecS.strAll + ".py");
  //  std::cerr << "eFilePY=" << eFilePY << "\n";
  //  std::cerr << "min/max(eRecVar.F)=" << eRecVar.F.minCoeff() << " / " <<
  //  eRecVar.F.maxCoeff() << "\n";
  DEFINE_PCOLOR_NC(eFileNC, GrdArr, eRecVar.F, eDrawArr.DrawContourBathy,
                   eDrawArr.ListLineSegment, eDrawArr.ListMarker);
  bool IsSpherical = GrdArr.IsSpherical;
  if (!IsSpherical) {
    std::cerr << "Error in PYTHON_PCOLOR_BASEMAP\n";
    std::cerr << "This python code works only in spherical coordinate\n";
    throw TerminalException{1};
  }
  int IsFE = GrdArr.IsFE;
  if (IsFE == 1) {
    std::cerr << "Error in PYTHON_PCOLOR_BASEMAP\n";
    std::cerr << "This python code works only for finite difference\n";
    throw TerminalException{1};
  }
  //
  std::ofstream os(eFilePY);
  os << "from mpl_toolkits.basemap import Basemap\n";
  os << "import os\n";
  os << "import matplotlib as mpl\n";
  os << "mpl.use('Agg')\n";
  os << "import matplotlib.pyplot as plt\n";
  os << "import numpy as np\n";
  os << "from collections import namedtuple\n";
  os << "import math\n";
  os << "import netCDF4\n";
  os << "\n";
  os << "InFile = \'" << FinalFileInScript(InPlaceRun, eFileNC) + "\'\n";
  os << "\n";
  os << "dataset = netCDF4.Dataset(InFile)\n";
  os << "lon = dataset.variables[\'lon\'][:,:]\n";
  os << "lat = dataset.variables[\'lat\'][:,:]\n";
  os << "z = dataset.variables['field'][:,:]\n";
  os << "\n";
  double minLON = eDrawArr.eQuadFrame.MinLon;
  double maxLON = eDrawArr.eQuadFrame.MaxLon;
  double midLON = (minLON + maxLON) / static_cast<double>(2);
  double minLAT = eDrawArr.eQuadFrame.MinLat;
  double maxLAT = eDrawArr.eQuadFrame.MaxLat;
  double midLAT = (minLAT + maxLAT) / static_cast<double>(2);
  std::cerr << "eDrawArr.eQuadFrame=" << eDrawArr.eQuadFrame.MinLon << " / "
            << eDrawArr.eQuadFrame.MaxLon << " / " << eDrawArr.eQuadFrame.MinLat
            << " / " << eDrawArr.eQuadFrame.MaxLat << "\n";
  std::string strResolution;
  if (eDrawArr.GridResolution == "HighRes")
    strResolution = "h";
  if (eDrawArr.GridResolution == "MediumRes")
    strResolution = "i";
  if (eDrawArr.GridResolution == "LowRes")
    strResolution = "l";

  os << "lon_0 = " << midLON << "\n";
  os << "lat_0 = " << midLAT << "\n";
  os << "minLAT = " << minLAT << "\n";
  os << "maxLAT = " << maxLAT << "\n";
  os << "minLON = " << minLON << "\n";
  os << "maxLON = " << maxLON << "\n";
  os << "my_dpi=96\n";
  os << "plt.figure(figsize=(1600/my_dpi, 1000/my_dpi), dpi=my_dpi)\n";
  os << "map = "
        "Basemap(projection='merc',lon_0=lon_0,lat_0=lat_0,lat_ts=lat_0,\\\n";
  os << "          llcrnrlat=minLAT,urcrnrlat=maxLAT,\\\n";
  os << "          llcrnrlon=minLON,urcrnrlon=maxLON,\\\n";
  os << "          rsphere=6371200.,resolution=\'" << strResolution
     << "\',area_thresh=1000)\n";
  os << "print \"After map\"";
  os << "\n";
  os << "map.drawcoastlines(linewidth=0.25)\n";
  os << "map.drawcountries(linewidth=0.25)\n";
  os << "map.fillcontinents(color=\'coral\',lake_color=\'aqua\')\n";
  os << "map.drawmapboundary(fill_color=\'aqua\')\n";
  os << "\n";
  os << "x, y = map(lon, lat)\n";
  double eMinVal = RecS.minval;
  double eMaxVal = RecS.maxval;
  int nbLevelSpa = eDrawArr.nbLevelSpa;
  double deltaVal = (eMaxVal - eMinVal) / static_cast<double>(nbLevelSpa - 1);
  std::vector<double> clevs(nbLevelSpa);
  for (int iLev = 0; iLev < nbLevelSpa; iLev++) {
    double eLev = eMinVal + iLev * deltaVal;
    clevs[iLev] = eLev;
  }
  os << "clevs = [";
  for (int iLev = 0; iLev < nbLevelSpa; iLev++) {
    if (iLev > 0)
      os << ", ";
    os << clevs[iLev];
  }
  os << "]\n";
  //  os << "cs = map.contour(x,y, z, 15,linewidths=1.5)\n";   // very old code
  //  os << "cs = map.contour(x, y, z, clevs, linewidths=1.5)\n";
  os << "cs = map.pcolormesh(x, y, z, linewidths=1.5)\n";
  os << "\n";
  os << "cbar = plt.colorbar(cs)\n";
  os << "cbar.ax.set_ylabel(\'" << RecS.Unit << "\')\n";
  if (eDrawArr.ListLineSegment.size() > 0) {
    os << "LLS_ListLon = dataset.variables[\'LLS_ListLon\'][:]\n";
    os << "LLS_ListLat = dataset.variables[\'LLS_ListLat\'][:]\n";
    os << "nbLine=len(LLS_ListLon)/2\n";
    os << "for iLine in range(nbLine):\n";
    os << "    lon1=LLS_ListLon[2*iLine]\n";
    os << "    lon2=LLS_ListLon[2*iLine+1]\n";
    os << "    lat1=LLS_ListLat[2*iLine]\n";
    os << "    lat2=LLS_ListLat[2*iLine+1]\n";
    os << "    "
          "map.drawgreatcircle(lon1,lat1,lon2,lat2,linewidth=2,color=\'b\')\n";
  }
  os << "\n";
  os << "plt.title(\'" << eDrawArr.TitleStr << "\')\n";
  std::string OutFile = FileName + "_storsave." + ePerm.Extension;
  os << "plt.savefig(\'" << FinalFileInScript(InPlaceRun, OutFile) << "\')\n";
  os << "print \"After savefig\"\n";
  os.close();
  //
  TripleNCL eTripl{"python", TargetFile, eFileNC, eFilePY, ePerm.eBChoice};
  GeneralType eGen(eTripl);
  eCall.SubmitJob(eGen);
}

void PLOT_PCOLOR_NCL(std::string const &FileName, GridArray const &GrdArr,
                     DrawArr const &eDrawArr, RecVar const &eRecVar,
                     NCLcaller<GeneralType> &eCall,
                     PermanentInfoDrawing const &ePerm) {
  //  std::cerr << "Beginning of PLOT_PCOLOR_NCL\n";
  RealInfoNclExtension eReal = GetRealInfoNclExtension(ePerm.Extension);
  std::string TargetFile = FileName + "_storsave." + eReal.eExtensionReal;
  RecSymbolic RecS = eRecVar.RecS;
  bool InPlaceRun = ePerm.eBChoice.InPlaceRun;
  //  std::cerr << "TargetFile=" << TargetFile << "\n";
  // td::cerr << "PrefixTemp.str()=" << ePerm.PrefixTemp.str() << "\n";
  std::string eFileNC =
      FinalFile(InPlaceRun, TargetFile,
                ePerm.PrefixTemp.str() + "DataPcolor_" + eDrawArr.VarNameUF +
                    "_" + RecS.strAll + ".nc");
  // std::cerr << "eFileNC=" << eFileNC << "\n";
  std::string eFileNCL =
      FinalFile(InPlaceRun, TargetFile,
                ePerm.PrefixTemp.str() + "ScriptPcolor_" + eDrawArr.VarNameUF +
                    "_" + RecS.strAll + ".ncl");
  //  std::cerr << "eFileNCL=" << eFileNCL << "\n";
  DEFINE_PCOLOR_NC(eFileNC, GrdArr, eRecVar.F, eDrawArr.DrawContourBathy,
                   eDrawArr.ListLineSegment, eDrawArr.ListMarker);
  int IsFE = GrdArr.IsFE;
  bool IsSpherical = GrdArr.IsSpherical;
  //
  std::ofstream OUTncl(eFileNCL);
  OUTncl << "load \"$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl\"\n";
  OUTncl << "load \"$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl\"\n";
  OUTncl << "begin\n";
  OUTncl << "  ;\n";
  OUTncl << "  ; Data reading\n";
  OUTncl << "  ;\n";
  OUTncl << "  f = addfile(\"" << FinalFileInScript(InPlaceRun, eFileNC)
         << "\", \"r\")\n";
  OUTncl << "  lat  = f->lat\n";
  OUTncl << "  lon  = f->lon\n";
  OUTncl << "  eVarF  = f->field\n";
  OUTncl << "  wks  = gsn_open_wks (\"" << eReal.eExtensionReal << "\",\""
         << FinalFileInScript(InPlaceRun, FILE_RemoveExtension(TargetFile))
         << "\")\n";
  if (eDrawArr.DrawContourBathy) {
    OUTncl << "  DEP  = f->dep\n";
    OUTncl << "  res1 = True\n";
    OUTncl << "  res1@gsnTickMarksOn   = False; no tickmarks\n";
    OUTncl << "  res1@gsnDraw          = False; don't draw\n";
    OUTncl << "  res1@gsnFrame         = False; don't advance frame\n";
    OUTncl << "  res1@gsnLeftString    = \"\"; no titles\n";
    OUTncl << "  res1@gsnRightString   = \"\"\n";
    OUTncl << "  res1@tiXAxisString    = \"\"\n";
    OUTncl << "  res1@tiYAxisString    = \"\"\n";
    OUTncl << "  res1@cnLineThicknessF = 1.5; thicker contours\n";
    OUTncl << "  res1@cnLineLabelsOn   = False; no line labels\n";
    OUTncl << "  plot2 = gsn_csm_contour(wks,DEP,res1)\n";
  }
  OUTncl << "  res2 = True               ; plot mods desired\n";
  OUTncl << "  res2@gsnDraw   = False\n";
  OUTncl << "  res2@gsnFrame  = False\n";
  OUTncl << "  res2@gsnMaximize     = True    ; Maximize plot in frame\n";
  //  OUTncl << "  res2@gsnPaperOrientation  = \"Portrait\"\n";
  OUTncl << "  res2@gsnPaperOrientation = \"" << eDrawArr.LandPortr << "\"\n";
  OUTncl << "  ;\n";
  OUTncl << "  ; General frame information\n";
  OUTncl << "  ;\n";
  if (IsSpherical) {
    OUTncl << "  res2@mpProjection = \"Mercator\"\n";
    OUTncl << "  res2@mpLimitMode         = \"Corners\"             ; choose "
              "range of map\n";
    OUTncl << "  res2@mpLeftCornerLatF    = " << eDrawArr.eQuadFrame.MinLat
           << "\n";
    OUTncl << "  res2@mpLeftCornerLonF    = " << eDrawArr.eQuadFrame.MinLon
           << "\n";
    OUTncl << "  res2@mpRightCornerLatF   = " << eDrawArr.eQuadFrame.MaxLat
           << "\n";
    OUTncl << "  res2@mpRightCornerLonF   = " << eDrawArr.eQuadFrame.MaxLon
           << "\n";
    if (eDrawArr.FillLand) {
      OUTncl << "  res2@mpFillOn      = True\n";
      OUTncl << "  res2@mpDataBaseVersion      = \"" << eDrawArr.GridResolution
             << "\"          ; use high resolution coast\n";
      OUTncl << "  res2@mpLandFillColor       = \"gray\"            ; set land "
                "to be gray\n";
    } else {
      OUTncl << "  res2@mpFillOn      = False\n";
    }
  } else {
    OUTncl << "  res2@sfXArray = lon\n";
    OUTncl << "  res2@sfYArray = lat\n";
    OUTncl << "  res2@trGridType = \"TriangularMesh\"\n";
  }
  OUTncl << "  res2@pmTickMarkDisplayMode  = \"Always\"           ; turn on "
            "tickmarks\n";
  OUTncl << "  ;\n";
  OUTncl << "  ; Contour map information\n";
  OUTncl << "  ;\n";
  OUTncl << "  res2@cnFillDrawOrder        = \"PreDraw\"\n";
  OUTncl << "  res2@cnFillOn             = " << NCL_bool(eDrawArr.cnFillOn)
         << "   ; turn on color for contours\n";
  OUTncl << "  res2@cnLinesOn            = " << NCL_bool(eDrawArr.cnLinesOn)
         << "\n";
  OUTncl << "  res2@cnLineLabelsOn       = "
         << NCL_bool(eDrawArr.cnLineLabelsOn)
         << "   ; turn off contour line labels\n";
  if (IsFE == 1 && eDrawArr.cnFillMode == "CellFill") {
    std::cerr << "The \"CellFill\" option can only be used in finite "
                 "difference grids\n";
    throw TerminalException{1};
  }
  if (IsFE == 1 && eDrawArr.cnFillMode == "AreaFill") {
    OUTncl << "; This is an option to emulate \"AreaFill\" method\n";
    OUTncl << "  res2@cnFillMode           = \"RasterFill\"\n";
    OUTncl << "  res2@cnRasterSmoothingOn = True\n";
  } else {
    OUTncl << "  res2@cnFillMode           = \"" << eDrawArr.cnFillMode
           << "\"\n";
  }
  OUTncl
      << "            ; AreaFill : slow and buggy but maybe more beautiful\n";
  OUTncl << "            ; RasterFill : fast and efficient\n";
  OUTncl << "            ; CellFill : similar to RasterFill but only for "
            "finite difference\n";
  OUTncl << "  ;  res2@cnRasterSmoothingOn  = True\n";
  OUTncl << "  res2@cnSmoothingOn = " << NCL_bool(eDrawArr.cnSmoothingOn)
         << "\n";
  OUTncl << "  ;  res2@cnSmoothingDistanceF  = 0.05\n";
  OUTncl << "  res2@cnSmoothingTensionF  = -1\n";
  OUTncl << "  res2@cnLevelSelectionMode = \"ManualLevels\"     ; set manual "
            "contour levels\n";
  OUTncl << "  res2@cnMinLevelValF       = " << RecS.minval
         << "  ; set min contour level\n";
  OUTncl << "  res2@cnMaxLevelValF       = " << RecS.maxval
         << "  ; set max contour level\n";
  //  std::cerr << "RecS(minval/maxval)=" << RecS.minval << " / " << RecS.maxval
  //  << "\n";
  int nbLevelSpa = eDrawArr.nbLevelSpa;
  //  std::cerr << "nbLevelSpa=" << nbLevelSpa << "\n";
  double TheLevelSpa =
      (RecS.maxval - RecS.minval) / static_cast<double>(nbLevelSpa);
  OUTncl << "  res2@cnLevelSpacingF      = " << TheLevelSpa
         << "     ; set contour spacing\n";
  int nbLabelStride = eDrawArr.nbLabelStride;
  OUTncl << "  res2@lbLabelStride            = " << nbLabelStride << "\n";
  OUTncl << "  ;  res2@gsnScalarContour     = False               ; contours "
            "desired\n";
  if (eDrawArr.DoTitle) {
    OUTncl << "  res2@tiMainString    = \"" << eDrawArr.TitleStr << "\"\n";
    OUTncl << "  res2@tiMainFont      = \"Helvetica\"\n";
    OUTncl << "  res2@tiMainFontHeightF = 0.015\n";
    OUTncl << "  ;  res2@cnTitlePosition  = \"Top\"\n";
  }
  OUTncl << "  res2@gsnSpreadColors      = True               ; use full color "
            "map\n";
  OUTncl << "  res2@gsnSpreadColorEnd     = -3\n";
  // LABEL BAR
  OUTncl << "  ;\n";
  OUTncl << "  ; Label bar plotting\n";
  OUTncl << "  ;\n";
  if (eDrawArr.DoColorBar && IsSpherical) {
    OUTncl << "  res2@lbLabelBarOn = True\n";
  } else {
    OUTncl << "  res2@lbLabelBarOn = False\n";
  }
  if (eDrawArr.DoTitleString) {
    OUTncl << "  res2@lbTitleString    = \"" << RecS.VarName1 << " ["
           << RecS.Unit << "]\"\n";
  } else {
    OUTncl << "  res2@lbTitleString    = \"\"  \n";
  }
  OUTncl << "  res2@lbTitleFont      = \"Helvetica\"\n";
  OUTncl << "  res2@lbTitleFontHeightF = 0.015\n";
  OUTncl << "  res2@lbTitleDirection     = \"Across\"\n";
  OUTncl << "  res2@lbTitlePosition = \"Right\"\n";
  OUTncl << "  res2@lbTitleAngleF = 90\n";
  OUTncl << "  res2@lbOrientation        = \"Vertical\"     ; Vertical label "
            "bar\n";
  OUTncl << "  res2@pmLabelBarOrthogonalPosF = 0.025          ; move label bar "
            "closer\n";
  OUTncl << "  ;  res2@lbHeightF               = 0.7          ; move label bar "
            "closer\n";
  OUTncl << "  res2@pmLabelBarDisplayMode = \"Always\"          ; Turn on a "
            "label bar.\n";
  OUTncl << "  res2@lbPerimOn             = False             ; no box around "
            "it\n";
  OUTncl << "  res2@lbBoxLinesOn         = False               ; Yes/No "
            "labelbar box lines\n";
  if (IsSpherical) {
    OUTncl << "  res2@pmLabelBarWidthF = 0.03\n";
  }
  OUTncl << "  ; res2@gsnRightString  = \"Sea surface elevation\"\n";
  OUTncl << "  ; res2@gsnLeftString    = \"Difference\"\n";
  // COLORMAP
  OUTncl << "  ;\n";
  OUTncl << "  ; Colormap assignation\n";
  OUTncl << "  ;\n";
  OUTncl << "  gsn_define_colormap (wks,\"" << eDrawArr.ColorMap << "\")\n";
  OUTncl << "  ;     other possibilities: hotres, rainbow, ViBlGrWhYeOrRe, "
            "BlWhRe, GrayWhiteGray, BlGrYeOrReVi200\n";
  OUTncl << "  i = NhlNewColor(wks,0.8,0.8,0.8)      ; add gray to colormap\n";
  OUTncl << "  i = NhlNewColor(wks,0.9,0.9,0.9)      ; add gray to colormap\n";
  // THE PLOTTING ITSELF
  OUTncl << "  ;\n";
  OUTncl << "  ; Pcolor kind of plot\n";
  OUTncl << "  ;\n";
  if (IsFE == 1) {
    OUTncl << "  res2@sfXArray            = lon\n";
    OUTncl << "  res2@sfYArray            = lat\n";
    if (eDrawArr.UseNativeGrid) {
      OUTncl << "  res2@sfElementNodes      = f->ele\n";
      OUTncl << "  res2@sfFirstNodeIndex    = 0\n";
    }
  } else {
    OUTncl << "  eVarF@lat2d=lat\n";
    OUTncl << "  eVarF@lon2d=lon\n";
  }
  if (IsSpherical) {
    OUTncl << "  plot = gsn_csm_contour_map(wks,eVarF,res2)\n";
  } else {
    OUTncl << "  plot = gsn_csm_contour(wks,eVarF,res2)\n";
  }
  ADD_LISTLINESEGMENT(OUTncl, eDrawArr);
  ADD_LISTMARKER(OUTncl, eDrawArr);
  ADD_RIVER(OUTncl, eDrawArr);
  ADD_ANNOTATION_TEXT(OUTncl, eDrawArr.TheAnnot);
  for (auto &eLine : eDrawArr.ListInsertLines)
    OUTncl << eLine << "\n";
  OUTncl << "  draw(plot)\n";
  OUTncl << "  frame(wks)\n";
  OUTncl << "end\n";
  OUTncl.close();
  //
  TripleNCL eTripl{eReal.eProgReal, TargetFile, eFileNC, eFileNCL,
                   ePerm.eBChoice};
  GeneralType eGen(eTripl);
  eCall.SubmitJob(eGen);
}

void PLOT_PCOLOR(std::string const &FileName, GridArray const &GrdArr,
                 DrawArr const &eDrawArr, RecVar const &eRecVar,
                 NCLcaller<GeneralType> &eCall,
                 PermanentInfoDrawing const &ePerm) {
  if (!IsEqualSizeMatrices(GrdArr.GrdArrRho.LON, eRecVar.F)) {
    std::cerr
        << "Now GrdArr.GrdArrRho.LON and eRecVar.F should have the same size\n";
    std::cerr << "|GrdArr.GrdArrRho.LON|=" << GrdArr.GrdArrRho.LON.rows()
              << " / " << GrdArr.GrdArrRho.LON.cols() << "\n";
    std::cerr << "           |eRecVar.F|=" << eRecVar.F.rows() << " / "
              << eRecVar.F.cols() << "\n";
    throw TerminalException{1};
  }
  std::string eProg = ePerm.ChoiceProgram.at("Pcolor_method");
  //  std::cerr << "PLOT_PCOLOR eProg=" << eProg << " eRecVar.RecS.strAll=" <<
  //  eRecVar.RecS.strAll << "\n";
  bool IsMatch = false;
  if (eProg == "ncl") {
    PLOT_PCOLOR_NCL(FileName, GrdArr, eDrawArr, eRecVar, eCall, ePerm);
    IsMatch = true;
  }
  if (eProg == "python") {
    PLOT_PCOLOR_BASEMAP(FileName, GrdArr, eDrawArr, eRecVar, eCall, ePerm);
    IsMatch = true;
  }
  if (!IsMatch) {
    std::cerr << "eProg = " << eProg
              << " while allowed values are python and ncl\n";
    std::cerr << "If you want the crop, use pngcrop or pdfcrop in the type "
                 "instead of png, pdf\n";
    throw TerminalException{1};
  }
}

void PLOT_PCOLOR_OR_QUIVER(std::string const &FileName, GridArray const &GrdArr,
                           DrawArr const &eDrawArr, RecVar const &eRecVar,
                           NCLcaller<GeneralType> &eCall,
                           PermanentInfoDrawing const &ePerm) {
  std::string eNature = eRecVar.RecS.VarNature;
  if (eNature == "uv") {
    PLOT_QUIVER(FileName, GrdArr, eDrawArr, eRecVar, eCall, ePerm);
    return;
  }
  if (eNature == "rho") {
    PLOT_PCOLOR(FileName, GrdArr, eDrawArr, eRecVar, eCall, ePerm);
    return;
  }
  std::cerr << "Allowed eNature are uv and rho. Here eNature=" << eNature
            << "\n";
  throw TerminalException{1};
}

struct DrawLinesArr {
  bool DoTitle;
  std::string TitleStr;
  std::string VarName;
  AnnotationRec TheAnnot;
  // Drawing itself
  bool DoTitleLines;
  bool DoExplicitLabel;
  int nbLabel;
  std::string StyleDate;
  double TheMax;
  double TheMin;
  bool IsTimeSeries;
  bool PairComparison;
  bool DrawHorizVertLines;
  std::string XAxisString;
  std::string YAxisString;
  std::vector<std::string> ListName_plot;
  MyVector<double> ListX;
  std::vector<MyVector<double>> ListListVect;
};

void LINES_DEFINE_NC(std::string const &eFileNC, DrawLinesArr const &eDrawArr) {
  double eCritValue = GetStandardMissingValue();
  double dataMiss[1];
  dataMiss[0] = eCritValue;
  std::string MissVal = "_FillValue";
  //
  netCDF::NcFile dataFile(eFileNC, netCDF::NcFile::replace,
                          netCDF::NcFile::nc4);
  int nbArr = eDrawArr.ListListVect.size();
  int nbEntry = eDrawArr.ListListVect[0].size();
  netCDF::NcDim eDimArr = dataFile.addDim("nbArr", nbArr);
  netCDF::NcDim eDimEntry = dataFile.addDim("nbEntry", nbEntry);
  netCDF::NcDim eDimSix = dataFile.addDim("six", 6);
  std::vector<std::string> ListDim = {"nbArr", "nbEntry"};
  std::vector<std::string> ListDimX = {"nbEntry"};
  std::vector<std::string> ListDimSix = {"nbEntry", "six"};
  netCDF::NcVar eVar = dataFile.addVar("ListListVect", "double", ListDim);
  eVar.putAtt(MissVal, netCDF::NcType::nc_DOUBLE, 1, dataMiss);
  netCDF::NcVar eVarX = dataFile.addVar("ListX", "double", ListDimX);
  netCDF::NcVar eVarSix = dataFile.addVar("ListTimeSix", "int", ListDimSix);
  //
  std::vector<double> val(nbArr * nbEntry);
  int idx = 0;
  for (int iArr = 0; iArr < nbArr; iArr++) {
    MyVector<double> eVect = eDrawArr.ListListVect[iArr];
    for (int i = 0; i < nbEntry; i++) {
      double eVal = eVect(i);
      if (eVal < -9900)
        eVal = eCritValue;
      val[idx] = eVal;
      idx++;
    }
  }
  eVar.putVar(val.data());
  //
  std::vector<double> valX(nbEntry);
  for (int i = 0; i < nbEntry; i++)
    valX[i] = eDrawArr.ListX(i);
  eVarX.putVar(valX.data());
  //
  if (eDrawArr.IsTimeSeries) {
    std::vector<double> valYear(nbEntry);
    for (int i = 0; i < nbEntry; i++) {
      double eMJD = eDrawArr.ListX(i);
      std::vector<int> eVect = DATE_ConvertMjd2six(eMJD);
      double FirstDayYear = DATE_ConvertSix2mjd({eVect[0], 1, 1, 0, 0, 0});
      double FirstDayYearNext =
          DATE_ConvertSix2mjd({eVect[0] + 1, 1, 1, 0, 0, 0});
      double lenYear = FirstDayYearNext - FirstDayYear;
      double eX =
          static_cast<double>(eVect[0]) + (eMJD - FirstDayYear) / lenYear;
      valYear[i] = eX;
    }
    netCDF::NcVar eVarTime = dataFile.addVar("ListTime", "double", ListDimX);
    eVarTime.putVar(valYear.data());
    //
    //    netCDF::NcDim eDimSix=dataFile.addDim("six", 6);
    //    std::vector<std::string> ListDimSix={"nbEntry", "six"};
    //    netCDF::NcVar eVarSix=dataFile.addVar("ListTimeSix", "integer",
    //    ListDimSix);
    std::vector<int> valSix(nbEntry * 6);
    int idx = 0;
    for (int i = 0; i < nbEntry; i++) {
      double eMJD = eDrawArr.ListX(i);
      std::vector<int> eVect = DATE_ConvertMjd2six(eMJD);
      for (int u = 0; u < 6; u++) {
        valSix[idx] = eVect[u];
        idx++;
      }
    }
    eVarSix.putVar(valSix.data());
  }
}

void LINES_PLOT_NCL(std::string const &FileName, DrawLinesArr const &eDrawArr,
                    NCLcaller<GeneralType> &eCall,
                    PermanentInfoDrawing const &ePerm) {
  RealInfoNclExtension eReal = GetRealInfoNclExtension(ePerm.Extension);
  std::string TargetFile = FileName + "_storsave." + eReal.eExtensionReal;
  bool InPlaceRun = ePerm.eBChoice.InPlaceRun;
  std::string eFileNC = FinalFile(InPlaceRun, TargetFile,
                                  ePerm.PrefixTemp.str() + "DataLines_" +
                                      eDrawArr.VarName + ".nc");
  std::string eFileNCL = FinalFile(InPlaceRun, TargetFile,
                                   ePerm.PrefixTemp.str() + "ScriptLines_" +
                                       eDrawArr.VarName + ".ncl");
  LINES_DEFINE_NC(eFileNC, eDrawArr);
  int nbArr = eDrawArr.ListListVect.size();
  //
  std::ofstream OUTncl(eFileNCL);
  OUTncl << "load \"$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl\"\n";
  OUTncl << "load \"$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl\"\n";
  PrintMyScriptSubtitle(OUTncl);
  OUTncl << ";************************************************\n";
  OUTncl << "begin\n";
  OUTncl << "  ;\n";
  OUTncl << "  ; Data reading\n";
  OUTncl << "  ;\n";
  OUTncl << "  f = addfile(\"" << FinalFileInScript(InPlaceRun, eFileNC)
         << "\", \"r\")\n";
  if (eDrawArr.IsTimeSeries && !eDrawArr.DoExplicitLabel) {
    OUTncl << "  ListX=f->ListTime\n";
  } else {
    OUTncl << "  ListX=f->ListX\n";
  }
  OUTncl << "  Data=f->ListListVect\n";
  OUTncl << "  TheMax=max(Data)\n";
  OUTncl << "  TheMin=0\n";
  OUTncl << "  wks  = gsn_open_wks (\"" << eReal.eExtensionReal << "\",\""
         << FinalFileInScript(InPlaceRun, FILE_RemoveExtension(TargetFile))
         << "\")\n";
  OUTncl << "  res                  = True   ; plot mods desired\n";
  OUTncl << "  res@gsnFrame          = False\n";
  OUTncl << "  res@xyMonoDashPattern = False\n";
  OUTncl << "  res@xyDashPatterns = (/0/)\n";
  OUTncl << "  res@tiMainFontHeightF  = 0.015\n";
  if (eDrawArr.DrawHorizVertLines) {
    OUTncl << "  res@tmXMajorGrid = True\n";
    OUTncl << "  res@tmXMajorGridThicknessF      = 1.0           ; 2.0 is "
              "default\n";
    OUTncl << "  res@tmXMajorGridLineDashPattern = 2             ; select "
              "short dash lines\n";
    OUTncl << "  res@tmYMajorGrid                = True          ; implement y "
              "grid\n";
    OUTncl << "  res@tmYMajorGridThicknessF      = 1.0           ; 2.0 is "
              "default\n";
    OUTncl << "  res@tmYMajorGridLineDashPattern = 2             ; select "
              "short dash lines\n";
  }
  int nbCharX = eDrawArr.XAxisString.size();
  if (nbCharX > 0) {
    OUTncl << "  res@tiXAxisString   = \"x (m)\"\n";
    OUTncl << "  res@tiXAxisFontHeightF = 0.020\n";
  }
  int nbCharY = eDrawArr.YAxisString.size();
  if (nbCharY > 0) {
    OUTncl << "  res@tiYAxisString   = \"" << eDrawArr.YAxisString << "\"\n";
    OUTncl << "  res@tiYAxisFontHeightF = 0.015\n";
  }
  std::vector<std::string> ListColors = {
      "black", "red",    "blue",       "purple", "green",  "cyan",
      "azure", "marron", "aquamarine", "violet", "DarkRed"};
  int nbColor = ListColors.size();
  if (nbArr > nbColor) {
    std::cerr << "Number of colors is insufficient\n";
    throw TerminalException{1};
  }
  std::string TheTotal;
  for (int iArr = 0; iArr < nbArr; iArr++) {
    if (iArr > 0)
      TheTotal += ",";
    TheTotal += "\"" + ListColors[iArr] + "\"";
  }
  OUTncl << "  res@xyLineColors = (/" << TheTotal << "/)\n";
  OUTncl << "  res@trYMaxF=" << eDrawArr.TheMax << "\n";
  OUTncl << "  res@trYMinF=" << eDrawArr.TheMin << "\n";
  OUTncl << "  res@trXMaxF=max(ListX)\n";
  OUTncl << "  res@trXMinF=min(ListX)\n";
  if (eDrawArr.DoExplicitLabel) {
    OUTncl << "  res@tmXBMode=\"Explicit\" \n";
    int nbLabel = eDrawArr.nbLabel;
    int nbEnt = eDrawArr.ListX.size();
    std::string strDay, strLabel;
    for (int iLabel = 0; iLabel < nbLabel; iLabel++) {
      double pos_d = static_cast<double>(iLabel) *
                     static_cast<double>(nbEnt - 1) /
                     static_cast<double>(nbLabel - 1);
      int pos_i;
      NearestInteger_double_int(pos_d, pos_i);
      pos_i = std::max(0, pos_i);
      pos_i = std::min(nbEnt - 1, pos_i);
      std::cerr << "iLabel=" << iLabel << " nbLabel=" << nbLabel
                << " nbEnt=" << nbEnt << " pos_d=" << pos_d
                << " pos_i=" << pos_i << "\n";
      if (iLabel > 0) {
        strDay += ", ";
        strLabel += ", ";
      }
      if (iLabel == 0)
        pos_i = 0;
      if (iLabel == nbLabel - 1)
        pos_i = nbEnt - 1;
      double eMJD = eDrawArr.ListX(pos_i);
      strDay += DoubleToString(eMJD);
      std::vector<int> eDate = DATE_ConvertMjd2six(eMJD);
      strLabel +=
          "\"" + DATE_ConvertSix2string_style(eDate, eDrawArr.StyleDate) + "\"";
    }
    OUTncl << "  res@tmXBValues=(/" << strDay << "/)\n";
    OUTncl << "  res@tmXBLabels=(/" << strLabel << "/)\n";
  }
  OUTncl << "  plot=gsn_csm_xy(wks,ListX,Data,res)\n";
  if (eDrawArr.DoTitle) {
    OUTncl << "  txresB             = True\n";
    OUTncl << "  txresB@txFontHeightF = 0.02\n";
    OUTncl << "  txresB@txFontColor = \"black\"\n";
    OUTncl << "  strLeft=\"\"\n";
    OUTncl << "  strMid=\"" << eDrawArr.TitleStr << "\"\n";
    OUTncl << "  strRight=\"\"\n";
    OUTncl << "  subtitles(wks, plot, strLeft, strMid, strRight, txresB)\n";
  }
  if (eDrawArr.ListName_plot.size() > 0) {
    OUTncl << "  lgres                    = True\n";
    OUTncl << "  lgres@lgLineColors     = (/" << TheTotal << "/)\n";
    OUTncl << "  lgres@lgItemType         = \"Lines\"\n";
    OUTncl << "  lgres@lgDashIndexes = (/";
    for (int iArr = 0; iArr < nbArr; iArr++) {
      if (iArr > 0)
        OUTncl << ",";
      OUTncl << "0";
    }
    OUTncl << "/)\n";
    OUTncl << "  lgres@lgLabelFontHeightF = .07\n";
    OUTncl << "  lgres@vpWidthF           = 0.11\n";
    OUTncl << "  lgres@vpHeightF          = 0.12\n";
    OUTncl << "  lgres@lgPerimOn = False\n";
    OUTncl << "  ;   lgres@lgPerimColor       = \"orange\"\n";
    OUTncl << "  lgres@lgPerimThicknessF  = 5.0\n";
    OUTncl << "  ListLabels= (/";
    for (int iArr = 0; iArr < nbArr; iArr++) {
      if (iArr > 0)
        OUTncl << ",";
      OUTncl << "\"" << eDrawArr.ListName_plot[iArr] << "\"";
    }
    OUTncl << "/)\n";
    OUTncl << "  lbid = gsn_create_legend(wks," << nbArr
           << ",ListLabels,lgres)\n";
    OUTncl << "  amres = True\n";
    OUTncl << "  amres@amParallelPosF   = 0.35\n";
    OUTncl << "  amres@amOrthogonalPosF = 0.30\n";
    OUTncl << "  annoid1 = gsn_add_annotation(plot,lbid,amres)\n";
  }
  ADD_ANNOTATION_TEXT(OUTncl, eDrawArr.TheAnnot);
  OUTncl << "  draw(plot)\n";
  OUTncl << "  frame(wks)\n";
  OUTncl << "end\n";
  OUTncl.close();
  //
  TripleNCL eTripl{eReal.eProgReal, TargetFile, eFileNC, eFileNCL,
                   ePerm.eBChoice};
  GeneralType eGen(eTripl);
  eCall.SubmitJob(eGen);
}

void LINES_PLOT_PYTHON(std::string const &FileName,
                       DrawLinesArr const &eDrawArr,
                       NCLcaller<GeneralType> &eCall,
                       PermanentInfoDrawing const &ePerm) {
  std::string TargetFile = FileName + "_storsave." + ePerm.Extension;
  bool InPlaceRun = ePerm.eBChoice.InPlaceRun;
  std::string eFileNC = FinalFile(InPlaceRun, TargetFile,
                                  ePerm.PrefixTemp.str() + "DataLines_" +
                                      eDrawArr.VarName + ".nc");
  std::string eFilePY = FinalFile(InPlaceRun, TargetFile,
                                  ePerm.PrefixTemp.str() + "ScriptLines_" +
                                      eDrawArr.VarName + ".py");
  LINES_DEFINE_NC(eFileNC, eDrawArr);
  int nbArr = eDrawArr.ListListVect.size();
  //
  if (nbArr != 1) {
    std::cerr << "The code works only for one single time series\n";
    throw TerminalException{1};
  }
  //
  if (!eDrawArr.IsTimeSeries) {
    std::cerr << "The code works only for time series\n";
    throw TerminalException{1};
  }
  //
  std::ofstream os(eFilePY);
  os << "  ;\n";
  os << "  ; Import statements\n";
  os << "  ;\n";
  os << "import matplotlib.pyplot as plt\n";
  os << "import datetime\n";
  os << "import numpy as np\n";
  os << "  ;\n";
  os << "  ; Data reading\n";
  os << "  ;\n";
  os << "InFile = \'" << FinalFileInScript(InPlaceRun, eFileNC) + "\'\n";
  os << "\n";
  os << "dataset = netCDF4.Dataset(InFile)\n";
  os << "ListTimeSix  = dataset.variables[\'ListTimeSix\'][:,:]\n";
  os << "ListListVect = dataset.variables[\'ListTimeVect\'][:,:]\n";
  os << "  ;\n";
  os << "  ; Construction of the arrays\n";
  os << "  ;\n";
  os << "LDim = ListTimeSix.shape\n";
  os << "nbEntry = LDim[0]\n";
  os << "ListTime = np.zeros(nbEntry)\n";
  os << "for i in range(nbEntry):\n";
  os << "    eYear  = ListTimeSix[i,0]\n";
  os << "    eMonth = ListTimeSix[i,1]\n";
  os << "    eDay   = ListTimeSix[i,2]\n";
  os << "    eHour  = ListTimeSix[i,3]\n";
  os << "    eMin   = ListTimeSix[i,4]\n";
  os << "    eSec   = ListTimeSix[i,5]\n";
  os << "    ListTime[i] = datetime.datetime(eYear, eMonth, eDay, eHour, eMin, "
        "eSec)\n";
  os << "LDimVect = ListListVect.shape\n";
  os << "nbVect = LDimVect[1]\n";
  os << "ListVal = np.zeros(nbEntry)\n";
  os << "for i in range(nbEntry):\n";
  os << "    ListVal = ListListVect[i,0]\n";
  int nbCharX = eDrawArr.XAxisString.size();
  if (nbCharX > 0) {
    std::cerr << "XAxisString need to be put\n";
  }
  int nbCharY = eDrawArr.YAxisString.size();
  if (nbCharY > 0) {
    std::cerr << "YAxisString need to be put\n";
  }
  if (eDrawArr.DoTitle) {
    std::cerr << "Need to write something here\n";
  }
  os << "  ;\n";
  os << "  ; Plotting the data\n";
  os << "  ;\n";
  os << "plt.title(\'" << eDrawArr.TitleStr << "\')\n";
  std::string OutFile = FileName + "_storsave." + ePerm.Extension;
  os << "plt.savefig(\'" << FinalFileInScript(InPlaceRun, OutFile) << "\')\n";
  os << "print \"After savefig\"\n";
  os.close();
  //
  TripleNCL eTripl{"python", TargetFile, eFileNC, eFilePY, ePerm.eBChoice};
  GeneralType eGen(eTripl);
  eCall.SubmitJob(eGen);
}

void LINES_PLOT(std::string const &FileName, DrawLinesArr const &eDrawArr,
                NCLcaller<GeneralType> &eCall,
                PermanentInfoDrawing const &ePerm) {
  std::string eProg = ePerm.ChoiceProgram.at("Lines_method");
  bool IsMatch = false;
  if (eProg == "ncl") {
    LINES_PLOT_NCL(FileName, eDrawArr, eCall, ePerm);
    IsMatch = true;
  }
  if (eProg == "python") {
    LINES_PLOT_PYTHON(FileName, eDrawArr, eCall, ePerm);
    IsMatch = true;
  }
  if (!IsMatch) {
    std::cerr << "eProg = " << eProg
              << " while allowed values are python and ncl\n";
    std::cerr << "If you want the crop, use pngcrop or pdfcrop in the type "
                 "instead of png, pdf\n";
    throw TerminalException{1};
  }
}

// clang-format off
#endif  // SRC_OCEAN_BASIC_PLOT_H_
// clang-format on
