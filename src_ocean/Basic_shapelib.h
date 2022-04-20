#ifndef SRC_OCEAN_BASIC_SHAPELIB_H_
#define SRC_OCEAN_BASIC_SHAPELIB_H_

#include "Basic_Ocean_types.h"
#include "Basic_file.h"
#include "Basic_netcdf.h"
#include "Temp_common.h"
#include "Timings.h"
#include <string>
#include <vector>
#include <shapefil.h>

struct SHPquad {
  double x;
  double y;
  double z;
  double m;
};

struct SHPpart {
  std::vector<SHPquad> eVectQuad;
  int i;
  int iPart;
  int ePartType;
  int nSHPType;
  int nShapeId;
};

struct ShapefileData {
  std::string eFile;
  std::vector<double> minbound;
  std::vector<double> maxbound;
  int nshp;
  int tshp;
  std::vector<SHPpart> ListBlock;
};

ShapefileData ExtractShapefileFromUnstructuredMesh(
    GridArray const &GrdArr, std::vector<std::vector<int>> const &eListList) {
  double minLon = GrdArr.GrdArrRho.LON.minCoeff();
  double maxLon = GrdArr.GrdArrRho.LON.maxCoeff();
  double minLat = GrdArr.GrdArrRho.LAT.minCoeff();
  double maxLat = GrdArr.GrdArrRho.LAT.maxCoeff();
  std::vector<double> minbound{minLon, minLat, 0, 0};
  std::vector<double> maxbound{maxLon, maxLat, 0, 0};
  int tshp = 3;
  int nbBlock = eListList.size();
  int nshp = nbBlock;
  //
  std::vector<SHPpart> ListBlock;
  int iIns = 0;
  int iInsPart = 0;
  for (auto &eList : eListList) {
    std::vector<SHPquad> eVectQuad;
    for (auto &eVert : eList) {
      double eX = GrdArr.GrdArrRho.LON(eVert, 0);
      double eY = GrdArr.GrdArrRho.LAT(eVert, 0);
      eVectQuad.push_back({eX, eY, 0, 0});
    }
    int ePartType = 5;
    int nSHPType = 3;
    int nShapeId = iIns;
    ListBlock.push_back(
        {eVectQuad, iIns, iInsPart, ePartType, nSHPType, nShapeId});
    iIns++;
  }
  //
  ShapefileData eSHP;
  eSHP.eFile = "unset";
  eSHP.minbound = minbound;
  eSHP.maxbound = maxbound;
  eSHP.nshp = nshp;
  eSHP.tshp = tshp;
  eSHP.ListBlock = ListBlock;
  return eSHP;
}

void PrintShapefileData(std::ostream &os, ShapefileData const &eShape) {
  os << "Before printing PrintShapefileData\n";
  os << "eFile=" << eShape.eFile << "\n";
  os << " After printing PrintShapefileData\n";
  for (int i = 0; i < 4; i++) {
    os << "i=" << i << " minBound=" << eShape.minbound[i]
       << " maxbound=" << eShape.maxbound[i] << "\n";
  }
  os << "nshp=" << eShape.nshp << " tshp=" << eShape.tshp << "\n";
  int nbBlock = eShape.ListBlock.size();
  os << "nbBlock=" << nbBlock << "\n";
  for (int iBlock = 0; iBlock < nbBlock; iBlock++) {
    os << "iBlock=" << iBlock << " / " << nbBlock << "\n";
    SHPpart ePart = eShape.ListBlock[iBlock];
    int len = ePart.eVectQuad.size();
    os << "  ePart.iPart=" << ePart.iPart << " i=" << ePart.i << " len=" << len
       << "\n";
    int ePartType = ePart.ePartType;
    int nSHPType = ePart.nSHPType;
    int nShapeId = ePart.nShapeId;
    os << "  ePartType=" << ePartType << " nSHPType=" << nSHPType
       << " nShapeId=" << nShapeId << "\n";
    std::vector<double> ListX(len), ListY(len), ListZ(len), ListM(len);
    for (int i = 0; i < len; i++) {
      ListX[i] = ePart.eVectQuad[i].x;
      ListY[i] = ePart.eVectQuad[i].y;
      ListZ[i] = ePart.eVectQuad[i].z;
      ListM[i] = ePart.eVectQuad[i].m;
    }
    os << "  X(min/max) = " << VectorMin(ListX) << " / " << VectorMax(ListX)
       << "\n";
    os << "  Y(min/max) = " << VectorMin(ListY) << " / " << VectorMax(ListY)
       << "\n";
    os << "  Z(min/max) = " << VectorMin(ListZ) << " / " << VectorMax(ListZ)
       << "\n";
    os << "  M(min/max) = " << VectorMin(ListM) << " / " << VectorMax(ListM)
       << "\n";
    double ErrCycle =
        fabs(ListX[0] - ListX[len - 1]) + fabs(ListY[0] - ListY[len - 1]);
    os << "  ErrCycle=" << ErrCycle << "\n";
  }
}

bool SHP_FileIsNull(SHPHandle const &shphandle) {
  if (shphandle == NULL)
    return true;
  return false;
}

int SHP_GetShapeType(std::string const &eFile) {
  /*  if (!IsExistingFile(eFile)) {
    std::cerr << "The file eFile = " << eFile << " is missing\n";
    throw TerminalException{1};
    }*/
  std::cerr << "eFile=" << eFile << "\n";
  std::cerr << "Before SHPOpen\n";
  SHPHandle shphandle = SHPOpen(eFile.c_str(), "rb");
  std::cerr << "After SHPOPen\n";
  if (SHP_FileIsNull(shphandle)) {
    std::cerr << "Error opening " << eFile << " for reading\n";
    std::cerr << "shphandle is null\n";
    throw TerminalException{1};
  }
  int nshp = 0, tshp = 0;
  double minbound[4], maxbound[4];
  std::cerr << "Before SHPGetInfo\n";
  SHPGetInfo(shphandle, &nshp, &tshp, minbound, maxbound);
  std::cerr << "nshp=" << nshp << " tshp=" << tshp << "\n";
  for (int i = 0; i < 4; i++) {
    std::cerr << "i=" << i << " minbound=" << minbound[i]
              << " maxbound=" << maxbound[i] << "\n";
  }
  std::cerr << "Before SHPClose\n";
  SHPClose(shphandle);
  std::cerr << " After SHPClose\n";
  return tshp;
}

ShapefileData SHP_ReadFromNetcdf(std::string const &eFile) {
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
  netCDF::NcVar dataPanPartStart = dataFile.getVar("panPartStart");
  netCDF::NcVar dataPanPartType = dataFile.getVar("panPartType");
  netCDF::NcVar dataList_nSHPType = dataFile.getVar("List_nSHPType");
  netCDF::NcVar dataList_nShapeId = dataFile.getVar("List_nShapeId");
  netCDF::NcDim eDim = dataPanPartStart.getDim(0);
  int nParts = eDim.getSize();
  std::vector<int> eVAR1(nParts);
  std::vector<int> ListPartStart(nParts);
  std::vector<int> ListPartType(nParts);
  std::vector<int> List_nSHPType(nParts);
  std::vector<int> List_nShapeId(nParts);
  dataPanPartStart.getVar(ListPartStart.data());
  dataPanPartType.getVar(ListPartType.data());
  dataList_nSHPType.getVar(List_nSHPType.data());
  dataList_nShapeId.getVar(List_nShapeId.data());
  //
  netCDF::NcVar padfX = dataFile.getVar("padfX");
  netCDF::NcVar padfY = dataFile.getVar("padfY");
  netCDF::NcVar padfZ = dataFile.getVar("padfZ");
  netCDF::NcVar padfM = dataFile.getVar("padfM");
  netCDF::NcDim fDim = padfX.getDim(0);
  int nVertices = fDim.getSize();
  std::vector<SHPquad> eVectQuadTot(nVertices);
  std::vector<double> eVARd(nVertices);
  padfX.getVar(eVARd.data());
  for (int i = 0; i < nVertices; i++)
    eVectQuadTot[i].x = eVARd[i];
  padfY.getVar(eVARd.data());
  for (int i = 0; i < nVertices; i++)
    eVectQuadTot[i].y = eVARd[i];
  padfZ.getVar(eVARd.data());
  for (int i = 0; i < nVertices; i++)
    eVectQuadTot[i].z = eVARd[i];
  padfM.getVar(eVARd.data());
  for (int i = 0; i < nVertices; i++)
    eVectQuadTot[i].m = eVARd[i];
  //
  std::vector<SHPpart> ListBlock;
  for (int iPart = 0; iPart < nParts; iPart++) {
    int eStart = ListPartStart[iPart];
    int eEnd;
    if (iPart < nParts - 1) {
      eEnd = ListPartStart[iPart + 1];
    } else {
      eEnd = nVertices;
    }
    std::vector<SHPquad> eVectQuad;
    for (int i = eStart; i < eEnd; i++)
      eVectQuad.push_back(eVectQuadTot[i]);
    int iPartIns = 0;
    int iIns = iPart;
    int ePartType = ListPartType[iPart];
    int nSHPType = List_nSHPType[iPart];
    int nShapeId = List_nShapeId[iPart];
    ListBlock.push_back(
        {eVectQuad, iIns, iPartIns, ePartType, nSHPType, nShapeId});
  }
  std::vector<double> ListX(nVertices), ListY(nVertices), ListZ(nVertices),
      ListM(nVertices);
  for (int i = 0; i < nVertices; i++) {
    ListX[i] = eVectQuadTot[i].x;
    ListY[i] = eVectQuadTot[i].y;
    ListZ[i] = eVectQuadTot[i].z;
    ListM[i] = eVectQuadTot[i].m;
  }
  std::vector<double> vectMinBound{VectorMin(ListX), VectorMin(ListY),
                                   VectorMin(ListZ), VectorMin(ListM)};
  std::vector<double> vectMaxBound{VectorMax(ListX), VectorMax(ListY),
                                   VectorMax(ListZ), VectorMax(ListM)};
  int tshp = 3;
  int nshp = nParts;
  ShapefileData eSHP;
  eSHP.eFile = eFile;
  eSHP.minbound = vectMinBound;
  eSHP.maxbound = vectMaxBound;
  eSHP.nshp = nshp;
  eSHP.tshp = tshp;
  eSHP.ListBlock = ListBlock;
  return eSHP;
}

void SHP_WriteGmshInputFile(std::string const &eFile,
                            ShapefileData const &eSHP) {
  std::ofstream os(eFile);
  os << std::setprecision(12);
  int pid = 0;
  int lid = 0;
  int nbPath = eSHP.ListBlock.size();
  std::vector<int> ListIndexLoop(nbPath, 0);
  for (int iPath = 0; iPath < nbPath; iPath++) {
    std::vector<SHPquad> eVectQuad = eSHP.ListBlock[iPath].eVectQuad;
    int lenPath = eVectQuad.size();
    int pidSave = pid;
    for (int i = 0; i < lenPath; i++) {
      double eX = eVectQuad[i].x;
      double eY = eVectQuad[i].y;
      pid++;
      os << "Point(" << pid << ") = {" << eX << ", " << eY << ", 0};\n";
    }
    int lidSave = lid;
    for (int i = 0; i < lenPath; i++) {
      int iNext;
      if (i == lenPath - 1) {
        iNext = 0;
      } else {
        iNext = i + 1;
      }
      int idx1 = pidSave + 1 + i;
      int idx2 = pidSave + 1 + iNext;
      lid++;
      os << "Line(" << lid << ") = {" << idx1 << ", " << idx2 << "};\n";
    }
    lid++;
    os << "Line Loop(" << lid << ") = {";
    ListIndexLoop[iPath] = lid;
    for (int i = 0; i < lenPath; i++) {
      if (i > 0)
        os << ",";
      int idx = lidSave + 1 + i;
      os << idx;
    }
    os << "};\n";
  }
  lid++;
  os << "Plane Surface(" << lid << ") = {";
  int lidSurf = lid;
  for (int iPath = 0; iPath < nbPath; iPath++) {
    if (iPath > 0)
      os << ",";
    os << ListIndexLoop[iPath];
  }
  os << "};\n";
  os << "Recombine Surface {" << lidSurf << "};\n";
}

void SHP_WriteNetcdfFile(std::string const &eFile, ShapefileData const &eSHP) {
  int nbBlock = eSHP.ListBlock.size();
  int nVertices = 0;
  std::vector<double> panPartStart_vect;
  std::vector<double> panPartType_vect;
  std::vector<double> List_nSHPType_vect;
  std::vector<double> List_nShapeId_vect;
  std::vector<SHPquad> eVectQuadTotal;
  for (int iBlock = 0; iBlock < nbBlock; iBlock++) {
    SHPpart eBlock = eSHP.ListBlock[iBlock];
    int len = eBlock.eVectQuad.size();
    panPartStart_vect.push_back(nVertices);
    panPartType_vect.push_back(eBlock.ePartType);
    List_nSHPType_vect.push_back(eBlock.nSHPType);
    List_nShapeId_vect.push_back(eBlock.nShapeId);
    nVertices += len;
    for (auto &eVal : eBlock.eVectQuad)
      eVectQuadTotal.push_back(eVal);
  }
  int nParts = panPartStart_vect.size();
  std::vector<double> padfX(nVertices), padfY(nVertices), padfZ(nVertices),
      padfM(nVertices);
  for (int i = 0; i < nVertices; i++) {
    padfX[i] = eVectQuadTotal[i].x;
    padfY[i] = eVectQuadTotal[i].y;
    padfZ[i] = eVectQuadTotal[i].z;
    padfM[i] = eVectQuadTotal[i].m;
  }

  if (!FILE_IsFileMakeable(eFile)) {
    std::cerr << "Request to create file eFile=" << eFile << "\n";
    std::cerr << "but the directory does not exist\n";
    throw TerminalException{1};
  }
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::replace, netCDF::NcFile::nc4);
  netCDF::NcDim dimPart = dataFile.addDim("nParts", nParts);
  netCDF::NcDim dimVert = dataFile.addDim("nVertices", nVertices);

  std::vector<std::string> ListDim1{"nParts"};
  std::vector<std::string> ListDim2{"nVertices"};

  netCDF::NcVar dataPanPartStart =
      dataFile.addVar("panPartStart", "int", ListDim1);
  netCDF::NcVar dataPanPartType =
      dataFile.addVar("panPartType", "int", ListDim1);
  netCDF::NcVar dataList_nSHPType =
      dataFile.addVar("List_nSHPType", "int", ListDim1);
  netCDF::NcVar dataList_nShapeId =
      dataFile.addVar("List_nShapeId", "int", ListDim1);
  //
  netCDF::NcVar padfX_vect = dataFile.addVar("padfX", "double", ListDim2);
  netCDF::NcVar padfY_vect = dataFile.addVar("padfY", "double", ListDim2);
  netCDF::NcVar padfZ_vect = dataFile.addVar("padfZ", "double", ListDim2);
  netCDF::NcVar padfM_vect = dataFile.addVar("padfM", "double", ListDim2);
  //
  dataPanPartStart.putVar(panPartStart_vect.data());
  dataPanPartType.putVar(panPartType_vect.data());
  dataList_nSHPType.putVar(List_nSHPType_vect.data());
  dataList_nShapeId.putVar(List_nShapeId_vect.data());
  //
  padfX_vect.putVar(padfX.data());
  padfY_vect.putVar(padfY.data());
  padfZ_vect.putVar(padfZ.data());
  padfM_vect.putVar(padfM.data());
}

ShapefileData SHP_ReadShapefile(std::string const &eFile) {
  if (!IsExistingFile(eFile)) {
    std::cerr << "The file eFile = " << eFile << " is missing\n";
    throw TerminalException{1};
  }
  SHPHandle shphandle = SHPOpen(eFile.c_str(), "rb");
  if (SHP_FileIsNull(shphandle)) {
    std::cerr << "Error opening " << eFile << " for reading\n";
    std::cerr << "shphandle is null\n";
    throw TerminalException{1};
  }
  int nshp, tshp;
  double minbound[4], maxbound[4];
  SHPGetInfo(shphandle, &nshp, &tshp, minbound, maxbound);
  if (nshp < 1) {
    std::cerr << "Error in shpgetinfo, wrong number of shapes nshp=" << nshp
              << "\n";
    throw TerminalException{1};
  }
  if (tshp == SHPT_NULL || tshp == SHPT_POINT || tshp == SHPT_MULTIPOINT ||
      tshp == SHPT_POINTZ || tshp == SHPT_MULTIPOINTZ || tshp == SHPT_POINTM ||
      tshp == SHPT_MULTIPOINTM || tshp == SHPT_MULTIPATCH) {
    std::cerr << "Wrong kind of patch\n";
    throw TerminalException{1};
  }
  std::vector<SHPpart> ListBlock;
  for (int i = 0; i < nshp; i++) {
    SHPObject *shpobj = SHPReadObject(shphandle, i);
    int nParts = shpobj->nParts;
    int nSHPType = shpobj->nSHPType;
    int nShapeId = shpobj->nShapeId;
    for (int iPart = 0; iPart < nParts; iPart++) {
      int startIndex = shpobj->panPartStart[iPart];
      int ePartType = shpobj->panPartType[iPart];
      int endIndex;
      if (iPart == nParts - 1) {
        endIndex = shpobj->nVertices;
      } else {
        endIndex = shpobj->panPartStart[iPart + 1];
      }
      std::vector<SHPquad> eVectQuad;
      for (int ielt = startIndex; ielt < endIndex; ielt++) {
        double eX = shpobj->padfX[ielt];
        double eY = shpobj->padfY[ielt];
        double eZ = shpobj->padfZ[ielt];
        double eM = shpobj->padfZ[ielt];
        eVectQuad.push_back({eX, eY, eZ, eM});
      }
      ListBlock.push_back({eVectQuad, i, iPart, ePartType, nSHPType, nShapeId});
    }
  }
  std::vector<double> vectMinBound(4), vectMaxBound(4);
  for (int i = 0; i < 4; i++) {
    vectMinBound[i] = minbound[i];
    vectMaxBound[i] = maxbound[i];
  }
  SHPClose(shphandle);
  ShapefileData eSHP;
  eSHP.eFile = eFile;
  eSHP.minbound = vectMinBound;
  eSHP.maxbound = vectMaxBound;
  eSHP.nshp = nshp;
  eSHP.tshp = tshp;
  eSHP.ListBlock = ListBlock;
  return eSHP;
}

void WriteShapefile(std::string const &eFile, ShapefileData const &eSHP) {
  SHPHandle shphandle = SHPCreate(eFile.c_str(), eSHP.tshp);
  if (SHP_FileIsNull(shphandle)) {
    std::cerr << "Error opening " << eFile << " for reading\n";
    std::cerr << "shphandle is null\n";
    throw TerminalException{1};
  }
  int nbBlock = eSHP.ListBlock.size();
  for (int ish = 0; ish < eSHP.nshp; ish++) {
    int nVertices = 0;
    std::vector<double> panPartStart_vect;
    std::vector<double> panPartType_vect;
    std::vector<double> List_nSHPType;
    std::vector<double> List_nShapeId;
    std::vector<SHPquad> eVectQuadTotal;
    for (int iBlock = 0; iBlock < nbBlock; iBlock++) {
      if (eSHP.ListBlock[iBlock].i == ish) {
        SHPpart eBlock = eSHP.ListBlock[iBlock];
        int len = eBlock.eVectQuad.size();
        panPartStart_vect.push_back(nVertices);
        panPartType_vect.push_back(eBlock.ePartType);
        List_nSHPType.push_back(eBlock.nSHPType);
        List_nShapeId.push_back(eBlock.nShapeId);
        nVertices += len;
        for (auto &eVal : eBlock.eVectQuad)
          eVectQuadTotal.push_back(eVal);
      }
    }
    int nParts = panPartStart_vect.size();
    std::vector<double> padfX(nVertices), padfY(nVertices), padfZ(nVertices),
        padfM(nVertices);
    for (int i = 0; i < nVertices; i++) {
      padfX[i] = eVectQuadTotal[i].x;
      padfY[i] = eVectQuadTotal[i].y;
      padfZ[i] = eVectQuadTotal[i].z;
      padfM[i] = eVectQuadTotal[i].m;
    }
    if (!IsVectorConstant(List_nSHPType)) {
      std::cerr << "List_nSHPType is not constant\n";
      throw TerminalException{1};
    }
    if (!IsVectorConstant(List_nShapeId)) {
      std::cerr << "List_nShapeId is not constant\n";
      throw TerminalException{1};
    }
    std::vector<int> panPartStart_i(nParts);
    std::vector<int> panPartType_i(nParts);
    for (int iPart = 0; iPart < nParts; iPart++) {
      panPartStart_i[iPart] = panPartStart_vect[iPart];
      panPartType_i[iPart] = panPartType_vect[iPart];
    }
    int nSHPType = List_nSHPType[0];
    int nShapeId = List_nShapeId[0];
    std::cerr << "ish=" << ish << "\n";
    std::cerr << "nShapeId=" << nShapeId << "  nSHPType=" << nSHPType << "\n";
    SHPObject *shpobj = SHPCreateObject(
        nSHPType, nShapeId, nParts, panPartStart_i.data(), panPartType_i.data(),
        nVertices, padfX.data(), padfY.data(), padfZ.data(), padfM.data());
    std::cerr << "After SHPCreateObject\n";
    int nShapeId_append = -1;
    SHPWriteObject(shphandle, nShapeId_append, shpobj);
    std::cerr << "After SHPWriteObject\n";
  }
  SHPClose(shphandle);
}

void DBF_PrintDatabaseInfo(std::ostream &os, std::string const &eFile) {
  DBFHandle dbfhandle = DBFOpen(eFile.c_str(), "rb");
  int nField = DBFGetFieldCount(dbfhandle);
  int nRecord = DBFGetRecordCount(dbfhandle);
  os << "nField=" << nField << " nRecord=" << nRecord << "\n";
  std::vector<std::string> ListType, ListName;
  for (int iField = 0; iField < nField; iField++) {
    char FieldName[30];
    int pnWidth, pnDecimals;
    DBFFieldType eType =
        DBFGetFieldInfo(dbfhandle, iField, FieldName, &pnWidth, &pnDecimals);
    std::string type;
    if (eType == FTString)
      type = "string";
    if (eType == FTInteger)
      type = "integer";
    if (eType == FTDouble)
      type = "double";
    if (eType == FTLogical)
      type = "logical";
    if (eType == FTInvalid)
      type = "invalid";
    std::string strFieldName = FieldName;
    os << "iField=" << iField << " name=" << strFieldName << " type=" << type
       << " pnWidth=" << pnWidth << " pnDecimals=" << pnDecimals << "\n";
    ListType.push_back(type);
    ListName.push_back(strFieldName);
  }
  for (int iRecord = 0; iRecord < nRecord; iRecord++) {
    os << "iRecord=" << iRecord << "/" << nRecord << "\n";
    for (int iField = 0; iField < nField; iField++) {
      std::string type = ListType[iField];
      std::string str;
      if (type == "string") {
        const char *ptr = DBFReadStringAttribute(dbfhandle, iRecord, iField);
        str = ptr;
      }
      if (type == "integer") {
        int val = DBFReadIntegerAttribute(dbfhandle, iRecord, iField);
        str = std::to_string(val);
      }
      if (type == "double") {
        double val = DBFReadDoubleAttribute(dbfhandle, iRecord, iField);
        str = std::to_string(val);
      }
      if (type == "logical") {
        const char *ptr = DBFReadLogicalAttribute(dbfhandle, iRecord, iField);
        str = ptr;
      }
      os << "  iField=" << iField << " type=" << type
         << " name=" << ListName[iField] << " val=" << str << "\n";
    }
  }
  DBFClose(dbfhandle);
}

struct DBFdataset {
  std::vector<std::string> ListIntVar;
  std::vector<std::string> ListDoubleVar;
  std::vector<std::string> ListBoolVar;
  std::vector<std::string> ListStringVar;
  //
  std::vector<int> ListIntValue;
  std::vector<double> ListDoubleValue;
  std::vector<bool> ListBoolValue;
  std::vector<std::string> ListStringValue;
};

DBFdataset DBF_ReadDataset(std::string const &eFile) {
  DBFHandle dbfhandle = DBFOpen(eFile.c_str(), "rb");
  int nField = DBFGetFieldCount(dbfhandle);
  int nRecord = DBFGetRecordCount(dbfhandle);
  std::cerr << "nField=" << nField << " nRecord=" << nRecord << "\n";
  std::vector<std::string> ListType;
  std::vector<std::string> ListIntVar, ListDoubleVar, ListBoolVar,
      ListStringVar;

  for (int iField = 0; iField < nField; iField++) {
    char FieldName[30];
    int pnWidth, pnDecimals;
    DBFFieldType eType =
        DBFGetFieldInfo(dbfhandle, iField, FieldName, &pnWidth, &pnDecimals);
    std::string strFieldName = FieldName;
    std::string type;
    if (eType == FTString) {
      type = "string";
      ListStringVar.push_back(strFieldName);
    }
    if (eType == FTInteger) {
      type = "integer";
      ListIntVar.push_back(strFieldName);
    }
    if (eType == FTDouble) {
      type = "double";
      ListDoubleVar.push_back(strFieldName);
    }
    if (eType == FTLogical) {
      type = "logical";
      ListDoubleVar.push_back(strFieldName);
    }
    if (eType == FTInvalid)
      type = "invalid";
    std::cerr << "iField=" << iField << " name=" << strFieldName
              << " type=" << type << " pnWidth=" << pnWidth
              << " pnDecimals=" << pnDecimals << "\n";
    ListType.push_back(type);
  }
  std::vector<int> ListIntValue;
  std::vector<double> ListDoubleValue;
  std::vector<bool> ListBoolValue;
  std::vector<std::string> ListStringValue;
  for (int iRecord = 0; iRecord < nRecord; iRecord++) {
    std::cerr << "iRecord=" << iRecord << "/" << nRecord << "\n";
    for (int iField = 0; iField < nField; iField++) {
      std::string type = ListType[iField];
      std::string str;
      if (type == "string") {
        const char *ptr = DBFReadStringAttribute(dbfhandle, iRecord, iField);
        ListStringValue.push_back(std::string(ptr));
      }
      if (type == "integer") {
        int val = DBFReadIntegerAttribute(dbfhandle, iRecord, iField);
        ListIntValue.push_back(val);
      }
      if (type == "double") {
        double val = DBFReadDoubleAttribute(dbfhandle, iRecord, iField);
        ListDoubleValue.push_back(val);
      }
      if (type == "logical") {
        const char *ptr = DBFReadLogicalAttribute(dbfhandle, iRecord, iField);
        std::cerr << "We need to write the corresponding code\n";
        std::cerr
            << "In dearth of an example, that is unfortunately quite hard\n";
        std::cerr << "But now you can because there is an exception\n";
        std::cerr << "str=" << std::string(ptr) << "\n";
        throw TerminalException{1};
      }
    }
  }
  DBFClose(dbfhandle);
  DBFdataset dds;
  dds.ListIntVar = ListIntVar;
  dds.ListDoubleVar = ListDoubleVar;
  dds.ListBoolVar = ListBoolVar;
  dds.ListStringVar = ListStringVar;
  //
  dds.ListIntValue = ListIntValue;
  dds.ListDoubleValue = ListDoubleValue;
  dds.ListBoolValue = ListBoolValue;
  dds.ListStringValue = ListStringValue;
  return dds;
}

void DBF_WriteDataset(std::string const &eFile, DBFdataset const &dds) {
  DBFHandle dbfhandle = DBFCreate(eFile.c_str());
  for (auto &eVar : dds.ListIntVar)
    DBFAddField(dbfhandle, eVar.c_str(), FTInteger, 10, 0);
  for (auto &eVar : dds.ListDoubleVar)
    DBFAddField(dbfhandle, eVar.c_str(), FTDouble, 11, 0);
  for (auto &eVar : dds.ListStringVar)
    DBFAddField(dbfhandle, eVar.c_str(), FTString, 80, 0);
  for (auto &eVar : dds.ListBoolVar)
    DBFAddField(dbfhandle, eVar.c_str(), FTLogical, 11,
                0); // We do not really know actually
  //
  std::vector<int> ListPossibleLength;
  int siz1, siz2, siz3;
  //
  siz1 = dds.ListIntVar.size();
  if (siz1 > 0) {
    siz2 = dds.ListIntValue.size();
    siz3 = siz2 / siz1;
    ListPossibleLength.push_back(siz3);
  }
  //
  siz1 = dds.ListDoubleVar.size();
  if (siz1 > 0) {
    siz2 = dds.ListDoubleValue.size();
    siz3 = siz2 / siz1;
    ListPossibleLength.push_back(siz3);
  }
  //
  siz1 = dds.ListBoolVar.size();
  if (siz1 > 0) {
    siz2 = dds.ListBoolValue.size();
    siz3 = siz2 / siz1;
    ListPossibleLength.push_back(siz3);
  }
  //
  siz1 = dds.ListStringVar.size();
  if (siz1 > 0) {
    siz2 = dds.ListStringValue.size();
    siz3 = siz2 / siz1;
    ListPossibleLength.push_back(siz3);
  }
  //
  if (ListPossibleLength.size() == 0) {
    std::cerr << "The set ListPossibleLength should be nontrivial\n";
    throw TerminalException{1};
  }
  if (VectorMin(ListPossibleLength) != VectorMax(ListPossibleLength)) {
    std::cerr << "The length are incoherents\n";
    throw TerminalException{1};
  }
  int nRecord = ListPossibleLength[0];
  //
  int idxInt = 0;
  int idxDouble = 0;
  int idxBool = 0;
  int idxString = 0;
  int nbIntVar = dds.ListIntVar.size();
  int nbDoubleVar = dds.ListDoubleVar.size();
  int nbStringVar = dds.ListStringVar.size();
  int nbBoolVar = dds.ListBoolVar.size();
  for (int iRecord = 0; iRecord < nRecord; iRecord++) {
    int pos = 0;
    for (int i = 0; i < nbIntVar; i++) {
      int ins = DBFWriteIntegerAttribute(dbfhandle, iRecord, pos,
                                         dds.ListIntValue[idxInt]);
      if (ins == 0) {
        std::cerr << "Error with DBFWriteIntegerAttribute\n";
        throw TerminalException{1};
      }
      idxInt++;
      pos++;
    }
    for (int i = 0; i < nbDoubleVar; i++) {
      int ins = DBFWriteDoubleAttribute(dbfhandle, iRecord, pos,
                                        dds.ListDoubleValue[idxDouble]);
      if (ins == 0) {
        std::cerr << "Error with DBFWriteDoubleAttribute\n";
        throw TerminalException{1};
      }
      idxDouble++;
      pos++;
    }
    for (int i = 0; i < nbStringVar; i++) {
      int ins = DBFWriteStringAttribute(dbfhandle, iRecord, pos,
                                        dds.ListStringValue[idxDouble].c_str());
      if (ins == 0) {
        std::cerr << "Error with DBFWriteStringAttribute\n";
        throw TerminalException{1};
      }
      idxString++;
      pos++;
    }
    for (int i = 0; i < nbBoolVar; i++) {
      idxBool++;
      pos++;
      std::cerr << "Need to write the code\n";
      throw TerminalException{1};
    }
  }
  DBFClose(dbfhandle);
}

ShapefileData ReadCoastlineInformation(std::string const &eFile) {
  std::string eExtension = FILE_GetExtension(eFile);
  if (eExtension == "nc") {
    return SHP_ReadFromNetcdf(eFile);
  }
  if (eExtension == "shp") {
    return SHP_ReadShapefile(eFile);
  }
  std::cerr << "Error in ReadCoastlineInformation\n";
  std::cerr << "No matching extension\n";
  std::cerr << "eExtension=" << eExtension << "\n";
  throw TerminalException{1};
}

void WriteCoastlineInformation(std::string const &eFile,
                               ShapefileData const &eSHP) {
  std::string eExtension = FILE_GetExtension(eFile);
  if (eExtension == "shp") {
    std::string eFileRed = FILE_RemoveEndingExtension(eFile, "shp");
    std::string eFileRed_shp = eFileRed + ".shp";
    std::string eFileRed_dbf = eFileRed + ".dbf";
    WriteShapefile(eFileRed_shp, eSHP);
    int nbBlock = eSHP.ListBlock.size();
    DBFdataset dds;
    dds.ListIntVar = {"physical", "entity"};
    std::vector<int> ListIntValue(2 * nbBlock, 1);
    dds.ListIntValue = ListIntValue;
    DBF_WriteDataset(eFileRed_dbf, dds);
    return;
  }
  if (eExtension == "nc") {
    SHP_WriteNetcdfFile(eFile, eSHP);
    return;
  }
  if (eExtension == "geo") {
    SHP_WriteGmshInputFile(eFile, eSHP);
    return;
  }
  std::cerr << "Error in WriteCoastlineInformation\n";
  std::cerr << "No matching extension\n";
  std::cerr << "eExtension=" << eExtension << "\n";
  throw TerminalException{1};
}

#endif  // SRC_OCEAN_BASIC_SHAPELIB_H_
