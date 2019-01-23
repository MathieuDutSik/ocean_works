//#include "Model_interpolation.h"

#include "Namelist.h"
#include "Basic_netcdf.h"
#include "Basic_string.h"
#include "Basic_plot.h"
#include "NamelistExampleOcean.h"
#include "SphericalGeom.h"
#include "Model_grids.h"
//#include "Floats.h"

int main(int argc, char *argv[])
{
  try {
    FullNamelist eFull=NAMELIST_GetStandard_CREATE_TracerSourceTerm();
    if (argc != 2) {
      std::cerr << "CREATE_TracerSourceTerm is used as\n";
      std::cerr << "CREATE_TracerSourceTerm [file.nml]\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull);
      return -1;
    }
    std::string eFileName=argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    SingleBlock eBlPROC=eFull.ListBlock.at("PROC");
    std::string OutFile=eBlPROC.ListStringValues.at("OutFile");
    std::string eModelName=eBlPROC.ListStringValues.at("MODELNAME");
    CHECK_Model_Allowedness(eModelName);
    std::string GridFile=eBlPROC.ListStringValues.at("GridFile");
    TripleModelDesc eTriple{eModelName, GridFile, "unset", "unset", {}};
    GridArray GrdArr=RETRIEVE_GRID_ARRAY(eTriple);
    bool IsSpherical = GrdArr.IsSpherical;
    std::vector<double> ListLON=eBlPROC.ListListDoubleValues.at("ListLON");
    std::vector<double> ListLAT=eBlPROC.ListListDoubleValues.at("ListLAT");
    std::vector<double> ListDistM=eBlPROC.ListListDoubleValues.at("ListDistM");
    int nbTracer=ListLON.size();
    std::cerr << "nbTracer=" << nbTracer << "\n";
    //
    int nbNode=GrdArr.GrdArrRho.LON.size();
    std::cerr << "nbNode=" << nbNode << "\n";
    //
    std::vector<int> SOURCE_NODE;
    std::vector<int> ListNbMatch(nbTracer, 0);
    for (int iNode=0; iNode<nbNode; iNode++) {
      bool IsCorr=false;
      for (int iTracer=0; iTracer<nbTracer; iTracer++) {
	double eLon1=GrdArr.GrdArrRho.LON(iNode,0);
	double eLat1=GrdArr.GrdArrRho.LAT(iNode,0);
	double eLon2=ListLON[iTracer];
	double eLat2=ListLAT[iTracer];
	//      std::cerr << "------ iNode=" << iNode << " iTracer=" << iTracer << "\n";
	//      std::cerr << "eLon12=" << eLon1 << "," << eLon2 << " eLat12=" << eLat1 << "," << eLat2 << "\n";
	double eDistM=GeodesicDistanceM_General(eLon1, eLat1, eLon2, eLat2, IsSpherical);
	//      std::cerr << "eDistKM=" << eDistKM << " ListDist=" << ListDistKM[iTracer] << "\n";
	if (eDistM < ListDistM[iTracer]) {
	  IsCorr=true;
	  ListNbMatch[iTracer]++;
	}
      }
      if (IsCorr)
	SOURCE_NODE.push_back(iNode);
    }
    for (int iTracer=0; iTracer<nbTracer; iTracer++) {
    std::cerr << "iTracer=" << iTracer << " nbMatch=" << ListNbMatch[iTracer] << "\n";
    }
    int NBI=SOURCE_NODE.size();
    std::cerr << "NBI=" << NBI << "\n";
    int *Anode;
    Anode = new int[NBI];
    for (int iB=0; iB<NBI; iB++)
      Anode[iB]=SOURCE_NODE[iB]+1;
    //
    MyMatrix<double> eMat(nbTracer, NBI);
    for (int iTracer=0; iTracer<nbTracer; iTracer++) {
      for (int iB=0; iB<NBI; iB++) {
	int iNode=SOURCE_NODE[iB];
	double eLon1=GrdArr.GrdArrRho.LON(iNode,0);
	double eLat1=GrdArr.GrdArrRho.LAT(iNode,0);
	double eLon2=ListLON[iTracer];
	double eLat2=ListLAT[iTracer];
	double eDistM=GeodesicDistanceM_General(eLon1, eLat1, eLon2, eLat2, IsSpherical);
	double eVal;
	if (eDistM < ListDistM[iTracer]) {
	  eVal=1;
	}
	else {
	  eVal=0;
	}
	eMat(iTracer, iB)=eVal;
      }
    }
    //
    double *Aconc;
    Aconc = new double[nbTracer*NBI];
    int idx=0;
    for (int iB=0; iB<NBI; iB++) {
      for (int iTracer=0; iTracer<nbTracer; iTracer++) {
	Aconc[idx]=eMat(iTracer, iB);
	//      std::cerr << "iTracer=" << iTracer << "  iB=" << iB << " eCons=" << eMat(iTracer,iB) << "\n";
	idx++;
      }
    }
    //
    // Now the netcdf proper
    //
    netCDF::NcFile dataFile(OutFile, netCDF::NcFile::replace, netCDF::NcFile::nc4);
    netCDF::NcDim eDimNBI=dataFile.addDim("NBI", NBI);
    netCDF::NcDim eDimTracer=dataFile.addDim("NTR", nbTracer);
    //
    std::vector<std::string> ListDim1{"NBI"};
    std::vector<std::string> ListDim2{"NBI", "NTR"};
    //
    netCDF::NcVar eVar1=dataFile.addVar("SOURCE_NODE", "int", ListDim1);
    netCDF::NcVar eVar2=dataFile.addVar("SOURCE_CONC", "double", ListDim2);
    eVar1.putVar(Anode);
    eVar2.putVar(Aconc);
    delete [] Anode;
    delete [] Aconc;
    std::cerr << "Normal termination of CREATE_TracerSourceTerm\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Error in CREATE_TracerSourceTerm\n";
    exit(e.eVal);
  }
}
