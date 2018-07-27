#include "River.h"
int main(int argc, char *argv[])
{
  try {
    if (argc <= 3) {
      std::cerr << "MERGE_RiverFiles is used as\n";
      std::cerr << "MERGE_RiverFiles [TargetFile] [file1] [file2] .... [fileN]\n";
      std::cerr << "with file.nml the file describing the chosen options\n";
      return -1;
    }
    std::string RiverFile=argv[1];
    std::vector<std::string> ListRiverFile;
    for (int i=2; i<argc; i++) {
      std::string eFile=argv[i];
      ListRiverFile.push_back(eFile);
    }
    MergeRiverFile(RiverFile, ListRiverFile);
    std::cerr << "Normal termination of MERGE_RiverFiles\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Error in MERGE_RiverFiles\n";
    exit(e.eVal);
  }
}
