#ifndef SRC_OCEAN_EMAILBASHSYSTEM_H_
#define SRC_OCEAN_EMAILBASHSYSTEM_H_

#include <string>
#include <vector>

struct SendMailOper {
  std::string DestEmail;
  std::string FromEmail;
  std::string Subject;
  std::string Text;
  std::vector<std::string> ListFile;
};

void SendingEmail(SendMailOper const &eSendOper) {
  std::string TheFile = "/tmp/ScriptEmail_" + random_string(20) + ".perl";
  std::cerr << "TheFile = " << TheFile << "\n";
  {
    std::ofstream os(TheFile);
    os << "use MIME::Lite;\n";
    os << "my $msg = MIME::Lite->new(\n";
    os << "    From    => '" << eSendOper.FromEmail << "',\n";
    os << "    To      => '" << eSendOper.DestEmail << "',\n";
    os << "    Subject => '" << eSendOper.Subject << "',\n";
    os << "    Type    => 'multipart/mixed',\n";
    os << ");\n";
    os << "\n\n";
    os << "$msg->attach(\n";
    os << "    Type     => 'TEXT',\n";
    os << "    Data     => \"" << eSendOper.Text << "\",\n";
    os << ");\n";
    os << "\n\n";
    for (auto &eFile : eSendOper.ListFile) {
      std::string eExtension = FILE_GetExtension(eFile);
      std::string eType = "unset";
      if (eExtension == "png")
        eType = "image/png";
      std::string eFileName = FILE_GetNakedFilename(eFile);
      os << "$msg->attach(\n";
      os << "    Type     => '" << eType << "',\n";
      os << "    Path     => '" << eFile << "',\n";
      os << "    Filename => '" << eFileName << "',\n";
      os << ");\n";
    }
    os << "$msg->send;\n";
  }
  std::string TheCommand = "perl " + TheFile;
  std::cerr << "TheCommand = " << TheCommand << "\n";
  int iret = system(TheCommand.c_str());
  if (iret != 0) {
    std::cerr << "unable to run the TheCommand=" << TheCommand << "\n";
    throw TerminalException{1};
  }
  RemoveFileIfExist(TheFile);
  std::cerr << "Exiting the Email submission system\n";
}

void CALL_SEND_ATTACHMENT(SendMailOper const &eSendOper) {
  std::cerr << "Beginning of CALL_SEND_ATTACHMENT\n";
  while (true) {
    bool IsOK = true;
    for (auto &eFile : eSendOper.ListFile)
      if (!IsExistingFile(eFile))
        IsOK = false;
    std::cerr << "IsOK=" << IsOK << "\n";
    for (auto &eFile : eSendOper.ListFile)
      std::cerr << "  eFile=" << eFile << " test=" << IsExistingFile(eFile)
                << "\n";
    sleep(5);
    if (IsOK) {
      std::cerr << "Before call to SendingEmail\n";
      SendingEmail(eSendOper);
      break;
    }
  }
  std::cerr << "Exiting of CALL_SEND_ATTACHMENT\n";
}

struct BashOper {
  std::string TheCommand;
  std::vector<std::string> ListFile;
  int CorrectRetVal;
};

void CALL_BASH_OPERATION(BashOper const &eBashOper) {
  std::cerr << "Beginning of CALL_BASH_OPERATION\n";
  while (true) {
    bool IsOK = true;
    for (auto &eFile : eBashOper.ListFile) {
      bool test = IsExistingFile(eFile);
      std::cerr << "  eFile=" << eFile << " test=" << test << "\n";
      if (!test)
        IsOK = false;
    }
    std::cerr << "IsOK=" << IsOK << "\n";
    if (IsOK) {
      std::cerr << "TheCommand=" << eBashOper.TheCommand << "\n";
      int iret = system(eBashOper.TheCommand.c_str());
      if (iret != eBashOper.CorrectRetVal) {
        std::cerr << "Error with return value iret=" << iret << "\n";
        std::cerr << "eBashOper.CorrectRetVal=" << eBashOper.CorrectRetVal
                  << "\n";
        throw TerminalException{1};
      }
      break;
    }
    sleep(5);
  }
  std::cerr << "Exiting of CALL_BASH_OPERATION\n";
}

#endif // SRC_OCEAN_EMAILBASHSYSTEM_H_
