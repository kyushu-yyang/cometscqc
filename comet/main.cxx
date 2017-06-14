#include <iostream>
#include <TString.h>

#include "XQuenchExcept.hpp"
#include "XQuenchLogger.hpp"
#include "XCOMETConstruction.h"

using namespace Quench;

int main(int argc, char** argv)
{
 
  XQuenchLogger* log = XQuenchLogger::GetInstance();
  log->Start(XQuenchLogger::DEBUG, "quenchlogger.log");
  //

  XCOMETConstruction* comet = new XCOMETConstruction();

  try {
    comet->SetTime(0.*sec, 80.*sec, 4.e-3*msec);
    //comet->SetDisplayStep(10);
    comet->SetOperationTime(60*day);
    comet->SetHotSpot(45,1,2,6e+3);
    comet->SetCurrent(2700.*Amp);
    comet->SetDumpResistor(0.185*Ohm);
    comet->SetInductance(12.69);
    //comet->SetDiode(0.7);
    comet->SetThreshold(0.1);
    comet->SetDetectTime(0.1*sec);

    comet->ConstructCS0( Form("%s/170512cs0x552000.txt",argv[1]), Form("%s/tempCS0.dat",argv[1]) );
    comet->ConstructCS1( Form("%s/170512cs1x552000.txt",argv[1]), Form("%s/tempCS1.dat",argv[1]) );
    comet->ConstructMS1( Form("%s/170512ms1x552000.txt",argv[1]), Form("%s/tempMS1.dat",argv[1]) );
    comet->ConstructMS2( Form("%s/170512ms2x552000.txt",argv[1]), Form("%s/tempMS2.dat",argv[1]) );
    comet->Begin();
    comet->Run();
    comet->End();
  }
  catch (XQuenchExcept except) {
    std::cerr << "Error: " << except.what() << std::endl;
    delete comet;
  }

  return 0;
}
