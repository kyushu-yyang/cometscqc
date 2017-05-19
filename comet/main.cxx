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
    comet->SetOperationTime(90*day);
    comet->SetHotSpot(24,1,2);
    comet->SetCurrent(2700.*Amp);
    comet->SetDumpResistor(0.185*Ohm);
    comet->SetInductance(12.69);
    //comet->SetDiode(0.7);
    comet->SetThreshold(0.1);
    comet->SetDetectTime(0.1*sec);

    comet->ConstructCS0( Form("%s/phits288/161029CS0Track.dat",argv[1]), Form("%s/tempdis90/cs0.txt",argv[1]) );
    comet->ConstructCS1( Form("%s/phits288/161029CS1Track.dat",argv[1]), Form("%s/tempdis90/cs1.txt",argv[1]) );
    comet->ConstructMS1( Form("%s/phits288/161029MS1Track.dat",argv[1]), Form("%s/tempdis90/ms1.txt",argv[1]) );
    comet->ConstructMS2( Form("%s/phits288/161029MS2Track.dat",argv[1]), Form("%s/tempdis90/ms2.txt",argv[1]) );
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
