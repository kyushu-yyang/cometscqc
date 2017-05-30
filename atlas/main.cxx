#include <iostream>
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
    comet->SetTime(0.*sec, 60.*sec, 4.e-5*msec);
    //comet->SetDisplayStep(10);
    comet->SetCurrent(7600.*Amp);
    comet->SetHotSpot(800/2, 1, 2);
    comet->SetDumpResistor(0.*Ohm);
    comet->SetInductance(1.35);
    comet->SetDiode(0.);
    comet->SetThreshold(0.05);
    comet->SetDetectTime(0.1*sec);

    comet->ConstructAtlas();
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
