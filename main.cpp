#include <iostream>
#include <TApplication.h>
#include "XQuenchLogger.hpp"
#include "XQuenchExcept.hpp"
//#include "XQuenchContainer.hpp"
#include "XMatAluminium.hpp"
#include "XFieldHandle.hpp"
#include "XPostField.hpp"
#include "XCoilHandle.hpp"
#include "XCoilConductor.hpp"
#include "XPlotMaterial.hpp"
#include "IFdmUnits.hpp"

using namespace Quench;

void Test() {
  XQuenchLogger* log = XQuenchLogger::GetInstance();
  log->Start(XQuenchLogger::DEBUG, "quench.log");

  //XPostField* post = new XPostField();
  XPlotMaterial* mat = new XPlotMaterial();
  mat->SetOption("R");
  mat->SetRange(4.5, 300.);
  XMatAluminium* al = new XMatAluminium();
  al->SetMaterialProperty(4.3, 500, 5);

  mat->SetMaterial(al);
  mat->Add("B", 4.);
  mat->Add("B", 3.);
  mat->Add("B", 2.);
  mat->Add("B", 1.);
  mat->Plot();

  XFieldHandle* fld = new XFieldHandle();
  XCoilHandle* coil = new XCoilHandle();
  coil->SetMesh(270, 2, 9);

  try {
    fld->AddCoil("CS1", -79.525*cm, 59.525*cm, 672.*mm, 823.65*mm);
    fld->SetCurrent(2700.);

    //post->SetFieldHandler(fld);
    //post->Plot();
  }
  catch (XQuenchExcept except) {
    delete fld;
    delete coil;
    std::cerr << " ERROR: " << except.what() << std::endl;
  }

  log->Stop();
}

int main(int argc, char** argv)
{
  TApplication* app = new TApplication("app", &argc, argv);

  Test();

  app->Run();

  return 0;
}
