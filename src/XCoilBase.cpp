#include "XCoilBase.hpp"

using Quench::XCoilBase;

std::string XCoilBase :: GetTypeName()
{
  std::string name = "";
  const Geometry mat = this->GetType();

  switch (mat) {
    case kConductor: name = "Conductor"; break;
    case     kStrip: name = "Strip";     break;
    case     kShell: name = "Shell";     break;
    case       kG10: name = "G10";       break;
    case     kA5083: name = "A5083";     break;
    default: break;
  }

  return name;
}

