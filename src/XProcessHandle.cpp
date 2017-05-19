#include <cmath>
#include "XQuenchLogger.hpp"
#include "XQuenchExcept.hpp"
#include "XProcessHandle.hpp"

using Quench::XQuenchLogger;
using Quench::XProcessHandle;

XProcessHandle :: XProcessHandle()
    : Quench::XProcessManager(),
      fCurrent(2700.)
{}

XProcessHandle :: ~XProcessHandle()
{}

void XProcessHandle :: SetDecayParameter(const double I, const double B)
{
  fCurrent = I;
  fField   = B;
}
