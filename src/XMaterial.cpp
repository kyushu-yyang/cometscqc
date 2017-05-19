#include <iostream>
#include "XQuenchExcept.hpp"
#include "XQuenchLogger.hpp"
#include "XMaterial.hpp"

using Quench::XQuenchLogger;

void XMaterial :: SetMaterialProperty(const double T, const double RRR, const double B)
{
  fTemp = T;
  fRRR  = RRR;
  fFld  = B;
}


void XMaterial :: Print()
{
  QuenchError( XQuenchLogger::INFO, "Temperature:" << fTemp << " , RRR: " << fRRR
                                    << " , Field: " << fFld );
  QuenchError( XQuenchLogger::INFO, "Density: " << this->GetDensity() );
  QuenchError( XQuenchLogger::INFO, "Electric Resistivity: " << this->GetResistivity() );
  QuenchError( XQuenchLogger::INFO, "Thermal Conductivity: " << this->GetConductivity() );
  QuenchError( XQuenchLogger::INFO, "Heat Capacity: " << this->GetCapacity() );
}


double XMaterial :: GetSerialk(const double l1, const double k1,
                               const double l2, const double k2) const
{
  const double k = (l1 + l2) / (l1/k1 + l2/k2);
  return k;
}


double XMaterial :: GetParallelk(const double A1, const double k1,
                                 const double A2, const double k2) const
{
  const double k = (A1*k1 + A2*k2) / (A1 + A2);
  return k;
}
