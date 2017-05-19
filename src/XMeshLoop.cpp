#include "XQuenchExcept.hpp"
#include "XQuenchLogger.hpp"
#include "XMeshLoop.hpp"

using Quench::XQuenchLogger;

XMeshLoop :: XMeshLoop()
    : fMshZ(-1),
      fMshP(-1),
      fMshR(-1)
{}


XMeshLoop :: XMeshLoop(const int mz, const int mp, const int mr)
    : fMshZ(mz),
      fMshP(mp),
      fMshR(mr)
{}


XMeshLoop :: ~XMeshLoop() {}


void XMeshLoop :: SetMesh(const int mz, const int mp, const int mr)
{
  fMshZ = mz;
  fMshP = mp;
  fMshR = mr;
}


size_t XMeshLoop :: GetMesh(const Coil dim) const
{
  switch (dim) {
    case iZ:   return fMshZ;
    case iPhi: return fMshP;
    case iR:   return fMshR;
    default: return 0;
  }
}


bool XMeshLoop :: IsOverRange(const int id) const
{
  bool over_range = false;
  const int size = (fMshZ+2) * (fMshP+2) * (fMshR+2);
  if (id>=size)  over_range = true;

  return over_range;
}


const int XMeshLoop :: Id(const int i, const int j, const int k)
{
  const int id = k*(fMshZ+2)*(fMshP+2) + j*(fMshZ+2) + i;

  if ( IsOverRange(id) ) {
    QuenchError( XQuenchLogger::ERROR,"over range: " << 
                                      "id = " << id << ", i = " << i <<
                                      ", j = " << j  << ", k = " << k );
    XQuenchExcept except("id number is over range.");
    throw except;
  }
    
  return id;
}


void XMeshLoop :: StartLoop()
{
   int i = 0;  int j = 0;  int k = 0;

   for ( ; k<fMshR+2; k++) {
     // user defined function
     InLoopR(i, j, k);
     for ( ; j<fMshP+2; j++) {
       // user defined function
       InLoopPhi(i, j, k);
       for ( ; i<fMshZ+2; i++) {
         // user defined function
         InLoopZ(i, j, k);
       }
     }
   }
   // finished mesh loop
}
