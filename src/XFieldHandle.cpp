#include <iostream>
#include "XQuenchLogger.hpp"
#include "XFieldHandle.hpp"

using Quench::XQuenchLogger;
using Quench::XFieldHandle;
using Quench::XMagnetInfoContainer;
using Quench::XFieldContainer;

XFieldHandle :: XFieldHandle()
    : fTarget(), fCurrent(2700.)
{}


XFieldHandle :: XFieldHandle(const std::string& name)
    : fTarget(name), fCurrent(2700.)
{
  QuenchError( XQuenchLogger::INFO, "calculate the field of magnet: " << name );
}


XFieldHandle :: ~XFieldHandle()
{}


bool XFieldHandle :: is_exist(const std::string& name) const
{
  bool exist = false;

  for (std::vector<int>::size_type i=0; i<fHC.size(); i++) {
    if ( fHC[i]->GetName() == name )  {
      exist = true;
      break;
    }
  }

  return exist;
}


void XFieldHandle :: AddCoil(const std::string& name, const double z0, const double z1,
                                                      const double r0, const double r1)
{
  if ( is_exist(name) )
    QuenchError( XQuenchLogger::WARNING, "magnet: " << name << " is already defined." );
        
  XMagnetInfoContainer* mag = new XMagnetInfoContainer();
  mag->SetName( name );
  mag->SetDimension( z0, z1, r0, r1 );
  fHC.push_back(mag);
  //delete mag;
}


void XFieldHandle :: AddCoil(XMagnetInfoContainer* mag)
{
  if ( !is_exist(mag->GetName()) )
    fHC.push_back(mag);
}


void XFieldHandle :: SetMesh(const std::string& name, const int mz, const int mr)
{
  if ( !is_exist(name) )
    QuenchError( XQuenchLogger::WARNING, "magnet: " << name << " did not exist." );

  for (std::vector<int>::size_type i=0; i<fHC.size(); i++) {
    if (fHC.at(i)->GetName() == name) {
      fHC.at(i)->SetMesh(mz, mr);
      break;
    }
  }
} 


XMagnetInfoContainer* XFieldHandle :: GetInfoEntry(const std::string& name)
{
  for (std::vector<int>::size_type i=0; i<fHC.size(); i++) {
    if ( fHC.at(i)->GetName()==name )
      return fHC.at(i);
  }
  
  QuenchError( XQuenchLogger::WARNING, "magnet: " << name << " did not exist." );

  return NULL;
}


void XFieldHandle :: calfield(const std::string& name)
{
  XMagneticField* fld = new XBiotSavart();  

  /// find this magnet
  for (std::vector<int>::size_type i=0; i<fHC.size(); i++) {
    if ( fHC.at(i)->GetName()==name ) {
      fld->SetMapRange( fHC.at(i)->GetDimension()[2], fHC.at(i)->GetDimension()[3],
                        fHC.at(i)->GetDimension()[0], fHC.at(i)->GetDimension()[1] );
      fld->SetMapMesh( fHC.at(i)->GetMesh()[1], fHC.at(i)->GetMesh()[0] );
      break;
    }
  }

  /// calculate self field and mutal field
  double Bx, By;

  for (std::vector<int>::size_type i=0; i<fHC.size(); i++) {
    QuenchError( XQuenchLogger::DEBUG, "calculating magnet: " << fHC.at(i)->GetName() );
    
    fld->SetSolenoid( fHC.at(i)->GetDimension()[2], fHC.at(i)->GetDimension()[3],
                      fHC.at(i)->GetDimension()[0], fHC.at(i)->GetDimension()[1] );
    fld->SetSolenoidMesh( fHC.at(i)->GetMesh()[1], 50, fHC.at(i)->GetMesh()[0] );
    fld->SetCurrent( fCurrent );
    fld->Run();

    if (i==0)
      /// self field
      fCollect = fld->GetFieldContainer();
    else {
      /// mutal field
      for (std::vector<int>::size_type j=0; j<fCollect.size(); j++) {
        Bx = fCollect.at(j)->GetField().at(0) + fld->GetFieldContainer().at(j)->GetField().at(0);
        By = fCollect.at(j)->GetField().at(1) + fld->GetFieldContainer().at(j)->GetField().at(1);
        fCollect.at(j)->SetField( Bx, By );
      }
    }
  }

}

void XFieldHandle :: Run()
{
  if ( is_exist(fTarget) )
    calfield(fTarget);
}


std::vector<Quench::XFieldContainer*> XFieldHandle :: GetFieldCollection(const std::string& name)
{
  fTarget = name;

  return fieldcollection(fTarget);
}


Quench::XFieldContainer* XFieldHandle :: GetFieldEntry(const std::string& name, const int id)
{
  fTarget = name;

  return fieldentry(fTarget, id);
}


Quench::XFieldContainer* XFieldHandle :: GetFieldEntry(const int z, const int r)
{
  const int mz = GetInfoEntry(fTarget)->GetMesh()[1];
  const int id = r*mz + z;

  return fieldentry(fTarget, id);
}


Quench::XFieldContainer* XFieldHandle :: fieldentry(const std::string& name, const int id)
{
  //if ( is_exist(name) && fCollect.size()==0 )
  //  calfield(name);

  return fCollect.at(id);
}

std::vector<Quench::XFieldContainer*> XFieldHandle :: fieldcollection(const std::string& name)
{
  //if ( is_exist(name) && fCollect.size()==0 )
  //  calfield(name);

  return fCollect;
}
