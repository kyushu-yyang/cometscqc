#include <iostream>
#include <cmath>
#include <map>

#include "XMatAluminium.hpp"
#include "XMatCopper.hpp"
#include "XMatNbTi.hpp"
#include "XMatKapton.hpp"
#include "XMatG10.hpp"
#include "XMatAl5083.hpp"

#include "XQuenchExcept.hpp"
#include "XQuenchLogger.hpp"
#include "XProcessManager.hpp"

using Quench::XQuenchLogger;
using Quench::XProcessManager;


XProcessManager :: XProcessManager()
    : XMeshLoop(), 
      fCoil(NULL),
      fInsFactor(1.),
      fName(""),
      fIc(14.14e+3)
{}


XProcessManager :: ~XProcessManager()
{
  if (fCoil)  delete fCoil;
}


void XProcessManager :: SetCoilHandler(XCoilHandle* handler)
{
  fCoil = handler;
  // update mesh
  const int z = fCoil->GetMesh(iZ);
  const int p = fCoil->GetMesh(iPhi);
  const int r = fCoil->GetMesh(iR);
  SetMesh(z, p, r);

  fName = fCoil->GetName();
}


void XProcessManager :: SetUniformField(const double fld)
{
  for (int k=1; k<fMshR+1; k++) {
    for (int j=1; j<fMshP+1; j++) {
      for (int i=1; i<fMshZ+1; i++)
        fMC.at( Id(i,j,k) )->SetField(fld);
    }
  }

  QuenchInfo( "set the uniform magnetic field: " << fld << " Tesla." );
}


void XProcessManager :: SetFieldHandler(XFieldHandle* hand)
{
  if (!hand) {
    QuenchError( XQuenchLogger::ERROR, "null field handler." );
    XQuenchExcept except("null field handler.");
    throw except;
  }
  
  if (!fCoil) {
    QuenchError( XQuenchLogger::ERROR, "no coil handler set." );
    XQuenchExcept except("please set the coil handler first.");
    throw except;
  }

  double Bz, Br, Btot;

// modified on 5.17
/*
  // fill the magnetic field into the container
  for (int k=1; k<fMshR+1; k++) {
    for (int j=1; j<fMshP+1; j++) {
      for (int i=1; i<fMshZ+1; i++) {
        if (fDC.at( Id(i,j,k) )->GetGeometry()!=kShell ) {
          Bz   = hand->GetFieldEntry(i-1, k-1)->GetField().at(0);
          Br   = hand->GetFieldEntry(i-1, k-1)->GetField().at(1);
          Btot = sqrt( pow(Bz,2) + pow(Br,2) );
          fMC.at( Id(i,j,k) )->SetField( Btot );
        }
      }
    }
  }
*/

  std::vector<int> conductor = fCoil->GetLayerId(kConductor);
  std::vector<int> strip     = fCoil->GetLayerId(kStrip);

  for (std::vector<int>::size_type k=0; k<conductor.size(); k++) {
    for (int j=1; j<fMshP+1; j++) {
      for (int i=1; i<fMshZ+1; i++) {
        Bz   = hand->GetFieldEntry(i-1, k)->GetField().at(0);
        Br   = hand->GetFieldEntry(i-1, k)->GetField().at(1);
        Btot = sqrt( pow(Bz,2) + pow(Br,2) );

        fMC.at( Id(i,j,conductor.at(k)) )->SetField( Btot );
      }
    }
  }

  for (std::vector<int>::size_type k=0; k<strip.size(); k++) {
    for (int j=1; j<fMshP+1; j++) {
      for (int i=1; i<fMshZ+1; i++) {
        Btot = fMC.at( Id(i,j,strip.at(k)+1) )->GetField();
        fMC.at( Id(i,j,strip.at(k)) )->SetField( Btot );
      }
    }
  } 

  QuenchInfo( "set the calculated magnetic field." );
}


void XProcessManager :: SetRadiationHandler(XRadiationHandle* hand)
{
  if (!hand) {
    QuenchError( XQuenchLogger::ERROR, "radiation handler is null." );
    XQuenchExcept except("radiation handler is null.");
    throw except;
  }

  const int mshz = hand->GetMesh(iZ);
  //const int mshp = hand->GetMesh(iPhi);
  const int mshr = hand->GetMesh(iR);

  const int radmshz = fMshZ / mshz;
  //const int radmshp = fMshP / mshp;
  const int radmshr = fMshR / mshr;

  QuenchError( XQuenchLogger::INFO, "check radiation input file mesh -> z: " 
                                    << hand->GetMesh(iZ)   << " phi: "
                                    << hand->GetMesh(iPhi) << " r: " 
                                    << hand->GetMesh(iR) );

  int iz = 0; int jp = 0; int kr = hand->GetMesh(iR);
  int id = 0;
  double fluence, RRR, deposit;

  std::vector<int> conduct = fCoil->GetLayerId(kConductor);
  std::vector<int> strip = fCoil->GetLayerId(kStrip);
  std::vector<int> shell = fCoil->GetLayerId(kShell);

/*
  for (int k=0; k<fMshR+2; k++) {
    //kr = mshr - static_cast<int>(k/radmshr);   // modified to reverse the r direction for loading radiation
    kr = static_cast<int>((k-1)/radmshr);
    if (kr>=mshr)
      kr = mshr-1;
    if (kr<=0)
      kr = 0;
*/
  for (std::vector<int>::size_type k=0; k<conduct.size(); k++) {
    kr = conduct.at(k);

    for (int j=0; j<fMshP+2; j++) {
      // only for the case when phi mesh is equal to the rad phi size
      if (j>0 && j<fMshP+1) jp = j-1;
      else jp = 0;

      for (int i=0; i<fMshZ+2; i++) {
        iz = static_cast<int>(i/radmshz);
        if (iz<=0)
          iz = 0;
        if (iz>=mshz)
          iz = mshz-1;
        // get local id for radiation
        id = hand->Id( iz, jp, k );
        // get neutron fluence and deposit
        fluence = hand->GetEntry(id)->GetNeutron();
        deposit = hand->GetEntry(id)->GetDose();
        // convert neutron fluence to RRR 
        //RRR = hand->GetRRR( fDC.at(Id(i,j,k))->GetGeometry(), fluence );
        RRR = hand->GetRRR( kConductor, fluence );
        // fill RRR value into container
        fMC.at( Id(i,j,conduct.at(k)) )->SetDeposit( deposit );
        fMC.at( Id(i,j,conduct.at(k)) )->SetRRR(RRR);
        fMC.at( Id(i,j,conduct.at(k)) )->SetStack(fluence);
      }
    }
  }

// modified on 2017.05.16
  for (std::vector<int>::size_type k=0; k<strip.size(); k++) {
    for (int j=0; j<fMshP+2; j++) {
      for (int i=0; i<fMshZ+2; i++) {
        // get the neutron fluence from next layer
        fluence = fMC.at( Id(i,j,strip.at(k)+1) )->GetStack();
        RRR = hand->GetRRR( kStrip, fluence );
        deposit = fMC.at( Id(i,j,strip.at(k)+1) )->GetDeposit();
        // fill RRR value into container
        fMC.at( Id(i,j,strip.at(k)) )->SetDeposit( deposit );
        fMC.at( Id(i,j,strip.at(k)) )->SetRRR(RRR);
      }
    }
  }

  for (std::vector<int>::size_type k=0; k<shell.size(); k++) {
    for (int j=0; j<fMshP+2; j++) {
      for (int i=0; i<fMshZ+2; i++) {
        // get the neutron fluence from next layer
        //RRR = fMC.at( Id(i,j,shell.at(k)-1) )->GetRRR();
        RRR = 0.;
        deposit = fMC.at( Id(i,j,shell.at(k)-1) )->GetDeposit();
        // fill RRR value into container
        fMC.at( Id(i,j,shell.at(k)) )->SetDeposit( deposit );
        fMC.at( Id(i,j,shell.at(k)) )->SetRRR(RRR);
      }
    }
  }
  
}


void XProcessManager :: SetInitialTemperature(XInitialTemperature* temp)
{
  int id = 0;
  int cnt = 0;
  std::vector<XInitTempContainer*> Tint = temp->GetContainer();

  for (int k=0; k<fMshR+2; k++) {
    for (int j=0; j<fMshP+2; j++) {
      for (int i=0; i<fMshZ+2; i++) {
        if (k>0 && k<fMshR+1 && j>0 && j<fMshP+1 && i>0 && i<fMshZ+1) {
          id = Id(i, j, k);
          if (Tint.at(cnt)->GetNodeId()!=id) {
            QuenchFatal("id from initial temperature profile " << Tint.at(cnt)->GetNodeId() << " is not equal this id " << id << ".");
            XQuenchExcept except("id from initial temperature profile is not equal.");
            throw except;
          }
          fMC.at( id )->SetTemperature( Tint.at(cnt)->GetTemperature() );
          cnt++;
        }
      }
    }
  }
}


void XProcessManager :: SetUniformRRR(const Geometry part, const double RRR)
{
  std::vector<int> id = fCoil->GetLayerId(part);

  if ( id.size()==0 ) {
    QuenchError( XQuenchLogger::WARNING, "layer id size is zero." );
    throw;
  }

  int kr = 0;
  for (std::vector<int>::size_type k=0; k<id.size(); k++) {
    kr = id.at(k);
    for (int j=0; j<fMshP+2; j++) {
      for (int i=0; i<fMshZ+2; i++) 
        fMC.at( Id(i,j,kr) )->SetRRR( RRR );
    }
  }

  QuenchInfo( "set the uniform RRR value: " << RRR << " for " << fCoil->GetGeometryName(part) );
}


void XProcessManager :: SetUniformHeatGen(const double gen)
{
  for (int k=1; k<fMshR+1; k++) {
    for (int j=1; j<fMshP+1; j++) {
      for (int i=1; i<fMshZ+1; i++)
        fMC.at( Id(i,j,k) )->SetDeposit(gen);
    }
  }

  QuenchInfo( "set the uniform heat generation: " << gen );
}

void XProcessManager :: Initialize()
{
  if (!fCoil) {
    QuenchError( XQuenchLogger::ERROR, "coil handler is not set." );
    XQuenchExcept except("please set the coil handler first.");
    throw except;
  }

  const double T = 4.5;

  // allocate the data saving place
  init();
  InitTemp(T);
  InitPosition();
}


void XProcessManager :: Initialize(XCoilHandle* coil, XFieldHandle* fld)
{
  const double T = 4.5;

  SetCoilHandler(coil);
  init();
  InitTemp(T);
  InitPosition();
  SetFieldHandler(fld);
}


void XProcessManager :: ModifyGeometry(const int z, const int phi, const int r, const Geometry geo)
{
  const int id = Id(z, phi, r);
  fDC.at(id)->SetGeometry(geo);

  QuenchError( XQuenchLogger::DEBUG, "cell at (" << z << "," << phi << "," << r << ") is modified to geometry: " << geo );
}


void XProcessManager :: SetMaterial()
{
  if (fMC.size()<3) {
    QuenchError( XQuenchLogger::ERROR, "container is not set." );
    XQuenchExcept except("please check the handler setup.");
    throw except;
  }

  double T, RRR, B;
  Geometry geo;
  int id;

  for (int k=0; k<fMshR+2; k++) {
    for (int j=0; j<fMshP+2; j++) {
      for (int i=0; i<fMshZ+2; i++) {
        id  = Id(i,j,k);
        T   = fMC.at(id)->GetTemperature();
        RRR = fMC.at(id)->GetRRR();
        B   = fMC.at(id)->GetField();
        geo = fDC.at(id)->GetGeometry();

        if ( (i>0 && i<fMshZ+1) && (j>0 && j<fMshP+1) && (k>0 && k<fMshR+1) ) {
          switch (geo) {
            case kConductor: 
              fMC.at(id)->SetHeat( fMC.at(id)->GetDeposit()*4000. );
              SetConductorMat( id, fCoil->GetCoilType(k), T, RRR, B ); 
              break;
            case kStrip: 
              fMC.at(id)->SetHeat( fMC.at(id)->GetDeposit()*2700. );
              SetStripMat( id, fCoil->GetCoilType(k), T, RRR, B ); 
              break;
            case kShell: 
              fMC.at(id)->SetHeat( fMC.at(id)->GetDeposit()*2700. );
              SetShellMat( id, fCoil->GetCoilType(k), T, RRR, B );
              break;
            case kG10:
              fMC.at(id)->SetHeat( fMC.at(id)->GetDeposit()*1000. );
              SetG10Mat( id, T );
              break;
            case kA5083:
              fMC.at(id)->SetHeat( fMC.at(id)->GetDeposit()*2700. );
              SetA5083Mat( id, T );
              break;
            default:
              QuenchError( XQuenchLogger::WARNING, "geometry " << fCoil->GetGeometryName(geo) << " did not exist." );
              break;
          }
        }
        
        //fMC.at(id)->SetHeat( fMC.at(id)->GetDeposit()*4000. );
      }
    }
  }

}


size_t XProcessManager :: GetEntries() const
{
  if (fMC.size()!=fDC.size()) {
    QuenchError( XQuenchLogger::ERROR, "material container size is equal to the size of dimension container." );
    XQuenchExcept except("material container size is equal to the size of dimension container.");
    throw except;
  }

  return fMC.size();
}


void XProcessManager :: init()
{
  // initialize the container vector
  for (int k=0; k<fMshR+2; k++) {
    for (int j=0; j<fMshP+2; j++) {
      for (int i=0; i<fMshZ+2; i++) {
        fMC.push_back( new XMaterialInfo() );
        fDC.push_back( new XDimensionInfo );
      }
    }
  }

  // initialize id number
  int id = 0;
  for (int k=0; k<fMshR+2; k++) {
    for (int j=0; j<fMshP+2; j++) {
      for (int i=0; i<fMshZ+2; i++) {
        id = Id(i,j,k);
        fDC.at(id)->SetId(i, j, k);
        fDC.at(id)->SetNodeId( Id(i,j,k) );
        fDC.at(id)->SetPosition(0., 0., 0.);
        fDC.at(id)->SetPrePosition(0., 0., 0.);
        fDC.at(id)->SetPostPosition(0., 0., 0.);
        fDC.at(id)->SetCellSize(1., 1., 1.);
        fDC.at(id)->SetGeometry(kStrip);
      }
    }
  }

  // initialize geometry
  std::map<const int, const Geometry> coil = fCoil->GetCoilLayout();
  int k;  Geometry geo;

  for ( std::map<const int, const Geometry>::const_iterator it=coil.begin(); it!=coil.end(); ++it ) {
    k   = it->first;
    geo = it->second;
    for (int j=0; j<fMshP+2; j++) {
      for (int i=0; i<fMshZ+2; i++) {
        fDC.at( Id(i,j,k) )->SetGeometry(geo);
      }
    }
  }

}


void XProcessManager :: InitTemp(const double T)
{
  for (int k=0; k<fMshR+2; k++) {
    for (int j=0; j<fMshP+2; j++) {
      for (int i=0; i<fMshZ+2; i++) 
        fMC.at( Id(i,j,k) )->SetTemperature(T);
    }
  }
}


void XProcessManager :: InitPosition()
{
  if ( fMshR!=fCoil->GetLayoutEntries() ) {
    QuenchError( XQuenchLogger::ERROR, "total layout: " << fCoil->GetLayoutEntries()
                                       << ", layout " << fMshR << " is needed." );
    XQuenchExcept except("layout is not enough.");
    throw except;
  }

  double preZ  = 0.;  double preP  = 0.;  double preR  = 0.;        // pre-position
  double postZ = 0.;  double postP = 0.;  double postR = 0.;        // post-position

  double z  = 0.;  double p  = 0.;  double r  = 0.;        // position
  double lz = 0.;  double lp = 0.;  double lr = 0.;        // length

  int id = 0;

  // length edge of strip along direction z
  double edge = fCoil->GetStripEdge();
  // approach z factor
  const double factor = fCoil->GetApproachZ();
  QuenchError( XQuenchLogger::INFO, "factor along the z direction: " << factor );
  // length of cell
  lz = fCoil->GetCoilType(2)->GetTotalLength(iZ) * factor;
  lp = fCoil->GetCoilSize(iPhi) / fMshP;

  for (int k=1; k<fMshR+1; k++) {
    //lr = fCoil->GetCoilLayout(k)->GetTotalLength(iR);
    lr = fCoil->GetCoilType(k)->GetTotalLength(iR);
    r += lr/2.;
    //postR = k==fMshR ? r+lr/2. : r+lr/2.+fCoil->GetCoilLayout(k+1)->GetTotalLength(iR)/2.;
    postR = k==fMshR ? r+lr/2. : r+lr/2.+fCoil->GetCoilType(k+1)->GetTotalLength(iR)/2.;

    QuenchError( XQuenchLogger::CONFIG, "id: " << k << ", r position:" << r/mm << " [mm]"
                                        << " , cell size_r: " << lr/mm << " [mm]");

    p  = 0.;
    preP = 0.;

    for (int j=1; j<fMshP+1; j++) {
      p += lp/2.;
      postP = p+lp;

      z  = 0.;
      preZ = 0.;
      for (int i=1; i<fMshZ+1; i++) {
        id = Id(i,j,k);

        z += lz/2.;
        postZ = i==fMshZ ? z+lz/2. : z+lz;

        fDC.at(id)->SetPrePosition( preZ, preP, preR );
        fDC.at(id)->SetPosition( z, p, r );
        fDC.at(id)->SetPostPosition( postZ, postP, postR );
        fDC.at(id)->SetCellSize( lz, lp, lr );

        // re-setup the pre-position for the cell at the boundary
        //if ( fDC.at(id)->GetGeometry()==kStrip && i==1 ) {
        if ( i==1 ) {
          if ( fCoil->GetCoolingConfigure(k)==kSide || fCoil->GetCoolingConfigure(k)==kLeft ) {
            edge = fCoil->GetStripEdgeSize(k);
            fDC.at(id)->SetPrePosition( preZ-edge, preP, preR );
            fDC.at(id)->SetCellSize( lz+edge, lp, lr );
          }
        }

        // re-setup the post-position for the cell at the boundary
        //if ( fDC.at(id)->GetGeometry()==kStrip && i==fMshZ ) {
        if ( i==fMshZ ) {
          if ( fCoil->GetCoolingConfigure(k)==kSide || fCoil->GetCoolingConfigure(k)==kRight ) {
            edge = fCoil->GetStripEdgeSize(k);
            fDC.at(id)->SetPostPosition( postZ+edge, postP, postR );
            fDC.at(id)->SetCellSize( lz+edge, lp, lr );
          }
        }

        preZ = z;
        z += lz/2.;
      }

      preP = p;
      p += lp/2.;
    }

    preR = r;
    r += lr/2.;
  }

}


void XProcessManager :: SetConductorMat(const int id, const XCoilBase* cdt, const double T, const double RRR, const double B)
{
  XMatCopper    cu;
  XMatAluminium al;
  XMatNbTi      sc;
  //XMatKapton    kap;
  XMatG10       kap;

  sc.SetIcAt5Tesla(fIc);

  cu.SetMaterialProperty(T, 50., B);     // copper RRR: 50
  al.SetMaterialProperty(T, RRR, B);
  sc.SetMaterialProperty(T, RRR, B);
  kap.SetMaterialProperty(T, RRR, B);

  // get material ratio
  const double ratio_Al = fCoil->GetMaterialRatio(iAluminium);
  const double ratio_Cu = fCoil->GetMaterialRatio(iCopper);
  const double ratio_Sc = fCoil->GetMaterialRatio(iNbTi);

  // get density
  const double rho_Cu = cu.GetDensity();
  const double rho_Al = al.GetDensity();
  const double rho_Sc = sc.GetDensity();
  const double rho_avg = rho_Cu*ratio_Cu + rho_Al*ratio_Al + rho_Sc*ratio_Sc;
  //const double rho_avg = 4000.;

  fMC.at(id)->SetDensity( rho_avg );

  // calculate average capacity for conductor
  const double C_Al = ratio_Al * rho_Al * al.GetCapacity() / rho_avg;
  const double C_Cu = ratio_Cu * rho_Cu * cu.GetCapacity() / rho_avg;
  const double C_Sc = ratio_Sc * rho_Sc * sc.GetCapacity() / rho_avg;
  const double C_avg = C_Al + C_Cu + C_Sc;

  fMC.at(id)->SetCapacity( C_avg );

  // calculate average thermal conductivity
  //const double k_Al = al.GetConductivity();
  const double k_Al = al.hust_eq_therm(T,RRR,B);
  const double k_ins = kap.GetConductivity() * fInsFactor;

  //const double lz_ins = 2. * fCoil->GetCoilParts(kConductor)->GetInsSize(iZ);
  const double lz_ins = 2. * cdt->GetInsSize(iZ);
  const double lr_ins = 2. * cdt->GetInsSize(iR);
  const double lz_cdt = cdt->GetDimension(iZ);
  const double lr_cdt = cdt->GetDimension(iR);

  const double A_cdt  = cdt->GetArea();
  const double A_ins  = cdt->GetInsArea();
  
  const double kz = al.GetSerialk( lz_ins, k_ins, lz_cdt, k_Al );
  const double kr = al.GetSerialk( lr_ins, k_ins, lr_cdt, k_Al );
  const double kp = al.GetParallelk( A_ins, k_ins, A_cdt, k_Al );

  fMC.at(id)->SetConductivity( kz, kp, kr );

  // calculate average resistance
  /*
  const double dl_phi = fCoil->GetCoilSize(iPhi) / fMshP;
  const double A_Cu = A_cdt * ratio_Cu;
  const double A_Al = A_cdt * ratio_Al;

  const double R_Al = al.GetResistivity() * dl_phi / A_Al;
  const double R_Cu = cu.GetResistivity() * dl_phi / A_Cu;
  const double R_avg = pow( (1./R_Al + 1./R_Cu), -1. );

  if ( fMC.at(id)->GetStatus()!=kSuperconduct )
    fMC.at(id)->SetResistance( R_avg );
  */
}


void XProcessManager :: SetStripMat(const int id, const XCoilBase* strip, const double T, const double RRR, const double B)
{
  XMatAluminium al;
  //XMatKapton    kap;
  XMatG10    kap;

  al.SetMaterialProperty(T, RRR, B);
  kap.SetMaterialProperty(T, RRR, B);

  // setup density
  fMC.at(id)->SetDensity( al.GetDensity() );

  // setup heat capacity
  fMC.at(id)->SetCapacity( al.GetCapacity() );

  // setup thermal conductivity
  //const double k_Al  = al.GetConductivity();
  const double k_Al  = al.hust_eq_therm(T, RRR, B);
  const double k_ins = kap.GetConductivity() * fInsFactor;

  //const double lr_ins = 2. * fCoil->GetCoilParts(kStrip)->GetInsSize(iR);
  // both aluminium and insulation are using the total length
  const double lr_ins = 2. * strip->GetInsSize(iR);
  //const double lr_Al  = fCoil->GetCoilParts(kStrip)->GetDimension(iR);
  const double lr_Al  = strip->GetDimension(iR);

  const double kz = al.GetParallelk( lr_Al, k_Al, lr_ins, k_ins );
  const double kr = al.GetSerialk( lr_ins, k_ins, lr_Al, k_Al );
  const double kp = k_ins;

  fMC.at(id)->SetConductivity( kz, kp, kr );
  fMC.at(id)->SetResistance(0.);
}


void XProcessManager :: SetShellMat(const int id, const XCoilBase* shell, const double T, const double RRR, const double B)
{
  //XMatAluminium al;
  //XMatKapton    kap;
  XMatAl5083 al;
  XMatG10 ins;
  //XMatKapton ins;

  //al.SetMaterialProperty(T, RRR, B);
  al.SetTemperature(T);
  ins.SetTemperature(T);

  // setup density
  fMC.at(id)->SetDensity( al.GetDensity() );

  // setup heat capacity
  fMC.at(id)->SetCapacity( al.GetCapacity() );

  // setup thermal conductivity
  const double k_Al  = al.GetConductivity();
  const double k_ins = ins.GetConductivity();
  //const double k_ins = ins.GetConductivity() * fInsFactor;

  const double lr_ins = shell->GetInsSize(iR);
  const double lr_Al  = shell->GetDimension(iR);

  const double kz = k_Al;
  const double kr = al.GetSerialk( lr_ins, k_ins, lr_Al, k_Al );
  const double kp = k_Al;

  //fMC.at(id)->SetStack(k_Al);
  fMC.at(id)->SetConductivity( kz, kp, kr );
  fMC.at(id)->SetResistance(0.);
}

void XProcessManager :: SetG10Mat(const int id, const double T)
{
  XMatG10 g10;
  g10.SetTemperature(T);
  const double dens = g10.GetDensity();
  const double C    = g10.GetCapacity();
  const double k    = g10.GetConductivity();
  
  // fill data into entry
  fMC.at(id)->SetDensity( dens );
  fMC.at(id)->SetCapacity( C );
  fMC.at(id)->SetConductivity( k, k, k );
  fMC.at(id)->SetResistance(0.);
}

void XProcessManager :: SetA5083Mat(const int id, const double T)
{
  XMatAl5083 al;
  al.SetTemperature(T);
  const double dens = al.GetDensity();
  const double C    = al.GetCapacity();
  const double k    = al.GetConductivity();

  // fill data
  fMC.at(id)->SetDensity( dens );
  fMC.at(id)->SetCapacity( C );
  fMC.at(id)->SetConductivity( k, k, k );
  fMC.at(id)->SetResistance(0.);
}
