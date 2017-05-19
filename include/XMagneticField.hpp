/**
 *  @file   XMagneticField.hpp
 *  @author Y.Yang (Kyushu University)
 *  @date   3 Aug 2016
 */

#ifndef XMagneticField_HH
#define XMagneticField_HH

#include <vector>
#if __cplusplus <= 199711L
  #error This library needs C++11 at least
#endif
#include <array>

#ifndef IFdmUnits_HH
#include "IFdmUnits.hpp"
#endif

namespace Quench
{
  class XFieldContainer;
  class XMagneticField;
  class XBiotSavart;
}

/// class to contain the magnetic field data
//
class Quench::XFieldContainer
{
  public:
    /*! constructor */
    XFieldContainer() : fPos(NULL), fId(NULL) {}

    /*! deconstructor */
    ~XFieldContainer() {
      if (fPos)   delete [] fPos;
      if (fId)    delete [] fId;
    }

    /*! @brief setup id number */
    void SetId(const int idz, const int idr) {
      if (!fId)  fId = new int[2];
      fId[0] = idz;
      fId[1] = idr;
    }

    /*! @brief setup cell position */
    void SetPosition(const double z, const double r) {
      if (!fPos)  fPos = new double[2];
      fPos[0] = z;
      fPos[1] = r;
    }

    /*! @brief setup magnetic field */
    void SetField(const double Bz, const double Br) {
      fField.at(0) = Bz;
      fField.at(1) = Br;
    }

    /*! @brief setup vector potential */
    void SetPotential(const double vecA) { fVector = vecA; }

    /*! @brief return id number */
    int* GetId() const { return fId; }

    /*! @brief return cell position */
    double* GetPosition() const { return fPos; }

    /// @brief return magnetic field array
    std::array<double,2> GetField() const { return fField; }

    /*! @brief return magnetic vector potential */
    double GetPotential() const { return fVector; }

    
  private:
    double* fPos;
    int*    fId;
    double  fVector;
    std::array<double, 2> fField;
};


/// base class to handle magnetic field
//
class Quench::XMagneticField
{
  public:
    /*! deconstructor */
    virtual ~XMagneticField() {}

    /*! @brief setup solenoid size */
    virtual void SetSolenoid(const double z0, const double z1, const double r0, const double r1) = 0;

    /*! @brief setup solenoid mesh */
    virtual void SetSolenoidMesh(const int mz, const int mp, const int mr) = 0;

    /*! @brief setup calculation region */
    virtual void SetMapRange(const double z0, const double z1, const double r0, const double r1) = 0;

    /*! @brief setup regin mesh */
    virtual void SetMapMesh(const int mz, const int mr) = 0;

    /*! @brief setup current */
    virtual void SetCurrent(const double I) = 0;

    /*! @brief setup conductor size */
    virtual void SetConductorArea(const double A) = 0;

    /*! @brief run field calculation */
    virtual void Run() = 0;

    /*! @brief update the magnetic field */
    virtual Quench::XFieldContainer* GetFieldEntry(const int iz, const int jr) = 0;

    /*! @brief return field container vector */
    virtual std::vector<Quench::XFieldContainer*> GetFieldContainer() = 0;
};


/// class to calculate the magentic field
//
class Quench::XBiotSavart : public Quench::XMagneticField
{
  friend class XFieldHandle;

  public:
    /*! constructor */
    XBiotSavart();

    /*! deconstructor */
    virtual ~XBiotSavart();

    /*! @brief setup solenoid size */
    virtual void SetSolenoid(const double z0, const double z1, const double r0, const double r1);

    /*! @brief setup solenoid mesh */
    virtual void SetSolenoidMesh(const int mz, const int mp, const int mr);

    /*! @brief setup calculation region */
    virtual void SetMapRange(const double z0, const double z1, const double r0, const double r1);

    /*! @brief setup regin mesh */
    virtual void SetMapMesh(const int mz, const int mr);

    /*! @brief setup current */
    virtual void SetCurrent(const double I) { fCurrent = I / fCdtA; }

    /*! @brief setup conductor size */
    virtual void SetConductorArea(const double A) { fCdtA = A; }

    /*! @brief calculation field */
    virtual void Run() { calfield(); }

    /*! @brief return field container */
    virtual Quench::XFieldContainer* GetFieldEntry(const int iz, const int jr);

    /*! @brief return field container vector */
    virtual std::vector<Quench::XFieldContainer*> GetFieldContainer() { return fHC; }


  protected:
    /*! @brief calculate magnetic field */
    void calfield();

    /*! @brief calculate vector potential at this position */
    double calpotential(const double z, const double r) const;

  private:
    int*    fMapMesh;
    double* fMapRange;
    double  fCurrent;  /// current density
    double  fCdtA;
    int* fAMsh;       /// mesh to calculate potential
    double* fSolenoid;
    std::vector <Quench::XFieldContainer*> fHC;
};

#endif
