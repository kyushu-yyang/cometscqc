/**
 *  @file   XFieldHandle.hpp
 *  @author Y.Yang (Kyushu University)
 *  @date   3 Aug 2016
 */

#ifndef XFieldHandle_HH
#define XFieldHandle_HH

#include <vector>
#include <string>

#ifndef XMagneticField_HH
#include "XMagneticField.hpp"
#endif

namespace Quench
{ 
  class XMagnetInfoContainer;
  class XFieldHandle; 
}

/// class to contain the magnet information
//
class Quench::XMagnetInfoContainer
{
  public:
    /// @brief constructor
    XMagnetInfoContainer() : fName(""), fSolenoid(NULL), fMesh(NULL) {}

    /*! deconstructor */
    ~XMagnetInfoContainer() {
      if (fSolenoid)  delete[] fSolenoid;
      if (fMesh)      delete[] fMesh;
    }

    /*! save solenoid name */
    void SetName(const std::string& name) { fName = name; }

    /*! save solenoid dimension */
    void SetDimension(const double z0, const double z1, const double r0, const double r1) {
      if (!fSolenoid)  fSolenoid = new double[4];
      fSolenoid[0] = r0;
      fSolenoid[1] = r1;
      fSolenoid[2] = z0;
      fSolenoid[3] = z1;
    }

    /*! save solenoid mesh */
    void SetMesh(const int mz, const int mr) {
      if (!fMesh)  fMesh = new int[2];
      fMesh[0] = mr;
      fMesh[1] = mz;
    }

    /*! return solenoid name */
    inline std::string GetName() const { return fName; }

    /*! return solenoid dimension */
    inline double* GetDimension() const { return fSolenoid; }

    /*! return solenoid mesh */
    inline int* GetMesh() const { return fMesh; }

  private:
    std::string fName;
    double* fSolenoid;
    int*    fMesh;
};


/// class to handle the result of magnetic field
//
class Quench::XFieldHandle
{
  public:
    /*! constructor */
    XFieldHandle();

    /*! constructor */
    XFieldHandle(const std::string& name);

    /*! deconstructor */
    ~XFieldHandle();

    /*! @brief set target magnet */
    void SetTarget(const std::string& name) { fTarget = name; }

    /*! @brief get target magnet */
    std::string GetTarget() const { return fTarget; }

    /*! @brief check this magnet exists or not */
    bool is_exist(const std::string& name) const;

    /*! @brief setup current */
    void SetCurrent(const double I) { fCurrent = I; }

    /*! @brief add coil */
    void AddCoil(const std::string& name, const double z0, const double z1, 
                                          const double r0, const double r1);

    /*! @brief add coil */
    void AddCoil(Quench::XMagnetInfoContainer* mag);

    /*! @brief setup solenoid mesh */
    void SetMesh(const std::string& name, const int mz, const int mr);

    /*! @brief run field simulation */
    void Run();

    /*! @brief return the magnet info container */
    std::vector<Quench::XMagnetInfoContainer*> GetInfoCollection() { return fHC; }

    /*! @brief return the magnet with this name */
    Quench::XMagnetInfoContainer* GetInfoEntry(const std::string& name);

    /*! @brief return the field container */
    std::vector<Quench::XFieldContainer*> GetFieldCollection(const std::string& name);
    
    /*! @brief return the field container */
    std::vector<Quench::XFieldContainer*> GetFieldCollection() { return fieldcollection(fTarget); }

    /*! @brief return the field container with the input id */
    Quench::XFieldContainer* GetFieldEntry(const std::string& name, const int id);
    
    /*! @brief return the field container with the input id */
    Quench::XFieldContainer* GetFieldEntry(const int id) { return fieldentry(fTarget, id); }

    /*! @brief return the field entry */
    Quench::XFieldContainer* GetFieldEntry(const int z, const int r);


  protected:
    /*! @brief calculate field map */
    void calfield(const std::string& name);
    
    /*! @brief return id-th entry of field container */
    Quench::XFieldContainer* fieldentry(const std::string& name, const int id);

    /*! @brief return the field collection */
    std::vector<Quench::XFieldContainer*> fieldcollection(const std::string& name);

  private:
    std::string fTarget;
    double fCurrent;
    std::vector<Quench::XMagnetInfoContainer*> fHC;
    std::vector<Quench::XFieldContainer*> fCollect;
};

#endif
