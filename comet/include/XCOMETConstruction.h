/**
 * @file   XCOMETConstruction.h
 * @author Y.Yang (Kyushu University)
 * @date   7 Oct 2016
 **/

#ifndef XCOMETConstruction_HH

#include "XThermalSolver.hpp"
#include "XFieldHandle.hpp"
#include "XCoilHandle.hpp"
#include "XQuenchTransient.hpp"

using namespace Quench;

class XCOMETConstruction : public XQuenchTransient
{
  public:
    /// constructor
    XCOMETConstruction();

    /// deconstructor
    ~XCOMETConstruction();

    /// @brief set operation time
    void SetOperationTime(const double time) { fDay = time; }

    /// @brief calculate the coil total resistance
    double GetCoilResistance(XThermalSolver* solver);

    /// @brief construct cs0 coil
    void ConstructCS0(const std::string& radfile, const char* tempfile);

    /// @brief construct cs1 coil
    void ConstructCS1(const std::string& radfile, const char* tempfile);

    /// @brief construct ms1 coil
    void ConstructMS1(const std::string& radfile, const char* tempfile);

    /// @brief construct ms1 coil
    void ConstructMS2(const std::string& radfile, const char* tempfile);

    /// @brief update the resistance of conductor
    /// @detail please use this function after SetMaterial() and the current decay calculation
    void UpdateQuench(XThermalSolver* solve);

    /// @brief begin of run
    virtual void Begin();

    /// @brief run
    virtual void Run();

    /// @brief end of run
    virtual void End();

    /// @brief connect two magnets
    void ConnectMagnet(XThermalSolver* mag1, XThermalSolver* mag2, const double l);
    void ConnectShell(XThermalSolver* mag1, XThermalSolver* mag2, const double l);


  protected:
    /// @brief setup conductor size
    XCoilBase* GetConductor();

    /// @brief setup strip size
    XCoilBase* GetStrip(const double thick=1.*mm);

    /// @brief setup shell size
    XCoilBase* GetShell(const double thick=80.*mm);

    /// @brief construct magnetic field
    void ConstructField(XFieldHandle* fld);
    

  private:
    double fDay;
    XFieldHandle* fFld;
    XThermalSolver* fCS0;
    XThermalSolver* fCS1;
    XThermalSolver* fMS1;
    XThermalSolver* fMS2;
};

#endif
