/**
 *  @file   XProcessManager.hpp
 *  @author Y.Yang (Kyushu University)
 *  @date   18 Aug 2016
 **/

#ifndef XProcessManager_HH
#define XProcessManager_HH

#include <vector>
#include <string>
#include "XQuenchInfo.hpp"
#include "XMeshLoop.hpp"
#include "XCoilHandle.hpp"
#include "XRadiationHandle.hpp"
#include "XFieldHandle.hpp"
#include "XInitialTemperature.hpp"

namespace Quench
{ class XProcessManager; }

/// class to manage the data initialization and update
//
class Quench::XProcessManager : public XMeshLoop
{
  friend class XThermalSolver;

  public:
    /// @brief constructor
    XProcessManager();

    /// @brief deconstructor
    virtual ~XProcessManager();

    /// @brief setup the class contained coil information
    void SetCoilHandler(XCoilHandle* handler);

    /// @brief setup the field handler
    void SetFieldHandler(XFieldHandle* hand);

    /// @brief setup radiation handler
    void SetRadiationHandler(XRadiationHandle* hand);

    /// @brief setup initial temperature handler
    void SetInitialTemperature(XInitialTemperature* temp);

    /// @brief return the coil handler
    XCoilHandle* GetCoilHandler() { return fCoil; }

    /// @brief setup uniform RRR for strip
    /// @param part enumeration of geometry: kShell/kConductor/kStrip
    void SetUniformRRR(const Geometry part, const double RRR);

    /// @brief setup the uniform field
    void SetUniformField(const double fld);

    /// @brief setup the uniform heat generation
    void SetUniformHeatGen(const double gen);

    /// @brief multiply a factor for insulation tape thermal conductivity
    void SetInsFactor(const double factor=1.) { fInsFactor = factor; }

    /// @brief return the factor for the insulation tape thermal conductivity
    double GetInsFactor() const { return fInsFactor; }

    /// @brief initialization
    void Initialize();

    /// @brief initialization
    void Initialize(XCoilHandle* coil, XFieldHandle* fld);

    /// @brief modify the material property
    void ModifyGeometry(const int z, const int phi, const int r, const Geometry geo);

    /// @brief initialize the material information
    void SetMaterial();

    /// @brief returns the dimensional container
    std::vector<XDimensionInfo*> GetDimensionContainer() { return fDC; }

    /// @brief returns the size of dimensional container
    size_t GetDimensionEntries() const { return fDC.size(); }

    /// @brief returns the dimensional container's entry
    XDimensionInfo* GetDimensionEntry(const int entry) { return fDC.at(entry); }

    /// @brief returns the material container
    std::vector<XMaterialInfo*> GetMaterialContainer() { return fMC; }

    /// @brief returns the size of material container
    size_t GetMaterialEntries() const { return fMC.size(); }

    /// @brief return this material container
    XMaterialInfo* GetMaterialEntry(const int entry) { return fMC.at(entry); }

    /// @brief returns the entries of both dimension and material
    size_t GetEntries() const;


  protected:
    /// @brief initialize the container vector
    void init();

    /// @brief initialize the temperature
    void InitTemp(const double T);

    /// @brief calculate the position and fill it into container
    void InitPosition();

    /// @brief fill conductor/strip/shell material property
    void SetConductorMat(const int id, const XCoilBase* cdt, const double T, const double RRR, const double B);
    void SetStripMat(const int id, const XCoilBase* strip, const double T, const double RRR, const double B);
    void SetShellMat(const int id, const XCoilBase* shell, const double T, const double RRR, const double B);
    void SetG10Mat(const int id, const double T);
    void SetA5083Mat(const int id, const double T);


  protected:
    XCoilHandle* fCoil;
    std::vector<XDimensionInfo*> fDC;
    std::vector<XMaterialInfo*> fMC;
    double fInsFactor;

  private:
    std::string fName;
};

#endif
