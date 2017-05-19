/**
 *  @file   XMatG10.hpp
 *  @author Y.Yang (Kyushu University)
 *  @date   31 Oct 2016
 **/

#ifndef XMatG10_HH
#define XMatG10_HH

#include <string>
#include "XMaterial.hpp"

/// G10 material library
//
class XMatG10 : public XMaterial
{
  public:
    /// constructor
    XMatG10();

    /// deconstructor
    virtual ~XMatG10();

    /// @brief set insulation direction
    void SetDirection(const std::string &dir);

    /// @brief return the material density
    virtual double GetDensity() const { return 1850.; }

    /// @brief return the material thermal conductivity
    virtual double GetConductivity();

    /// @brief return the material electric resistivity
    virtual double GetResistivity() { return 0.; }

    /// @brief return the material heat capacity
    virtual double GetCapacity();


  protected:
    /// @brief get thermal conductivity from nist library
    double GetNistConductivity(const double T);

    /// @brief get heat capacity from nist library
    double GetNistCapacity(const double T);


  private:
    std::string fDirection;

};

#endif
