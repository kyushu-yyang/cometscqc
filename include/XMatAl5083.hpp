/**
 *  @file   XMatAl5083.hpp
 *  @author Y.Yang (Kyushu University)
 *  @date   31 Oct 2016
 **/

#ifndef XMatAl5083_HH
#define XMatAl5083_HH

#include "XMaterial.hpp"

/// Aluminium 5083-O material library
//
class XMatAl5083 : public XMaterial
{
  public:
    /// constructor
    XMatAl5083();

    /// deconstructor
    virtual ~XMatAl5083();

    /// @brief return the material density
    virtual double GetDensity() const { return 2700.; }

    /// @brief return the material thermal conductivity
    virtual double GetConductivity();

    /// @brief return the material electric resistivity
    virtual double GetResistivity();

    /// @brief return the material heat capacity
    virtual double GetCapacity();

  protected:
    /// @brief get the thermal conductivity from nist library
    double GetNistConductivity(const double T);

    /// @brief get the heat capacity from nist library
    double GetNistCapacity(const double T);

};

#endif
