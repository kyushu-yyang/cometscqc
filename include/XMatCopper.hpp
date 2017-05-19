/**
 *  @file   XMatCopper.hpp
 *  @author Y.Yang (Kyushu University)
 *  @date   5 Aug 2016
 */

#ifndef XMatCopper_HH
#define XMatCopper_HH

#include "XMaterial.hpp"

/// class to handle the material of copper
//
class XMatCopper : public XMaterial
{
  public:
    /*! constructor */
    XMatCopper();

    /*! deconstructor */
    virtual ~XMatCopper();

    /*! @brief setupp field */
    virtual void SetField(const double B) { fFld = B; }

    /*! @brief return material density */
    virtual double GetDensity() const { return 8960.; }

    /*! @brief return material resistivity */
    virtual double GetResistivity();

    /*! @brief return material thermal conductivity */
    virtual double GetConductivity();

    /*! @brief return material heat capacity */
    virtual double GetCapacity();


  protected:
    /*! @brief calculate resistivity */
    double calresist(const double T) const;

    /*! @brief calculate magnetoresist */
    double calmagresist(double rho, double B) const;

    /*! @brief calculate thermal conductivity */
    double calconduct(const double T) const;

    /*! @brief calculate magneto-thermal conductivity */
    double calmagcdt(double k) const;

    /*! @brief calculate heat capacity */
    double calcapacity(const double T) const;

  private:
    double fRhoRT;

};

#endif
