/**
 *  @file   XMatNbTi.hpp
 *  @author Y.Yang (Kyushu University)
 *  @date   12 Aug 2016
 */

#ifndef XMatNbTi_HH
#define XMatNbTi_HH

#include "XMaterial.hpp"

/// class for superconducting material: NbTi
//
class XMatNbTi : public XMaterial
{
  public:
    /*! constructor */
    XMatNbTi();

    /*! deconstructor */
    virtual ~XMatNbTi() {}

    /*! @brief return the density */
    virtual double GetDensity() const { return 6538.; }

    /*! @brief return the resistivity */
    virtual double GetResistivity() { return 1.; }

    /*! @brief return the conductivity */
    virtual double GetConductivity() { return 1.; }

    /*! @brief return the heat capacity */
    virtual double GetCapacity();

    /*! @brief return the critial current */
    virtual double GetCriticalI();

    /*! @brief return the critical temperature */
    virtual double GetCriticalT();

    /*! @brief return the current sharing temperature */
    virtual double GetSharingT(const double I);


  protected:
    /*! @brief calculate the heat capacity */
    double calcapacity(const double T, const double B) const;

    /*! @brief calculate the critical current */
    double calcriticalcurrent(const double T, const double B, const double I0=3000.) const;

    /*! @brief calculate the critical temperature */
    double calTc(const double B) const;

};

#endif
