/**
 *  @file   XMatKapton.hpp
 *  @author Y.Yang (Kyushu University)
 *  @date   6 Aug 2016
 */

#ifndef XMatKapton_HH
#define XMatKapton_HH

#include "XMaterial.hpp"

class XMatKapton : public XMaterial
{
  public:
    /*! constructor */
    XMatKapton();

    /*! deconstructor */
    ~XMatKapton() {}

    /*! return the heat capacity */
    virtual double GetCapacity();

    /*! return the thermal conductivity */
    virtual double GetConductivity();

    /*! @brief return the density */
    //virtual double GetDensity() const { return 1420.; }
    virtual double GetDensity() const { return 1740.; }

    /*! @brief return the resistivity */
    virtual double GetResistivity() { return 1.; }


  protected:
    /*! calculate thermal conductivity */
    void calconductivity(double &k);

    /*! calculate capacity */
    void calcapacity(double &C);
};

#endif
