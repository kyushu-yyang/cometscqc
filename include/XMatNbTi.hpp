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

    /*! setup the critical current at 5 Tesla, 4.2 K */
    void SetIcAt5Tesla(const double Ic);

    /*! setup the Ic parameter */
    void SetIcParameter(const int par=1);

    /*! return the Ic parameter */
    int GetIcParameter() const { return fPar; }
    
    /*! return the critical current at 5 Tesla, 4.2 K */
    double GetIcAt5Tesla() const { return fIc5; }

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
    virtual double GetSharingT(const double I, const double T);

    /*! @brief return the fitting parameter */
    void GetIcPar(double &Tc0, double &Bc20, double &C0, double &alpha, double &beta, double &gamma);


  protected:
    /*! @brief calculate the heat capacity */
    double calcapacity(const double T, const double B) const;

    /*! @brief calculate the critical current */
    double calcriticalcurrent(const double T, const double B, const double I0=14.2e+3) const;

    /*! @brief calculate the critical temperature */
    double calTc(const double B) const;


  private:
    double fIc5;
    int    fPar;
};

#endif
