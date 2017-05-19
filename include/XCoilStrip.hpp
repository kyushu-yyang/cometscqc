/**
 *  @file   XCoilStrip.hpp
 *  @author Y.Yang (Kyushu University)
 *  @date   15 Aug 2016
 */

#ifndef XCoilStrip_HH
#define XCoilStrip_HH

#include "XCoilBase.hpp"

namespace Quench
{ class XCoilStrip; }

/// abstract class to define the quench propagator strip parameters */
//
class Quench::XCoilStrip : public Quench::XCoilBase
{
  public:
    /*! constructor */
    XCoilStrip();

    /*! deconstructor */
    virtual ~XCoilStrip();

    /// @brief setup dimension of strip only
    virtual void SetDimension(const double lz, const double lr);

    /// @brief returns the array of dimension of strip only
    virtual double* GetDimension() { return fSize; }
    virtual double  GetDimension(const Coil dim) const;

    /// @brief returns the total length
    virtual double GetTotalLength(const Coil dim) const;

    /*! @brief setup the insulation thickness */
    virtual void SetInsSize(const double lz, const double lr);

    /*! @brief returns the insulation thickness */
    virtual double* GetInsSize() { return fIns; }
    virtual double  GetInsSize(const Coil dim) const;

    /*! @brief calculate the total area of strip */
    virtual double GetTotalArea() const;

    /*! @brief calculate the strip area only */
    virtual double GetArea() const;

    /*! @brief calculate the area of insulation for strip */
    virtual double GetInsArea() const;

    /// @brief return the enumeration of the type
    virtual Geometry GetType() const { return fType; }


  private:
    double* fSize;
    double* fIns;
    const Geometry fType;
};

#endif
