/**
 *  @file   XCoilConductor.hpp
 *  @author Y.Yang (Kyushu University)
 *  @date   14 Aug 2016
 **/

#ifndef XCoilConductor_HH
#define XCoilConductor_HH

#include "XCoilBase.hpp"

namespace Quench
{ class XCoilConductor; }

/// abstract class to define the conductor parameters
//
class Quench::XCoilConductor : public Quench::XCoilBase
{
  public:
    /// @brief constructor
    XCoilConductor();

    /*! deconstructor */
    virtual ~XCoilConductor();

    /// @brief  setup the dimension of conductor.
    /// @detail in this case, the dimension of conductor is the dimension of the stabilizer
    virtual void SetDimension(const double lz, const double lr);

    /// @brief returns the array of dimension of conductor (stabilizer)
    virtual double* GetDimension() { return fSize; }
    virtual double  GetDimension(const Coil dim) const;

    /// @brief returns the total length along r and z direction
    virtual double GetTotalLength(const Coil dim) const;

    /// @brief setup the insulation size along the z and r direction respectively
    virtual void SetInsSize(const double lz, const double lr);

    /// @brief  return the array of insulation size.
    /// @return the first array is Lz, and the second array is Lr
    virtual double* GetInsSize() { return fTape; }
    virtual double  GetInsSize(const Coil dim) const;

    /// @brief  calculate the total area of conductor including the stabilizer and insulation
    /// @detail the total length along both z and r axis is calculated as \f$ l = l_{cdt} + 
    ///         2\times l_{ins} \f$
    virtual double GetTotalArea() const;

    /// @brief  calculate this conductor area only.
    /// @detail noting that the insulation is not included
    virtual double GetArea() const;

    /// @brief  calculate the area of only for the insulation
    /// @detail the insulation is calculated as \f$ A_{ins} = A_{tot} - A_{cdt} \f$
    virtual double GetInsArea() const;

    /// @brief returns the type
    virtual Geometry GetType() const { return fType; }


  private:
    double* fSize;
    double* fTape;
    const Geometry fType;
};

#endif
