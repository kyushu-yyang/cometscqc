/**
 *  @file   XCoilShell.hpp
 *  @author Y.Yang (Kyushu University)
 *  @date   16 Aug 2016
 *          18 Aug 2016 - modified by Y.Yang
 **/

#ifndef XCoilShell_HH
#define XCoilShell_HH

#include "XCoilBase.hpp"

namespace Quench
{ class XCoilShell; }

class Quench::XCoilShell : public Quench::XCoilBase
{
  public:
    /*! constructor */
    XCoilShell();

    /*! deconstructor */
    virtual ~XCoilShell();

    /*! @brief setup the dimension of shell */
    virtual void SetDimension(const double lz, const double lr);

    /*! @brief return the dimension of shell */
    virtual double* GetDimension() { return fSize; }
    virtual double  GetDimension(const Coil dim) const;

    /// @brief returns the total length
    virtual double GetTotalLength(const Coil dim) const;

    /*! @brief setup the insulation thickness */
    virtual void SetInsSize(const double lz, const double lr);

    /*! @brief return the insulation thickness */
    virtual double* GetInsSize() { return fIns; }
    virtual double  GetInsSize(const Coil dim) const;

    /*! @brief calculate the total area of shell */
    virtual double GetTotalArea() const;

    /*! @brief calculate the shell area only */
    virtual double GetArea() const;

    /*! @brief calculate the insulation area */
    virtual double GetInsArea() const;

    /// @brief return the type
    virtual Geometry GetType() const { return fType; }


  private:
    double* fSize;
    double* fIns;
    const Geometry fType;
};

#endif
