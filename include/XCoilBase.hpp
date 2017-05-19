/**
 *  @file   XCoilBase.hpp
 *  @author Y.Yang (Kyushu University)
 *  @date   18 Aug 2016
 **/

#ifndef XCoilBase_HH
#define XCoilBase_HH

#include <string>

/// @eum coil geometry
enum Geometry
{
  kConductor,    /// conductor
  kStrip,        /// Al strip
  kShell,        /// shell structure
  kG10,          /// G10
  kA5083         /// A5083
};

/// @enum dimensional flag
enum Coil
{
  iZ = 0,        /// z axis
  iPhi = 1,      /// phi axis
  iR = 2         /// r axis
};

namespace Quench
{ class XCoilBase; }

/// base class to handle the info of strip/conductor/shell
//
class Quench::XCoilBase
{
  public:
    /// @brief deconstructor
    virtual ~XCoilBase() {}

    /// @brief setup the dimension for conductor/strip/shell
    virtual void SetDimension(const double lz, const double lr) = 0;

    /// @brief returns the dimension of conductor/strip/shell
    /// @param dim enumeration Coil::iZ/iR, noting iPhi is disable
    virtual double  GetDimension(const Coil dim) const = 0;
    virtual double* GetDimension() = 0;

    /// @brief returns the total length along z and r direction
    virtual double GetTotalLength(const Coil dim) const = 0;

    /// @brief setup the size of insulation for conductor/strip/shell
    virtual void SetInsSize(const double lz, const double lr) = 0;

    /// @brief returns the size of insulation for conductor/strip/shell
    /// @param dim enumeration Coil::iZ/iR, noting iPhi is disable
    virtual double  GetInsSize(const Coil dim) const = 0;
    virtual double* GetInsSize() = 0;

    /// @brief calculate the area for conductor/strip/shell and its total
    ///        area and insulation area
    virtual double GetArea() const = 0;
    virtual double GetTotalArea() const = 0;
    virtual double GetInsArea() const = 0;

    /// @brief return the type of geometry of abstract class
    virtual Geometry GetType() const = 0;

    /// @brief convert enumeration to string
    std::string GetTypeName();

};

#endif
