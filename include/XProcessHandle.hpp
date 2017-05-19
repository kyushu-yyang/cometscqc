/**
 *  @file   XProcessHandle.hpp
 *  @author Y.Yang (Kyushu University)
 *  @date   30 Aug 2016
 **/

#ifndef XProcessHandle_HH
#define XProcessHandle_HH

#include "XProcessManager.hpp"

namespace Quench
{ class XProcessHandle; }

/// class to handle the data container and update the data
/// @detail this class is for the quench calculation
//
class Quench::XProcessHandle : public Quench::XProcessManager
{
  public:
    /// @brief constructor
    XProcessHandle();

    /// @brief deconstructor
    virtual ~XProcessHandle();

    /// @brief setup current decay
    void SetCurrent(const double I) { fCurrent = I; }

    /// @brief setup field decay
    void SetField(const double B) { fField = B; }

    /// @brief setup current and field decay
    void SetDecayParameter(const double I, const double B);


  private:
    double fCurrent;
    double fField;
};

#endif
