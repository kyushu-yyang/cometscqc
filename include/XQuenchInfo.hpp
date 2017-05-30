/**
 *  @file   XQuenchContainer.hpp
 *  @author Y.Yang (Kyushu University)
 *  @date   1st. Aug. 2016
 **/

#ifndef XQuenchInfo_HH
#define XQuenchInfo_HH

#if __cplusplus <= 199711L
 #error This library needs C++11 at least
#endif
#include <array>
#include "XCoilBase.hpp"

/// @enum Quench status
enum QuenchStatus {
  /// superconducting status: operating
  kSuperconduct,
  /// transition status: start current sharing
  kTransition,
  /// normal: start quenching
  kNormal
};


namespace Quench
{ 
  class XMaterialInfo; 
  class XDimensionInfo;
}


/// class to store the material parameter
//
class Quench::XMaterialInfo
{
  public:
    /*! constructor */
    XMaterialInfo();

    /*! deconstructor */
    ~XMaterialInfo() {}

    /*! setup magnetic field */
    void SetField(const double &fld);

    /*! return the magnetic field */
    double GetField() const { return fField; }

    /// @brief setup pre-temperature
    void SetPreTemp(const double T) { fPreTemp = T; }

    /// @brief returns pre-temperature
    double GetPreTemp() const { return fPreTemp; }

    /*! setup temperature */
    void SetTemperature(const double temp);

    /*! return temperature */
    double GetTemperature() const { return fTemp; }

    /*! setup heat capacity */
    void SetCapacity(const double C);

    /*! return heat capacity */
    double GetCapacity() const { return fCapacity; }

    /// @brief setup the heat flux 
    void SetHeatFlux(const double Qx, const double Qy, const double Qz);

    /// @brief return heat generation
    /// @param dim enumeration of dimension: iZ/iPhi/iR
    double GetHeatFlux(const Coil dim) const { return fHeat.at(dim); }
    std::array<double,3> GetHeatFlux() const { return fHeat; }

    /*! @brief setup the heat generation */
    void SetHeat(const double gen) { fGen = gen; }

    /*! @brief return the heat generation */
    double GetHeat() const { return fGen; }

    /// @brief setup something user want to save it
    void SetStack(const double stack) { fStack = stack; }

    /// @brief returns the stack
    double GetStack() const { return fStack; }

    /// setup thermal conductivity
    void SetConductivity(const double kx, const double ky, const double kz);

    /// @brief return thermal conductivity 
    double GetConductivity(const Coil dim) const { return fk.at(dim); }
    std::array<double,3> GetConductivity() const { return fk; }

    /*! setup local resistance */
    void SetResistance(const double R) { fR = R; }

    /*! return local resistance */
    double GetResistance() const { return fR; }

    /// setup conductor voltage
    void SetVoltage(const double V) { fVolt = V; }

    /// return conductor local voltage
    double GetVoltage() const { return fVolt; }

    /*! @brief setup material density */
    void SetDensity(const double rho) { fRho = rho; }

    /*! @brief return material density */
    double GetDensity() const { return fRho; }

    /*! @brief setup residual resistance ratio */
    void SetRRR(const double RRR) { fRRR = RRR; }

    /*! @brief return residual resistance ratio */
    double GetRRR() const { return fRRR; }

    /*! @brief setup the superconducting status */
    void SetStatus(const QuenchStatus status) { fStatus = status; }

    /*! @brief return the superconducting status */
    QuenchStatus GetStatus() const { return fStatus; }

    /// @brief setup the static data of heat generation from radiation
    void SetDeposit(const double dose) { fDose = dose; }

    /// @brief return the energy deposit data
    double GetDeposit() { return fDose; }

    /// @brief setup the time step
    void SetTimeStep(const double step) { fStep = step; }

    /// @brief return the time step
    double GetTimeStep() const { return fStep; }

    /// @brief setup the quenched time
    void SetQuenchTime(const double t) { fQcht = t; }

    /// @brief return the quenched time
    double GetQuenchTime() const { return fQcht; }


  private:
    std::array<double, 3> fHeat;
    std::array<double, 3> fk;
    double  fField;
    double  fPreTemp;
    double  fTemp;
    double  fCapacity;
    double  fR;
    double  fRho;
    double  fRRR;
    QuenchStatus fStatus;
    double  fGen;
    double  fDose;
    double  fStep;
    double  fVolt;
    double  fStack;
    double  fQcht;
};


/// class discription:
/// class to handle the dimensional container
//
class Quench::XDimensionInfo
{
  public:
    /// @brief deconstructor
    ~XDimensionInfo() {}

    /// @brief setup cell id
    void SetId(const int i, const int j, const int k);

    /// @brief return cell id
    /// @param dim input the enumeration of dimension: iZ/iPhi/iR
    int GetId(const Coil dim) const { return fId.at(dim); }
    std::array<int,3> GetId() const { return fId; }

    /// @brief setup pre position
    void SetPrePosition(const double x, const double y, const double z);

    /// @brief returns the pre position
    double GetPrePosition(const Coil dim) const { return fPrePos.at(dim); }
    std::array<double,3> GetPrePosition() const { return fPrePos; }

    /// @brief setup cell position
    void SetPosition(const double x, const double y, const double z);

    /// @brief return cell position
    double GetPosition(const Coil dim) const { return fPos.at(dim); }
    double GetPosition2(const Coil dim) const { return (fPrePos.at(dim)+fPostPos.at(dim))/2.; }
    std::array<double,3> GetPosition() const { return fPos; }

    /// @brief setup post position
    void SetPostPosition(const double x, const double y, const double z);

    /// @brief return post position
    double GetPostPosition(const Coil dim) const { return fPostPos.at(dim); }
    std::array<double,3> GetPostPosition() const { return fPostPos; }

    /// @brief setup node id
    void SetNodeId(const int node);

    /// @brief return the node id
    int GetNodeId() const { return fNode; }

    /// @brief setup the geometry of this cell
    /// @param geo geometry enumeration: kConductor/kStrip/kShell
    void SetGeometry(const Geometry geo);

    /// @brief return the geometry of this cell
    Geometry GetGeometry() const { return fGeo; }

    /// @brief  set cell size
    /// @detail cell size is equal to the size of conductor or strip.
    void SetCellSize(const double lx, const double ly, const double lz); 

    /// @brief return cell size
    double GetCellSize(const Coil dim) const { return fCell.at(dim); }
    std::array<double,3> GetCellSize() const { return fCell; }

    /// @brief setup distance between two node
    void SetDistance(const double dx, const double dy, const double dz);

    /// @brief return distance between two node
    double GetDistance(const Coil dim) const { return fDistance.at(dim); }
    std::array<double,3> GetDistance() const { return fDistance; }


  private:
    std::array<   int, 3> fId;
    std::array<double, 3> fPos;
    std::array<double, 3> fPrePos;
    std::array<double, 3> fPostPos;
    std::array<double, 3> fCell;
    Geometry fGeo;
    double fNode;
    std::array<double, 3> fDistance;
};


#endif
