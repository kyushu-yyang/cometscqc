/**
 *  @file   XThermalSolver.hpp
 *  @author Y.Yang (Kyushu University)
 *  @date   31 Aug 2016
 **/

#ifndef XThermalSolver_HH
#define XThermalSolver_HH

#include <vector>
#include <fstream>
#include "XProcessManager.hpp"

/*
enum Cooling
{
  kLeft,
  kRight,
  kSide
};
*/

/// class to solve the thermal transfer equation
//
class XThermalSolver
{
  public:
    /// @brief constructor
    XThermalSolver();

    /// @brief deconstructor
    ~XThermalSolver();

    /// @brief setup process manager
    void SetProcessHandle(Quench::XProcessManager* hand);

    /// @brief setup time interval
    void SetTimeInterval(const double dt);

    /// @brief setup accelerating factor
    void SetAccelerateFactor(double acc=1.);
    
    /// @brief initialization
    void Initial();

    /// @brief find minimum time step
    double FindTimeStep() const;

    /// @brief  solve the thermal problem
    /// @detail time step is calculated by the condition: \f$ r_{x} + r_{y} + r_{z} < \frac{1}{2} \f$
    void Solve(const double dt);

    /// @brief setup boundary
    void SetBoundary();

    /// @brief connect cell to cell
    void Connect(const int from_z, const int from_p, const int from_r,
                 const int to_z  , const int to_p  , const int to_r);
    void Connect(const int from_id, const int to_id);

    /// @brief setup cooling to given temperature
    /// @param r the quench code set the default cooling path is along r diection.
    ///          and default cooling path is at strip. if the r!=kStrip, then warning
    ///          will be output.
    /// @param opt option for selecting the way to cool down the magnet, kSide means
    ///            the both side cooling.
    /// @detail this function must be used after SetBoundary()
    void SetCoolingPath(const int r, const double T=4.5*K, const Cooling opt=kSide);

    /// @brief  setup cooling point on shell
    /// @detail this point is on the last layer of coil
    void SetLastCoolingPoint(const int z, const double T=4.5*K);

    /// @brief  setup cooling point on the inner most layer
    void SetFirstCoolingPoint(const int z, const double T=4.5*K);

    /// @brief return process handler
    Quench::XProcessManager* GetProcess() { return fProcess; }

    /// @brief print info
    void Print(const int z, const int phi, const int r);

    /// @brief add output material property
    void AddOutput(const int z, const int phi, const int r);

    /// @brief use the cylindrical connection at phi direction
    void UseCylinderConnect();

    /// @brief use the conductor connection for phi thermal conduction
    void UseConductorConnect();


  protected:
    /// @brief in the 3d loop
    void InTheLoop(const int i, const int j, const int k);

    /// @brief set phi cylinerical connection
    void SetCylinderPhi();

    /// @brief set phit conductor connection
    void SetConductorPhi();


  private:
    Quench::XProcessManager* fProcess;
    int fMshZ;
    int fMshP;
    int fMshR;
    double fAcce;
    double fdt;
    bool fCylinder;
    std::vector<int> fPrint;
};

#endif
