/**
 *  @file   XTransientLoop.hpp
 *  @author Y.Yang (Kyushu University)
 *  @date   1 Sep 2016
 **/

#ifndef XTransientLoop_HH
#define XTransientLoop_HH

#include "XProcessManager.hpp"

class XThermalSolver;

/// base time loop class -> QuenchLoop, ThermalLoop
class XTransientLoop
{
  public:
    /// @brief constructor
    XTransientLoop();

    /// @brief deconstructor
    ~XTransientLoop();

    /// @brief setup time mesh
    void SetTime(const double t0, const double tf, const double dt);

    /// @brief set the count step
    /// @detail the data will be displayed on the screen every <step> loop
    void SetDisplayStep(const int step);

    /// @brief return the solver
    XThermalSolver* GetSolver() { return fSolver; }

    /// @brief begin of run
    virtual void Begin();

    /// @brief setup process manager
    void SetProcess(Quench::XProcessManager* hand);

    /// @brief run
    virtual void Run();

    /// @brief end of run
    virtual void End();


  protected:
    XThermalSolver* fSolver;
    double fTime0;
    double fTimef;
    double fdt;
    int    fDisplay;
};

#endif
