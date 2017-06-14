/**
 *  @file   XQuenchTransient.hpp
 *  @author Y.Yang (Kyushu University)
 *  @date   7 Oct 2016
 **/

#ifndef XQuenchTransient_HH

#include "XThermalSolver.hpp"
#include "XTransientLoop.hpp"

class XQuenchTransient : public XTransientLoop
{
  public:
    /// constructor
    XQuenchTransient();

    /// deconstructor
    ~XQuenchTransient();

    /// @brief setup the hot spot and heat generation
    void SetHotSpot(const int z, const int phi, const int r, const double q=5000*100.);

    /// @brief setup dump resistor resistance
    void SetDumpResistor(const double R);

    /// @brief returns the dump resistor resistance
    double GetDumpResistor() const { return fDumpRes; }

    /// @brief setup the coil inductance
    void SetInductance(const double L);

    /// @brief returns the coil inductance
    double GetInductance() const { return fInduct; }

    /// @brief setup operating current
    void SetCurrent(const double I);

    /// @brief return the operating current
    double GetCurrent() const { return fCurr; }

    /// @brief  setup allowable voltage
    /// @detail calculate the dump resistance from allowable voltage.
    ///         to use this function, the current must be set before this function.
    void SetVoltage(const double V);

    /// @brief returns the decay time
    inline double GetDecaytime() { return fInduct/fDumpRes; }

    /// @brief setup threshold voltage
    void SetThreshold(const double Vth);

    /// @brief returns the threshold voltage
    double GetThreshold() const { return fVth; }

    /// @brief setup detection time
    void SetDetectTime(const double time);

    /// @brief returns detection time
    double GetDetectTime() const { return fDetTime; }

    /// @brief setup the diode turn on voltage
    void SetDiode(const double V);

    /// @brief return the diode turn on voltage
    double GetDiode() const { return fDiode; }

    /// @brief return the diode voltage from the given current
    double GetDiodeVoltage(const double I) const;

    /// @brief calculate the current decay
    double CalCurrentDecay(const double preI, const double res, const double dt);

    /// @brief calculate the magnetic field decay
    void CalFieldDecay(XThermalSolver* solver);

    /// @brief count the total conductor number
    int GetTotalConductor(XThermalSolver* solver);

    /// @brief count quenched conductor
    int GetQuenchConductor(XThermalSolver* solver, QuenchStatus qch=kTransition);

    /// @brief set beginning of the run
    virtual void Begin();

    /// @brief run transient loop
    virtual void Run();

    /// @brief end of transient loop
    virtual void End();


  protected:
    double fDumpRes;
    double fInduct;
    double fPreI;
    double fCurr;
    double fVth;
    double fDetTime;
    double fDiode;
    int    fHotZ;
    int    fHotPhi;
    int    fHotR;
    double fHotSpotHeat;
    double* fMagConnect;
    double* fShellConnect;
};

#endif
