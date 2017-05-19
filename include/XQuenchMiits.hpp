/**
 *  @file   XQuenchMiits.hpp
 *  @author Y.Yang (Kyushu University)
 *  @date   13 Sep 2016
 **/

#ifndef XQuenchMiits_HH
#define XQuenchMiits_HH

#include <map>
#include "XCoilHandle.hpp"

class XQuenchMiits
{
  public:
    /// @brief constructor
    XQuenchMiits();

    /// @brief deconstructor
    ~XQuenchMiits();

    /// @brief setup dump resistor resistance
    void SetDumpResistor(const double R);

    /// @brief returns the resistance of dump resistor
    double GetDumpResistor() const { return fResist; }

    /// @brief setup magnet total inductance
    void SetInductance(const double L);

    /// @brief returns the total inductance
    double GetInductance() const { return fIndct; }

    /// @brief  returns the time constance
    /// @detail time constance: \f$ \tau = \frac{L}{R} f$
    double GetTimeConstance() const;

    /// @brief setup the initial current
    void SetCurrent(const double I0);

    /// @brief setup material RRR
    void SetRRR(const double RRR);

    /// @brief setup magnetic field
    void SetField(const double B);

    /// @brief setup initial temperature
    void SetInitialTemp(const double T0);

    /// @brief set final temperature
    void SetFinalTemp(const double Tf);

    /// @brief setup temperature range
    void SetTempRange(const double T0, const double Tf);

    /// @brief setup coil handler
    void SetCoilHandler(Quench::XCoilHandle* coil) { fCoil = coil; }

    /// @brief  calculate miits from current
    /// @detail calculate miits by integrating the current \f$ \int_0^{\infty} I(t)^2 dt \f$,
    ///         assuming the current is decaied following the dump resistor and inductance,
    ///         \f$ Miits = \frac{L}{2R} I_0^2 \f$
    double GetMiits() const;

    /// @brief  calculate miits from current integration
    /// @param  nt total time step number
    /// @detail miits: \f$ \int_0^{\infty} J_0 exp(-2\frac{R}{L}t) dt  f$
    std::map<double,double> GetMiitsDecay(const double t0, const double tf, const double T);

    /// @brief  calculate miits from given temperature
    /// @detail miits at T: \f$ U(T_0) = \int_0^{T_0} \gamma \frac{C(T)}{\rho(T)}dT \f$
    double CalMiits2(const double Temp);

    /// @brief calculate miits from material property
    std::map<double, double> Eval(const double T0, const double Tf);

    /// @brief interpolate and calculate the maximum temperature from the given Miits
    double GetMaxTemperature(const double miits);

    /// @brief  return the minimal propagation zone
    /// @detail MPZ: \f$ l = \sqrt{\frac{2k(T_c - T_0)}{J_c^2 \rho}} \f$
    double GetMPZ() const;

    /// @brief  return the quench velocity approached from exact solution
    /// @detail Quench velocity: \f$ v_{ad} = \frac{J_{op}}{\gamma C} \sqrt{\frac{L_0 T_{cs}}{T_{cs} - T_0}} \f$
    double GetVelocity(const double Iop);

    /// @brief return the resistance with conductor self resistance
    double GetQuenchRes(const double t);
    
    /// @brief interpolation
    double interpolate(std::map<double,double> hc, const double x) const;

    /// @brief get the collection of time vs temperature
    std::map<double, double> GetTimeTemp(const double T0, const double Tf);

    /// @brief get the collection of time vs current
    std::map<double, double> GetTimeCurrent(const double T0, const double Tf);

    /// @brief get the collection of time vs resistance
    std::map<double, double> GetTimeResist(const double T0, const double Tf);

    /// @brief get the collection of time vs voltage
    std::map<double, double> GetTimeVolt(const double T0, const double Tf);
    

  protected:
    /// @brief calculate the average capacity for conductor
    double GetAvgCapacity(const double T, const double RRR, const double B);

    /// @brief calculate the average resistivity for conductor
    double GetAvgResistance(double T, double RRR, double B, double l=1.) const;

    /// @brief interpolate to find miits according to the miits
    double findtime(std::map<double,double> time, double miits);


  private:
    double fResist;
    double fIndct;
    double fCurr;
    double fField;
    double fRRR;
    double fTemp0;
    double fTempf;
    Quench::XCoilHandle* fCoil;
};

#endif
