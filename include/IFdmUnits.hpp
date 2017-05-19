#ifndef IFdmUnits_HH
#define IFdmUnits_HH

#include <cmath>

static const double m  = 1.;
static const double m2 = m * m;
static const double m3 = m * m * m;

static const double cm  = 0.01 * m;
static const double cm2 = cm * cm;
static const double cm3 = cm * cm * cm;

static const double mm  = 0.001 * m;
static const double mm2 = mm * mm;
static const double mm3 = mm * mm * mm;

static const double kg = 1.;
static const double g  = 0.001 * kg;

static const double K  = 1.;
static const double mK = 0.001 * K;

static const double sec  = 1.;
static const double msec = 0.001*sec;
static const double hour = 3600.*sec;
static const double day  = 24.*hour;
static const double year = 365*day;

static const double Amp  = 1.;
static const double mAmp = 1e-3*Amp;
static const double kAmp = 1e+3*Amp;

//static const double V  = 1.;
//static const double mV = 1e-3*V;
//static const double kV = 1e+3*V;

static const double Ohm = 1.;
static const double mu0 = 4. * M_PI * 1e-7;

static const double Lwf = 2.45e-8;

#endif
