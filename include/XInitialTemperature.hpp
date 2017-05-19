/**
 *  @file   XIntialTemperature.hpp
 *  @author Y.Yang (Kyushu University)
 *  @date   1 Dec 2016
 **/

#ifndef XInitialTemperature_HH
#define XInitialTemperature_HH

/// class to contain the initial temperature profile
//
class XInitTempContainer
{
  public:
    /// @brief constructor
    XInitTempContainer() 
      : fNode(0), fIdz(0), fIdp(0), fIdr(0), fTemp(4.5) {}

    /// @brief deconstructor
    ~XInitTempContainer() {}

    /// @brief setup node id
    void SetNodeId(const int id) { fNode = id; }

    /// @brief return the node id
    int GetNodeId() const { return fNode; }

    /// @brief setup the direction id number
    void SetId(const int idz, const int idp, const int idr) {
      fIdz = idz;
      fIdp = idp;
      fIdr = idr;
    }

    /// @brief return the direction id number by the given direciton enum
    int GetId(const Coil dim) const {
      switch (dim) {
        case iZ:   return fIdz; break;
        case iPhi: return fIdp; break;
        case iR:   return fIdr; break;
        default: break;
      }
    }

    int* GetId() const {
      int* id = new int[3];
      id[0] = fIdz;  id[1] = fIdp;  id[2] = fIdr;
      return id;
    }

    /// @brief setup the inital temperature
    void SetTemperature(const double T) { fTemp = T; }

    /// @brief return the initial temperautre
    double GetTemperature() const { return fTemp; }

  private:
    int    fNode;
    int    fIdz;
    int    fIdp;
    int    fIdr;
    double fTemp;
};

#include <vector>

/// class to load the intial temperature profile
//
class XInitialTemperature
{
  public:
    /// @brief constructor
    XInitialTemperature();

    /// @brief deconstructor
    ~XInitialTemperature();

    /// @brief setup factor
    void SetFactor(const double factor);

    /// @brief return the factor
    double GetFactor() const { return fFactor; }

    /// @brief load temperature profile
    void Load(const char* filename);

    /// @brief return the temperature container vector
    std::vector<XInitTempContainer*> GetContainer() { return fTemp; }


  private:
    double fIntialTemp;
    double fFactor;
    std::vector<XInitTempContainer*> fTemp;
};

#endif
