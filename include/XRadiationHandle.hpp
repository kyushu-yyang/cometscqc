/**
 *  @file   XRadiationHandle.hpp
 *  @author Y.Yang (Kyushu University)
 *  @date   30 Aug 2016
 **/

#ifndef XRadiationHandle_HH
#define XRadiationHandle_HH

#include <string>
#include <array>
#include "XCoilHandle.hpp"

enum Radiation
{
  kNeutron,
  kDPA,
  kDose
};

/// class to contain the info of radiation
//
class XRadiationContainer
{
  public:
    /// @brief deconstructor
    ~XRadiationContainer() {}

    /// @brief setup data's id
    void SetId(const int z, const int p, const int r);

    /// @brief return data id
    std::array<int,3> GetId() { return fId; }
    const int GetId(const Coil dim) { return fId.at(dim); }

    /// @brief setup position
    void SetPosition(const double z, const double p, const double r);

    /// @brief return positin
    std::array<double,3> GetPosition() { return fPos; }
    double GetPosition(const Coil dim) { return fPos.at(dim); }

    /// @brief setup neutron fluence (>1MeV)
    void SetNeutron(const double neu) { fNeutron = neu; }

    /// @brief return neutron fluence
    double GetNeutron() const { return fNeutron; }

    /// @brief setup dose ratio [W/kg]
    void SetDose(const double dose) { fDose = dose; }

    /// @brief return dose ratio
    double GetDose() const { return fDose; }

    /// @brief setup DPA
    void SetDPA(const double dpa) { fDpa = dpa; }

    /// @brief return DPA
    double GetDPA() const { return fDpa; }


  private:
    std::array<int, 3> fId;
    std::array<double,3> fPos;
    double fNeutron;
    double fDose;
    double fDpa;
};


/// class to read and handle the radiation data
//
class XRadiationHandle
{
  public:
    /// @brief constructor
    XRadiationHandle();
    XRadiationHandle(const std::string& filename);

    /// @brief deconstructor
    ~XRadiationHandle() {}

    /// @brief setup name of magnet
    void SetName(const std::string& name) { fName = name; }

    /// @brief return name of magnet
    std::string GetName() const { return fName; }

    /// @brief setup irradiation time
    void SetIrrTime(const double time);

    /// @brief check the id is over range or not
    bool IsOverRange(const int id) const;

    /// @brief find the local id from the given i, j, k
    const int Id(const int i, const int j, const int k);

    /// @brief load radiation file
    void Load(const std::string& filename);

    /// @brief return the mesh of radiation data
    const int GetMesh(const Coil dim);

    /// @brief return the radiation collect
    std::vector<XRadiationContainer*> GetCollection() { return fRC; }

    /// @brief get total entries of this collection
    size_t GetEntries() const { return fRC.size(); }

    /// @brief get this entry
    XRadiationContainer* GetEntry(const int entry) { return fRC.at(entry); }

    /// @brief calculate irradiation induced RRR
    double GetRRR(const Geometry geo, const double neu) const;


  private:
    /// @brief find the maximum of id
    int findmax(const Coil dim) const;


  private:
    std::string fName;
    double fIrrad;
    int fMshz;
    int fMshp;
    int fMshr;
    std::vector<XRadiationContainer*> fRC;
};

#endif
