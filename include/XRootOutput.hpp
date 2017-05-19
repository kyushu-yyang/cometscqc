/**
 *  @file   XRootOutput.hpp
 *  @author Y.Yang (Kyushu University)
 *  @date   14 Oct 2016
 **/

#ifndef XRootOutput_HH
#define XRootOutput_HH

#include <string>
#include <map>
#include "XProcessManager.hpp"

class TFile;
class TTree;

/// class to dump the data into the root tree file
//
class XRootOutput
{
  public:
    /// constructor
    XRootOutput();
    XRootOutput(const char* filename);

    /// deconstructor
    ~XRootOutput();

    /// @brief setup filename
    void SetFilename(const char* filename);

    /// @brief write the time dependence info into header file
    void SetHeader(const int i, const double t, const double I, const double R, const double V);

    /// @brief return the root file object
    TFile* GetFile() { return fFile; }

    /// @brief set root sub directory
    void SetSubDirectory(const char* dirname);

    /// @brief fill the data into root file
    void Fill(const char* name, Quench::XProcessManager* man);

    /// @brief closed the root file
    void Close();

    /// @brief check the given path exists or not 
    void SetPath(const char* path);


  private:
    TFile* fFile;
    TTree* fInfo;
    std::string fName;
    int fMshZ;
    int fMshP;
    int fMshR;
    std::map<const char*, TDirectory*> fDir;
};


#endif
