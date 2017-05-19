/**
 *  @file   XQuenchOutput.hpp
 *  @author Y.Yang (Kyushu University)
 *  @date   27 Aug 2016
 **/

#ifndef XQuenchOutput_HH
#define XQuenchOutput_HH

#include <fstream>
#include <string>
#include "XProcessManager.hpp"

enum FileOutput
{
  iROOT,
  iOfstream
};

namespace Quench
{ class XQuenchOutput; }

class TFile;

/// class to write the data information into a file
//
class Quench::XQuenchOutput
{
  public:
    /// @brief constructor
    XQuenchOutput(const std::string &filename, const FileOutput opt=iOfstream);

    /// @brief deconstructor
    ~XQuenchOutput();
    
    /// @brief write the data into file
    void Write(XProcessManager* man);

    /// @brief write geometry information
    void WriteGeometry(XProcessManager* man);

    /// @brief set the header stream
    void SetHeader(XProcessManager* man);

    /// @brief close the file
    void Close();


  protected:
    /// @brief initialization
    void init(const FileOutput opt);

  private:
    std::string fFilename;
    std::ofstream fNormfile;
    TFile* fRootfile;
};

#endif
