/**
 *  @file   XPostInfoPlot.hpp
 *  @author Y.Yang (Kyushu University)
 *  @date   24 Aug 2016
 **/

#ifndef XPostInfoPlot_HH
#define XPostInfoPlot_HH

#include <vector>
#include "XProcessManager.hpp"

class TGraph;
class TGraph2D;

/// class to plot the information from quench container
//
class XPostInfoPlot
{
  public:
    /// @brief deconstructor
    ~XPostInfoPlot();

    /// @brief pass the XProcessManager class
    void SetProcessManager(Quench::XProcessManager* pro);

    /// @brief plot the position of each node
    void PlotNode(const bool save=false);

    /// @brief plot 2d magnetic field map
    void PlotField(const bool save=false);


  protected:
    /// @brief return the graph of node position
    /// @param phi id number along the phi direction
    TGraph* GetNodeGraph(const int phi);

    /// @brief return the 2d graph of magnetic field
    TGraph2D* GetFieldGraph();

  private:
    Quench::XProcessManager* fProcess;
};

#endif
