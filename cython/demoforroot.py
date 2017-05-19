#!/usr/bin/env python

import os
import sys
sys.path.append("/home/yeyang/Documents/software/CometQuenchCode-/cython")
import matplotlib.pyplot as plt
import XPostRootLoad as pt
import XPostOutput
import XRootFiles

def usage():
    print "usage:"
    print "      /path/demoforroot -m magnet -g geofile -d datafile"


def Draw(ax, geofile, datafile, magnet, Tmin, Tmax):
    plot = pt.XPost2dPlot(geofile, datafile, magnet)
    #plot.SetMatInfo(pt.kField)
    plot.SetMatInfo(pt.kTemperature)
    #plot.SetMatInfo(pt.kVoltage)
    plot.SetDirection(pt.kZ)
    plot.SetPhi(1)
    #plot.SetColorMap("nipy_spectral")
    plot.SetColorMap("jet")
    plot.SetRange(Tmin, Tmax)
    plot.DrawThis(ax)
    ax.set_title(magnet)


if __name__=="__main__":
    magnet   = ["CS0", "CS1", "MS1", "MS2"]
    geofile  = ["geo%s.dat" %magnet[0], "geo%s.dat" %magnet[1], "geo%s.dat" %magnet[2], "geo%s.dat" %magnet[3]]
    datafile = sys.argv[1]

    rfile = XRootFiles.XRootFiles(datafile)
    rfile.SetSubDirectory(magnet)
    Tmax = rfile.FindMaxTemp()
    Tmin = rfile.FindMinTemp()
    print " TIME: %.2f sec" %(rfile.GetTime())

    path = os.path.abspath(os.getcwd()+"/output")
    out = XRootFiles.XRootOutput(path)
    out.SetCoilPosition("CS1", 24, 1, 2)
    out.ConstructFigure()
    out.Plot("T")
    plt.show()

    fig, ax = plt.subplots(2, 2, figsize=(12,6))

    for i in range(len(magnet)):
        Draw(ax[i/2][i%2], geofile[i], datafile, magnet[i], Tmin, Tmax)

    plt.tight_layout()
    #plt.savefig("temp.pdf")
    plt.show()

    """
    plot = pt.XPost2dPlot(geofile, datafile, magnet)
    plot.SetMatInfo(pt.kTemperature)
    plot.SetDirection(pt.kZ)
    plot.SetPhi(1)
    plot.SetRange(4.5, maxi)
    plot.Draw()
    """
