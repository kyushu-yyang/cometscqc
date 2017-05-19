#!/usr/bin/env python

import sys
import XPostLoading as pt
import XPostOutput
import matplotlib.pyplot as plt

plt.style.use('classic')

def usage():
    print "usage: "
    print " ./demo.py <option> [arguments]"


def run(argv):
    geo = False
    dat = False
    info = pt.kTemperature
    dim = pt.kZ

    for i in range(len(argv)):
        if argv[i]=="-h" or argv[i]=="--help":
            usage()
        elif argv[i]=="-g" or argv[i]=="--geo":
            geofile = argv[i+1]
            geo = True
        elif argv[i]=="-d" or argv[i]=="--data":
            datfile = argv[i+1]
            dat = True
        elif argv[i]=="-t" or argv[i]=="--temp":
            info = pt.kTemperature
        elif argv[i]=="-r" or argv[i]=="--rrr":
            info = pt.kRRR
        else:
            raise

    plot = pt.XPost2dPlot(geofile, datfile)
    plot.SetMatInfo(info)
    plot.SetDirection(dim)
    plot.Draw()



if __name__=="__main__":
    argv = sys.argv[1:]

    out = XPostOutput.XPostOutput(sys.argv[1], sys.argv[2])
    out.Print()

    plot = pt.XPost2dPlot(sys.argv[1], sys.argv[2])
    #plot.SetMatInfo(pt.kTemperature)
    #plot.SetMatInfo(pt.kDose)
    #plot.SetMatInfo(pt.kRRR)
    plot.SetMatInfo(pt.kField)
    plot.SetDirection(pt.kZ)
    plot.SetPhi(3)
    plot.Draw()
    plt.show()
