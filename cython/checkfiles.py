#!/usr/bin/env python

import matplotlib.pyplot as plt
import sys
import XRootFiles

if __name__=="__main__":
    files = XRootFiles.XRoot3dPlots(sys.argv[1])
    files.SetMagnet("CS1")
    print "TIME: %.2f sec" %files.GetTime()
    files.Plot()
    plt.show()

