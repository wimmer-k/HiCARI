#!/usr/bin/env python
import sys
import os
import fileinput

def main(argv):
    for r in range(448,483):
        if os.path.exists("/home/wimmer/programs/HiCARI/settings/runbyrun/all_run%04d_cal0411.dat" % (r)):
            os.system("cp set_base.dat set_run%04d.dat" % (r) )
            for line in fileinput.input('set_run%04d.dat' % r, inplace=True):
                if line.strip().startswith('HiCARI.Calibration.File:'):
                    line = 'HiCARI.Calibration.File:\t/home/wimmer/programs/HiCARI/settings/runbyrun/all_run%04d_cal0411.dat\n' % (r)
                sys.stdout.write(line)

if __name__=="__main__":
    sys.exit(main(sys.argv))
