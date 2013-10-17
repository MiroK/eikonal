import sys
import numpy as np
from math import log

if __name__ == "__main__":
  file_name = sys.argv[1]
  data = np.loadtxt(file_name)

  print "\t%s\t|\t%s\t|\t%s\t"  % ("L^{1}", "L^{2}", "C^{oo}")

  for i in range(1, len(data)):
    h = log(data[i, 0])/log(data[i-1, 0])
    for j in range(1, 4):
      e = log(data[i, j])/log(data[i-1, j])
      print "\t%.2g\t" % (e/h),
    print 


