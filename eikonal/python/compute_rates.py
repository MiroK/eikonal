from math import log
from numpy import loadtxt
from os import listdir
from os.path import isfile


if __name__ == "__main__":
  
  for file_name in [ f for f in listdir(".")\
                         if isfile(f) and f.endswith(".dat") ]:
    # get all data
    data = loadtxt(file_name)

    rates_file_name = ".".join([file_name[:file_name.find(".")], "rates"])
    print rates_file_name
    with open(rates_file_name, "w") as rates_file:
      for i in range(1, len(data)):
        h = log(data[i, 0])/log(data[i-1, 0])
        print "%.2E" % data[i, 0],
        rates_file.write("%.2E" % data[i, 0])
        for j in range(1, 4):
          try:
            e = log(data[i, j])/log(data[i-1, j])
          except ValueError:
            e = 0
          
          print "%.2g\t" % (e/h),
          rates_file.write("%.2g\t" % (e/h))
        print
        rates_file.write("\n")
