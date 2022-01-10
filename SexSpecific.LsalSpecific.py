##########################################################
### Import Necessary Modules

import argparse                        #provides options at the command line
import sys                             #take command line arguments and uses it in the script
import gzip                            #allows gzipped files to be read
import re                              #allows regular expressions to be used

##########################################################
### Command-line Arguments

parser = argparse.ArgumentParser(description="A script to identify sex-specific scaffolds (looking for high missing areas between individuals)")
parser.add_argument("-file", help = "The statistic file (output by VCF.stats.windowBased.v1.0.py script)", default=sys.stdin, required=True)
parser.add_argument("-ind1", help = "The list of individuals seperated by a comma (no spaces) to compare, default = NA", default="NA")
parser.add_argument("-ind2", help = "The list of individuals seperated by a comma (no spaces) to compare against, default = NA", default="NA")
parser.add_argument("-win", help = "What was the window size used to generate data? default=10000", default=10000)
parser.add_argument("-fold", help = "The fold difference to report between groups, default 2", default=2)
args = parser.parse_args()


class OpenFile():
    def __init__ (self, f, typ, fnum):
        """Opens a file (gzipped) accepted"""
        if re.search(".gz$", f):
            self.filename = gzip.open(f, 'rb')
        else:
            self.filename = open(f, 'r')
        if typ == "file":
            if int(fnum) == 1:
                sys.stderr.write("\n\tOpened stats file: {}\n\n".format(f))
                ReadFile(self.filename,fnum,f)

class ReadFile():
    def __init__ (self, f, fnum, fname):
        self.openFile = f
        self.keep1 = args.ind1.split(",")
        self.keep2 = args.ind2.split(",")
        self.header = "NA"
        self.total = {}
        self.count = {}
        self.totalOthers = {}
        self.countOthers = {}
        if re.search(".gz$", fname):
            self.header = self.openFile.readline().decode('utf-8').rstrip('\n')
        else:
            self.header = self.openFile.readline().rstrip('\n')
        for line in self.openFile:
            try:
                line = line.decode('utf-8')
            except:
                pass        
            line = line.rstrip('\n')   
            self.ind, self.mis, self.hRef, self.Het, self.hAlt, self.aDepth, self.window, self.scaffold = line.split()
            if self.ind in self.keep1:
                self.value = int(self.mis)/float(int(self.mis) + int(self.hRef) + int(self.Het) + int(self.hAlt))
                if self.scaffold in self.total:
                    if self.window in self.total[self.scaffold]:
                        self.total[self.scaffold][self.window] += self.value
                        self.count[self.scaffold][self.window] += 1                            
                    else:
                        self.total[self.scaffold][self.window] = self.value
                        self.count[self.scaffold][self.window] = 1                        
                else:
                    self.total[self.scaffold] = {}
                    self.count[self.scaffold] = {}
                    self.total[self.scaffold][self.window] = self.value
                    self.count[self.scaffold][self.window] = 1
            elif self.ind in self.keep2:
                self.value = int(self.mis)/float(int(self.mis) + int(self.hRef) + int(self.Het) + int(self.hAlt))
                if self.scaffold in self.totalOthers:
                    if self.window in self.totalOthers[self.scaffold]:
                        self.totalOthers[self.scaffold][self.window] += self.value
                        self.countOthers[self.scaffold][self.window] += 1                            
                    else:
                        self.totalOthers[self.scaffold][self.window] = self.value
                        self.countOthers[self.scaffold][self.window] = 1                        
                else:
                    self.totalOthers[self.scaffold] = {}
                    self.countOthers[self.scaffold] = {}
                    self.totalOthers[self.scaffold][self.window] = self.value
                    self.countOthers[self.scaffold][self.window] = 1            
                        
        for self.scaffold in sorted(self.total):
            for self.window in sorted(self.total[self.scaffold], key = int):
                self.misFraction1 = float(self.total[self.scaffold][self.window])/int(self.count[self.scaffold][self.window])
                self.misFraction2 = float(self.totalOthers[self.scaffold][self.window])/int(self.countOthers[self.scaffold][self.window])
                if (float(self.misFraction1)+0.001)/(float(self.misFraction2)+0.001) >= float(args.fold) or (float(self.misFraction1)+0.001)/(float(self.misFraction2)+0.001) <= 1/float(args.fold):
                    print ("{}\t{}\t{}\t{}\t{}".format(self.scaffold, int(self.window) - (int(args.win) - 1), int(self.window), self.misFraction1, self.misFraction2))
                             

if __name__ == '__main__':
    open_file = OpenFile(args.file, "file", 1)
