##########################################################
### Import Necessary Modules

import argparse                        #provides options at the command line
import sys                             #take command line arguments and uses it in the script
import gzip                            #allows gzipped files to be read
import re                              #allows regular expressions to be used

##########################################################
### Command-line Arguments

parser = argparse.ArgumentParser(description="A script to partition individuals form a rolling window statistic script (VCF.stats.windowBased.v1.0.py)")
parser.add_argument("-file", help = "The statistic file", default=sys.stdin, required=True)
parser.add_argument("-ind", help = "The list of individuals seperated by a comma (no spaces), default = NA", default="NA")
parser.add_argument("-avg", help = "Find the average for individuals and only return that (doesn't output the header and format is for circos plot), default = no, option = yes", default="no")
parser.add_argument("-out", help = "If finding the average, what field? default=count, options=depth, het, mis", default="count")
parser.add_argument("-win", help = "What was the window size used to generate data? default=10000", default=10000)
parser.add_argument("-chr", help = "Keep the chromosomes in the following comma delimited list (no spaces) and convert to number. default=no", default="no")
args = parser.parse_args()

#########################################################
### Variables


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
        self.keep = args.ind.split(",")
        self.header = "NA"
        self.total = {}
        self.count = {}
        if re.search(".gz$", fname):
            self.header = self.openFile.readline().decode('utf-8').rstrip('\n')
        else:
            self.header = self.openFile.readline().rstrip('\n')
        if args.avg == "no":
            print("{}".format(self.header))
        for line in self.openFile:
            try:
                line = line.decode('utf-8')
            except:
                pass        
            line = line.rstrip('\n')   
            self.ind, self.mis, self.hRef, self.Het, self.hAlt, self.aDepth, self.window, self.scaffold = line.split()
            if self.ind in self.keep:
                if args.avg == "no":
                    print("{}".format(line))
                else:
                    self.value = "NA"
                    if args.out == "count":
                        self.value = int(self.mis) + int(self.hRef) + int(self.Het) + int(self.hAlt)
                    elif args.out == "depth":
                        try:
                            self.value = float(self.aDepth)
                        except:
                            self.value = "NA"
                    elif args.out == "het":
                        self.value = int(self.Het)/float(int(self.mis) + int(self.hRef) + int(self.Het) + int(self.hAlt))
                    elif args.out == "mis":
                        self.value = int(self.mis)/float(int(self.mis) + int(self.hRef) + int(self.Het) + int(self.hAlt))
                    if self.value != "NA":
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
        if args.avg == "yes":
            for self.scaffold in sorted(self.total):
                for self.window in sorted(self.total[self.scaffold], key = int):
                    if args.chr == "no":                  
                        print ("{}\t{}\t{}\t{}".format(self.scaffold, int(self.window) - (int(args.win) - 1), int(self.window), float(self.total[self.scaffold][self.window])/int(self.count[self.scaffold][self.window])))
                    else:
                        self.allChroms = args.chr.split(",")
                        for self.chrNum, self.chr in enumerate(self.allChroms):
                            if self.chr == self.scaffold:
                                print ("{}\t{}\t{}\t{}".format(int(self.chrNum) + 1, int(self.window) - (int(args.win) - 1), int(self.window), float(self.total[self.scaffold][self.window])/int(self.count[self.scaffold][self.window])))               

if __name__ == '__main__':
    open_file = OpenFile(args.file, "file", 1)
