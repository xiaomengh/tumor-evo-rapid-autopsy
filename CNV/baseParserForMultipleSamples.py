#! /usr/bin/env python
###############################################################################
# parses pileup base string and returns the counts for all possible alleles
# for each position
# reads input (mpileup output) from sys.stdin
###############################################################################

import os
import sys

class parseString(object):

    def __init__(self, ref, string):
        self.ref = ref.upper()
        self.string = string.upper()
        self.types = {'ref':0, 'Balt':0}
        self.process()

    def process(self):
        # remove end of read character
        self.string = self.string.replace('$','')
        while self.string != '':
            if self.string[0] == '^':
                # skip two characters when encountering '^' as it indicates
                # a read start mark and the read mapping quality
                self.string = self.string[2:]
            elif self.string[0] == '*':
                #self.types['*'] += 1
                # skip to next character
                self.string = self.string[1:]

            elif self.string[0] in ['.',',']:
                if (len(self.string)== 1) or (self.string[1] not in ['+','-']):
                    # a reference base
                    self.types['ref'] += 1
                    self.string = self.string[1:]
                elif self.string[1] == '+':
                    insertionLength = int(self.string[2])
                    # insertionSeq = self.string[3:3+ insertionLength]
                    # self.types['+'].append(insertionSeq)
                    self.string = self.string[3+insertionLength:]
                elif self.string[1] == '-':
                    deletionLength = int(self.string[2])
                    # deletionSeq = self.string[3:3+deletionLength]
                    # self.types['-'].append(deletionSeq)
                    self.string = self.string[3+deletionLength:]

            elif self.string[0] in ['A','a','T','t','C','c','G','g'] and\
                 ((len(self.string)==1) or (self.string[1] not in ['-','+'])):
                # one of the four bases
                self.types['Balt'] += 1
                self.string = self.string[1:]
            else:
                # unrecognized character
                # or a read that reports a substitition followed by an insertion/deletion
                #self.types['X'].append(self.string[0])
                self.string = self.string[1:]
        return
    def __repr__(self):
        types = self.types
        return '\t'.join([str(types['ref']), str(types['Balt'])])


def main():
    for line in sys.stdin:
        toks = line.strip('\n').split('\t')
        ref = toks[2].upper()
        if ref != "N":
            samples=[]
            for i in range(int(len(toks)/3)-1):
                cov = toks[(i+1)*3]
                samples.append(cov+'\t'+parseString(ref,toks[(i+1)*3+1]).__repr__()) 
            print('\t'.join([toks[0], toks[1],str(int(toks[1])+1),ref]+samples))

if __name__ == '__main__':
    main()
