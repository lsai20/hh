#!/usr/bin/python

# pylmm is a python-based linear mixed-model solver with applications to GWAS
# Copyright (C) 2015  Nicholas A. Furlotte (nick.furlotte@gmail.com)

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


# TODO
# write the output in binary (done)
# given an LD file, thresh, and whether shared or local:
#   treat LD as 0 or 1 depending on whether it crosses thresh
#   exclude SNPs that don't meet threshold
# given region, only use snps from that region out of the b/tfile
# edit the usage/output to match current fxn

# TODO add an option whether to write as either
# normal GRM (no LD) or binary with threshold'd LD
# (non-threshold'd LD shelved atm, may reimplement version in draft later)

import sys
import pdb

from optparse import OptionParser,OptionGroup
usage = """usage: %prog [options] --[tfile | bfile] plinkFileBase outfileBase ldFile [shared | local] ldThresh
"""

parser = OptionParser(usage=usage)

basicGroup = OptionGroup(parser, "Basic Options")
#advancedGroup = OptionGroup(parser, "Advanced Options")

#basicGroup.add_option("--pfile", dest="pfile",
#                  help="The base for a PLINK ped file")
basicGroup.add_option("--tfile", dest="tfile",
                  help="The base for a PLINK tped file")
basicGroup.add_option("--bfile", dest="bfile",
                  help="The base for a PLINK binary ped file")
basicGroup.add_option("--emmaSNP", dest="emmaFile", default=None,
                  help="For backwards compatibility with emma, we allow for \"EMMA\" file formats.  This is just a text file with individuals on the rows and snps on the columns.")
basicGroup.add_option("--emmaNumSNPs", dest="numSNPs", type="int", default=0,
         help="When providing the emmaSNP file you need to specify how many snps are in the file")

basicGroup.add_option("-e", "--efile", dest="saveEig", help="Save eigendecomposition to this file.")
basicGroup.add_option("-n", default=1000,dest="computeSize", type="int", help="The maximum number of SNPs to read into memory at once (default 1000).  This is important when there is a large number of SNPs, because memory could be an issue.")

basicGroup.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=False,
                  help="Print extra info")

parser.add_option_group(basicGroup)
#parser.add_option_group(advancedGroup)

(options, args) = parser.parse_args()
if len(args) != 4: 
   parser.print_help()
   sys.exit()

###   additional arguments    ###
outFile = args[0]
ldFile = args[1]
ldFlagStr = args[2] # shared or local
ldThresh = float(args[3])
makeLocal = None

if ldFlagStr == "shared":
  makeLocal = False
elif ldFlagStr == "local":
  makeLocal = True
else:
  print "Error: invalid LD flag: %s" % ldFlagStr
  print "LD flag should specify GRM type as either 'shared' or 'local' " 
  exit(1)



print "pylmmKinship_binary_thresh.py running with:"
print "\toutFileBase = %s" % outFile
print "\tldFile = %s" % ldFile
print "\tldFlag = %s" % ldFlagStr
print "\tldThresh = %f" % ldThresh


import os
import numpy as np
from scipy import linalg
from pylmm.lmm import calculateKinship
from pylmm import input

if not options.tfile and not options.bfile and not options.emmaFile: 
   parser.error("You must provide at least one PLINK input file base (--tfile or --bfile) or an emma formatted file (--emmaSNP).")

if options.verbose: sys.stderr.write("Reading PLINK input...\n")
if options.bfile: IN = input.plink(options.bfile,type='b')
elif options.tfile: IN = input.plink(options.tfile,type='t')
#elif options.pfile: IN = input.plink(options.pfile,type='p')
elif options.emmaFile: 
   if not options.numSNPs: parser.error("You must provide the number of SNPs when specifying an emma formatted file.")
   IN = input.plink(options.emmaFile,type='emma')
else: parser.error("You must provide at least one PLINK input file base (--tfile or --bfile) or an emma formatted file (--emmaSNP).")

n = len(IN.indivs)
m = options.computeSize
W = np.ones((n,m)) * np.nan

print "reading max LD file %s" % ldFile
f = open(ldFile, "r")
ldMaxes = []
for line in f:
    ldMaxes.append(float(line))
f.close()


IN.getSNPIterator()
# Annoying hack to get around the fact that it is expensive to determine the number of SNPs in an emma file
if options.emmaFile: IN.numSNPs = options.numSNPs
i = 0
K = None
while i < IN.numSNPs:
   j = 0
   while j < options.computeSize and i < IN.numSNPs:
      snp,id = IN.next()
      if snp.var() == 0:
        i += 1
        continue

### main changes here ###
      ## count snps with low max LD toward local
      if ldFlagStr == "local":
        if ldMaxes[i] < ldThresh:
          W[:,j] = snp
          j += 1
          i += 1
        else: # exclude snp
          i += 1
          continue

      ## count snps with high max LD toward shared 
      elif ldFlagStr == "shared":
        if ldMaxes[i] >= ldThresh:
          W[:,j] = snp
          j += 1
          i += 1
        else: # exclude snp
          i += 1
          continue

      else:
        print "invalid LD flag"

   if j < options.computeSize: W = W[:,range(0,j)] 

### main changes above ###


   if options.verbose: sys.stderr.write("Processing first %d SNPs\n" % i)
   sz = W.shape()
   print("test")
   sys.stderr.write("Size of W:  %s \n" % sz)


   if K == None: 
      try: 
         K = linalg.fblas.dgemm(alpha=1.,a=W.T,b=W.T,trans_a=True,trans_b=False) # calculateKinship(W) * j
      except AttributeError: K = np.dot(W,W.T) 
   else:
      try: 
         K_j = linalg.fblas.dgemm(alpha=1.,a=W.T,b=W.T,trans_a=True,trans_b=False) # calculateKinship(W) * j
      except AttributeError: K_j = np.dot(W,W.T)
      K = K + K_j



K = K / float(IN.numSNPs)


###    Saving binary files .grm.N.bin, .grm.bin     ###

if options.verbose: sys.stderr.write("Creating binary file %s.grm.N.bin (contains the number of SNPs used to calculate the GRM) \n" % outFile)

from array import array
import struct
import shutil   # for copying tfam to grm.id

outFile_N_bin=outFile+".grm.N.bin"
ints=array('d',[])

for j in range(0, n): # for pairs i >= j
    for i in range(j, n):
        ints.append(IN.numSNPs)

s = struct.pack('f'*len(ints), *ints)

f = open(outFile_N_bin,'wb')
f.write(s)
f.close()

if options.verbose: sys.stderr.write("Saving binary Kinship file to %s.grm.bin\n" % outFile)

outFile_bin = outFile + ".grm.bin"
grmArray = array('d',[])

for j in range(0, n): # for i,j where i >= j
    for i in range(j, n):
        grmArray.append(K[i,j])

s = struct.pack('f'*len(grmArray), *grmArray)
f = open(outFile_bin,'wb')
f.write(s)
f.close()


###  TODO save grm.id  ###
#if options.verbose: sys.stderr.write("Saving Kinship ID file to %s.grm.ID\n" % outFile)

#outFile_id = outFile + ".grm.ID"
#shutil.copyfile(src_tfam_TODO+".tfam", outFile_id) # if tfile, tfam
#shutil.copyfile(src_tfam_TODO+".fam", outFile_id) # if bfile, fam



if options.saveEig:
   if options.verbose: sys.stderr.write("Obtaining Eigendecomposition\n")
   Kva,Kve = linalg.eigh(K)
   if options.verbose: sys.stderr.write("Saving eigendecomposition to %s.[kva | kve]\n" % outFile)
   np.savetxt(outFile+".kva",Kva)
   np.savetxt(outFile+".kve",Kve)
      