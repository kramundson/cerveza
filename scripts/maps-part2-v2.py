#! /usr/bin/env python

import os, sys, math 
from optparse import OptionParser

#Comai Lab, Ucdavis Genome Center
#Meric Lieberman, 2012
# This work is the property of UC Davis Genome Center - Comai Lab

# Use at your own risk. 
# We cannot provide support.
# All information obtained/inferred with this script is without any 
# implied warranty of fitness for any purpose or use whatsoever. 
#------------------------------------------------------------------------------
#
#MAPS - Part 2, maps-part2.py
#
#
#This program uses the list of potentially interesting positions generated by part1 and filters each of them based on a second set of more or less stringent criteria. Typically, users will run this second step multiple times to optimize the criteria, in order to obtain a final output that best represents their input data.
#
#INPUT: This program takes the list of potentially interesting position from the MAPS part 1 program as input.
#
#OUTPUT: 
#The user specifies an output file name for the output mutation / gentoyping file. 
#There is an non-assay-[input file name] file generated that is a table of the assayed positions that are lost do to cutoff parameters.
#
#PARAMETERS, Default value in []: 
#1. REQUIRED: 
#i. -f or --file, The input file is the output of the MAPS part 1 program [required] 
#ii. -o or --out, The output mutation / gentotyping file [required] 
#2. OPTIONAL: 
#i. --minCov or -v, Minimum total position coverage to be considered as a valid position. Should be the same as MAPS part 1. [1]
#
#The next four parameters analyze the basecalls at a particular position and for a particular library.
#
#In genotyping mode: each library is dealt with independently. 
#ii. --hetMinCov or -d, minimum coverage of each of the two observed major basecalls in order to be called heterozygous for that library. [5] 
#iii. --hetMinPer or -p, minimum percentage of each of the two observed major basecalls in order to be called heterozygous for that library. [20.00] 
#iv. --hetMinLibCov or -D, Minimum total coverage for consideration of heterozygous. This parameter is only used in genotyping mode. [5] 
#v. --homMinCov or -s, Minimum coverage in order for a position to be called homozygous for that library. [2]
#
#In mutation detection mode: a potential mutant allele is determined as present in only one library and not in any of the others (which all carry the WT allele). Whether or not it qualifies as a potential mutation depends on the following criteria: 
#ii. --hetMinCov or -d, minimum coverage of mutant allele to be called heterozygous for that library [5] 
#iii. --hetMinPer or -p, minimum percentage of the mutant and WT alleles in order to be called heterozygous for that library. [20.00] 
#v. --homMinCov or -s, Minimum coverage in order for a position to be called homozygous for that library. [2]
#
#In both modes again: 
#vi. -l or --MinLibs, Minimum number of libraries covered at least once to be considered a valid position. Should be the same as MAPS part 1. [3] 
#vii. --mode or -m, Output Mode: m == Mutation Detection, g == Genotyping. Should be the same as MAPS part 1. [m]






usage = "USAGE: %prog -f results-from-genotype.txt -o outputfile.txt [OPTIONS] [-h or --help]"
parser = OptionParser(usage=usage)
parser.add_option("-f", "--file", dest="f", default="sys.argv[1]", help="Input mpileup file.")
parser.add_option("-o", "--out", dest="o", default="sys.argv[2]", help="Output mpileup destination.")
parser.add_option("--minCov", "-c", dest="mincov", type = "int", default=1, help="Minimum position coverage")
parser.add_option("--hetMinCov", "-d", dest="hetMinCov", type = "int", default=5, help="Minimum coverage for consideration of het.")
parser.add_option("--hetMinPer", "-p", dest="hetMinPer", type = "float", default=20.00, help="Minimum percentagefor consideration of het.")
parser.add_option("--hetMinLibCov", "-D", dest="hetLibMinCov", type = "int", default=5, help="Minimum coverage for consideration of het.")
parser.add_option("--homMinCov", "-s", dest="homMinCov", type = "int", default=2, help="Minimum coverage for consideration of hom.")
parser.add_option("-l", "--MinLibs", dest="minlibs", type = "int", default=3, help="Part-1:Minimum number of libraries covered to be considered a valid position.")
parser.add_option("--mode", "-m", dest="mode",  type = "str", default='m', help="Output Mode:  m==Mutation Detection,  g==Genotyping")
(opt, args) = parser.parse_args()

try:
   #1. input mpileup file
   f = open(opt.f)
   #2. output file
   o = open(opt.o,'w')
   #3 Non-assay file, for libraries that get removed due to cutoff paramters
   # na = open("non-assay-"+opt.f,'w')
   na = open(os.path.join(os.path.dirname(opt.f), "non-assay-{}".format(os.path.basename(opt.f))))
except IOError:
   f = open(eval(opt.f))
   o = open(eval(opt.o),'w')
   # na = open("non-assay-"+eval(opt.f),'w')
   na = open(os.path.join(os.path.dirname(opt.f), "non-assay-{}".format(os.path.basename(opt.f))), 'w')
except:
   parser.error("Please check your command line paramters with -h or --help")

mincov = opt.mincov
minhet = opt.hetMinCov
minhetper = opt.hetMinPer
#only geno
minhetlibcov = opt.hetLibMinCov

minhom = opt.homMinCov
mode = opt.mode
minlibs = opt.minlibs

testtra = []
tc = 0
tcg = 0
tcb = 0
guide = ['A','C','G','T','+','*']
#split to length
def splitter(l, n):
    i = 0
    chunk = l[:n]
    while chunk:
        yield chunk
        i += n
        chunk = l[i:i+n]

# get and setup the header depending on the mode
head = f.readline()
header = head.split('\t')
newh = header[:4]+ ["#Libs"]

#liblist = map(lambda x: x.split('-')[-1],header[16::2])
liblist = map(lambda x: ''.join(x.split('-')[1:]),header[16::2])
if mode == 'g':
   o.write('\t'.join(newh+liblist)+'\n')
if mode == 'm':
   o.write('\t'.join(newh[:4]+["WT",'MA','Lib','Ho/He','WTCov','MACov', 'Type', 'LCov', '#libs', "InsertType"])+'\n')
#non assay table header
na.write('\t'.join(liblist)+'\n')

#for data line in the output file from part 1
for l in f:
   data = []
   bad = 0
   xline = l[:-1].split('\t')
   chrom = xline[0]
   pos = xline[1]
   ref = xline[2]
   totcova = xline[3]
   totcovt = int(xline[3])
   let = ['A','T','C','G','*','+']
   covs = xline[4:16][1::2]
   all = zip(covs, let)
   all.sort(lambda x, y: cmp(int(y[0]), int(x[0])))
   common = all[0][1]
   checker = []
   libs = list(splitter(xline[16:],2))
   #if mutation mode, figure out the most common base BY LIB COUNT
   # not by cov
   if mode == 'm':
      libs2 = [y[0] for y in libs]
      for t in libs2:
         if '-' in t:
            t1, t2 = t.split('-')
            checker.append(t1.split('_')[0]+'-'+t2.split('_')[0])
         else:
            checker.append(t.split('_')[0])
   
      types = list(set(checker))
      counts = [[checker.count(x), x] for x in types]
      counts.sort(lambda x, y: cmp( y[0], x[0] ))
      if '-' in counts[0][1] or len(counts) == 1:
         bad = 1
      else:
         common = counts[0][1]
         if common not in counts[1][1] and '-' in counts[1][1]:
            bad = 1
         else:
            common2 = counts[1][1]
   #if did not error during finding most common step
   if bad != 1:
      results = []     
      for lib in libs:
         #if mutation mode
         if mode == 'm':
            #if no lib data
            if lib[0] =='.':
               results.append('.')
            #elif a hom position, i.e. just a single base
            #if does not make min cov requirement, lib pos is removed from pos line
            elif '-' not in lib[0]:
               base, per = lib[0].split('_')
               if base != common:
                  if int(lib[1]) >= minhom:
                     results.append(base)
                  else:
                     results.append('.')
                     totcovt -= int(lib[1])
               else:
                  results.append(base)
            #else is a het lib position, i-e, two bases, such as C-T or A-G
            #if does not make min cov or min per requirement, lib pos is removed from pos line
            else:
               s1, s2 = lib[0].split('-')
               s1base, s1per = s1.split('_')
               s2base, s2per = s2.split('_')
               s1per = float(s1per)
               s2per = float(s2per)
               totcovl = int(lib[1])
               s1cov = round(float(s1per)/100.00*totcovl)
               s2cov = round(float(s2per)/100.00*totcovl)
               if s1base == common:
                  if s2cov >= minhet and s2per >= minhetper:
                     results.append(s1base+'-'+s2base)
                  else:
                     results.append('.')
                     totcovt -= int(lib[1])
               elif s2base == common:
                  if s1cov >= minhet and s1per >= minhetper:
                     results.append(s2base+'-'+s1base)
                  else:
                     results.append('.')
                     totcovt -= int(lib[1])   
               else:
                     results.append('.')
                     totcovt -= int(lib[1])
         #if gentoyping
         elif mode == 'g':
            #if no data for this lib in this positon
            if lib[0] =='.':
               results.append('.')
            #if a hom lib pos, i.e. single base
            #if does not make min cov requirement, lib pos is removed from pos line
            elif '-' not in lib[0]:
               base, per = lib[0].split('_')
               if int(lib[1]) >= minhom:
                  results.append(base)
               else:
                  results.append('.')
                  totcovt -= int(lib[1])
            #if a het lib pos, i.e 50/50 C/T
            #if does not make min cov or min per requirement, lib pos is removed from pos line
            else:
               s1, s2 = lib[0].split('-')
               s1base, s1per = s1.split('_')
               s2base, s2per = s2.split('_')
               totcovl = int(lib[1])
               s1cov = round(float(s1per)/100.00*totcovl)
               s2cov = round(float(s2per)/100.00*totcovl)
               if int(lib[1]) >= minhetlibcov:
                  if '+' in s1base or '*' in s1base:
                     sub = [s2base, s1base]
                  elif '+' in s2base or '*' in s2base:
                     sub = [s1base, s2base]
                  else:
                     sub = [s1base, s2base]
                     sub.sort(key = lambda x: guide.index(x))
                  if s1cov > s2cov:
                     if s2cov >= minhet and s2per >= minhetper:
                        results.append(sub[0]+'-'+sub[1])
                     else:
                        results.append('.')
                        totcovt -= int(lib[1])
                  elif s2cov > s1cov:
                     if s1cov >= minhet and s1per >= minhetper:
                        results.append(sub[0]+'-'+sub[1])
                     else:
                        results.append('.')
                        totcovt -= int(lib[1])
                  else:
                     if s1cov >= minhet and s2cov >= minhet:
                        results.append(sub[0]+'-'+sub[1])
                     else:
                        results.append('.')
                        totcovt -= int(lib[1])
               else:
                  results.append('.')
                  totcovt -= int(lib[1])
   #if has not become rejected yet, make sure meat min lib cutoff requirement                  
   if bad != 1:
      validlibnum = len(results) - results.count('.')   
      if validlibnum < minlibs:
         bad = 1
         na.write('\t'.join(xline[17::2])+'\n')
   
      test = list(set(results))
      try:
         test.remove('.')
      except:
         pass
   #if passed all tests so far, build output line
   #but first, make sure a C-G and C-G are treated as the same
   if bad != 1 and mode == 'g':
      newtest = []
      for con in test:
         if '-' in con:
            tcon = con.split('-')
            tcon.sort()
            tcon = '-'.join(tcon)
            newtest.append(tcon)
         else:
            newtest.append(con)
      newtest = list(set(newtest))
      if len(newtest) != len(test):
         test = newtest    
      
      
   if bad != 1:
      #genotyping output
      if len(test) != 1 and mode == 'g' and totcovt >= mincov:
         nline = [chrom, pos, ref, str(totcovt), validlibnum] + results
         nline = map(lambda x: str(x), nline)
         o.write('\t'.join(nline)+'\n')
      #mutation output
      #uses lib counts to determine wild type and mutant allele, errors/rejected if
      #not found confedently
      elif len(test) == 2 and mode == 'm' and totcovt >= mincov:
         line = [chrom, pos, ref, str(totcovt), common, common2] + results
         t1 = results.count(test[0])
         t2 = results.count(test[1])
         if t1 > t2 and t2 == 1:
            #if s1 is wild type
            wt = test[0]
            ma = test[1]
         elif t2>t1 and t1 == 1:
            wt = test[1]
            ma = test[0]
            #if t2 is wild type
         else:
            #these are those with multiple of the smaller
            bad = 1
         #WT can't be the het, MA has to share one base with the WT if a het
         if bad != 1:
            if '-' in wt:
               bad = 1
            if '-' in ma:
               if wt not in ma:
                  bad = 1

         if bad != 1:  
            mlib = liblist[results.index(ma)]
            mval = list(splitter(xline[16:],2))[results.index(ma)]
            mcov = mval[1]
         if bad != 1:         
            #take care of indel printing for orderly output format           
            if '-' not in ma:
               lt = "hom"
               ca1 = 0
               ca2 = int(mcov)
            else:
               lt = "het"
               if '-' in ma and '+' not in ma:
                  ma = ma.replace(wt,'')
                  ma = ma.replace('-','')
               if '+' in ma:
                  ma = ma.split('-')
                  ma.remove(wt)
                  ma = ma[0]
   
               s1, s2 = mval[0].split('-')
               s1base, s1per = s1.split('_')
               s2base, s2per = s2.split('_')
               ltotcov = int(mval[1])
               s1cov = int(round(float(s1per)/100.00*ltotcov))
               s2cov = int(round(float(s2per)/100.00*ltotcov))
               if s1base == wt and s2base == ma:
                  ca1 = s1cov
                  ca2 = s2cov
               elif s2base == wt and s1base == ma:
                  ca1 = s2cov
                  ca2 = s1cov
               else:
                  raise Exception("Check error base value, error in wt/ma selection for transition type field")
            line = [chrom, pos, ref, str(totcovt), wt, ma, mlib, lt, ca1, ca2]
            insertcount = '.'
            if '+' in ma:
               nma = '+'
               insertcount = ma
            else:
               nma = ma
            if '+' in wt:
               insertcount = wt
               wt = '+'
            if insertcount != '.':
               insertcount = insertcount.split('.+')[-1].replace('A','').replace('T','').replace('C','').replace('G','')
   
            
            #put result line together, and if meets all requirements(hasn't errored yet) then output
            line += [wt+nma, int(ca1)+int(ca2), validlibnum, insertcount]
            line = map(lambda x: str(x), line)
         if bad != 1:
         	o.write('\t'.join(line)+'\n')

f.close()
o.close()          
na.close()










