#!/usr/bin/python
# depth_summarizer.py
# Kirk Amundson
# 1 May 2019

import argparse
import re
import statistics

def get_options():
    
    """
    Parse command line options
    """
    
    parser = argparse.ArgumentParser(description="group per-base depth (BED3 format) by intervals specified in a BED3 file")
    parser.add_argument("-w", "--windows", dest = "window_file", help = "Nonoverlapping windows")
    parser.add_argument("-b", "--bed", dest = "bed_file", help = "Per-position coverage file")
    parser.add_argument("-o", "--out", dest = "output_file", help = "Output file")
    parser.add_argument("-g", "--global-out", dest = "global_output_file", help = "Output file for global mean and median")
    args = parser.parse_args()
    return args

def parse_depth(bedfile, global_out):
    
    """
    Go through BED3 file.
    For each position, look up associated gene
    
    Store position-specific depth for each gene
    """
    
    bd = {}
    alldp = []
    
    with open(bedfile, 'r') as f:
        
        for line in f:
            if line[0] == "#":
                continue
            
            l = line.rstrip('\n').split('\t')
            chrom = l[0]
            
            if chrom not in bd:
                bd[chrom] = {}
                
            pos = int(l[1])
            dp  = int(l[2])
                
            if pos not in bd[chrom]:
                bd[chrom][pos] = dp
                alldp.append(int(dp))
                
            else:
                sys.exit("duplicate positions")
                
        avg=statistics.mean(alldp)
        median=statistics.median(alldp)
        o=open(global_out, 'w')
        out="Mean: {}\nMedian: {}\n".format(str(avg), str(median))
        o.write(out)
        o.close()
                
    return bd

def parse_windows(bd, win, out):
    
    """
    Go through BED3 file
    For each entry in BED3, record positions associated with that entry
    """
    
    bindict = {}
    out = open(out, 'w')
    
    with open(win, 'r') as f:
        for line in f:
            
            if line[0] == "#":
                continue
                
            l = line.rstrip().split('\t')
            chrom = re.sub("ST4.03ch", "chr", l[0])
            start = int(l[1])
            end   = int(l[2])
            # start = int(l[3])-1
            # end   = int(l[4])
            # gene  = l[8].split(";")[0].split("=")[1]
            bin = "{}_{}_{}".format(chrom, start, end)
            bindict[bin] = []
            
            if chrom in bd:
                
                for i in range(start, end):
                    if i in bd[chrom]:
                        bindict[bin].append(bd[chrom][i])
                        
                if bindict[bin] != []:
                    avg = statistics.mean(bindict[bin])
                    med = statistics.median(bindict[bin])
                    
                else:
                    avg = "NA"
                    med = "NA"
                
                oline = "{}\t{}\t{}\n".format(line.rstrip('\n'), avg, med)
                out.write(oline)
                
    out.close()
    return bindict

def main():
    args = get_options()
    parsed_bed = parse_depth(args.bed_file, args.global_output_file)
    bindict   = parse_windows(parsed_bed, args.window_file, args.output_file)

if __name__ == "__main__":
    main()