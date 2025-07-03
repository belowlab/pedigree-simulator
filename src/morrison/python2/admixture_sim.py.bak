#!/usr/bin/python
#Modified May 4 to match recomb_simulation format
#Output:
#IID Generation MotherID FatherID Sex(0=F) ChromParent(0=maternal) Chromosome ChromosomeDiagram
#Input: 
    # N: Number of desired output individuals. will output two haplotypes per individual
    # alpha: Admixture proportion (minor pop)
    # lambda: generations since admixture
    # map: Map file with snp name, chromosome and cM distance: plink format
        #chr rs cmorgans bp 

import sys
import math
from random import randint
from random import uniform
from random import expovariate

def admix(alpha, lam, chr_length):
    #lam is generations since admixture
    #chr_length should be in centimorgans
    chr_length=float(chr_length)/100
    alpha = [0,]+[float(x) for x in alpha]
    for i in range(1, len(alpha)):
        alpha[i] = alpha[i]+alpha[i-1]
    if alpha[-1] > 1:
        sys.exit("Population proportions sum to more than 1")
    lam = float(lam)
    chr_diagram = []; #popid, stop, popid, stop....
   
    #Select starting population
    pop = uniform(0, 1)
    A = alpha+[pop]; A.sort()
    pop=A.index(pop)
    #get resampling points
    #P[Resample ancestry in interval] = 1-exp(-lambda*length)
    pos=0; breaks=[]
    while pos < chr_length: #Get recombination break points
        nb = pos + expovariate(lam)
        if nb < chr_length:
            breaks.append(nb*100)
        pos=nb
    chr_diagram=[pop]
    for b in breaks:
        pop = uniform(0, 1)
        A = alpha+[pop]; A.sort()
        pop=A.index(pop)
        if pop==chr_diagram[-1]:
            continue
        else:
            chr_diagram.append(b)
            chr_diagram.append(pop)
    chr_diagram.append(chr_length*100)
    for i in range(0, len(chr_diagram), 2):
        chr_diagram[i] = str(chr_diagram[i])
    #    chr_diagram[i+1] = "%0.4f" % chr_diagram[i+1]
    return(chr_diagram)

    

if __name__== '__main__':
    from optparse import OptionParser
    parser =OptionParser(usage = "%prog [options] alpha1 alpha2 ... alpha_npops\n")
    parser.add_option("-o", "--out", dest="out", default="out", \
            help="Name an output file.")
    parser.add_option("--chr-lengths", dest="chr_lengths", default="", \
            help="File containing list of chromosome lengths.")
    parser.add_option("--bp-to-cM", default=1000000, type=int, dest="bp_cM", \
            help="Recombination rate per basepair. Default is 1Mb = 1cM. If chromosome lengths are in cM use --bp-to-cM=0 or --bp-to-cM=1." )
    parser.add_option("-N", "--number", type="int", dest="N", \
            help="Number of admixed individuals to produce.", default=1)
    parser.add_option("--npop", type="int", dest="npop", help="Number of populations.", default=1)
    parser.add_option("--lambda", type="int", dest="lam", default=5, help="Number of generations since admixture.")
    (opts, args) = parser.parse_args()
    if len(args) != opts.npop-1:
        sys.exit("Please provide an admixture proportion for each population beyond the first.") 
    if not opts.chr_lengths:
        sys.exit("Chromosome lengths are required.")
    alpha = [float(a) for a in args]
    print("1: "+str(1-sum(alpha)))
    for i in range(len(alpha)):
        print(str(i+2)+": "+str(alpha[i]))
    try:
        inp = open(opts.chr_lengths, "rb")
    except:
        sys.exit("Failed to open "+opts.chr_lengths)
    if opts.bp_cM==0:
        opts.bp_cM = 1
    cms = inp.readlines()
    cms= [float(cm)/opts.bp_cM for cm in cms]
    inp.close()

    try:
        out = open(opts.out, "wb")
        #outfam=open(opts.out+".fam", "wb")
    except:
        sys.exit("Failed to open "+opts.out)
    for i in range(1, opts.N+1):
        for c in range(len(cms)):
            hap1 = admix(alpha, opts.lam, cms[c])
            hap2 = admix(alpha, opts.lam, cms[c])
            if opts.bp_cM > 1:
                for j in range(1, len(hap1), 2):
                    hap1[j] = str(int(round(opts.bp_cM * hap1[j])))
                for j in range(1, len(hap2), 2):
                    hap2[j] = str(int(round(opts.bp_cM * hap2[j])))
            ##FID IID PID MID SEX CHR_ORIGIN CHR CHR_DIAGRAM
            out.write(" ".join([str(i), str(i), "0", "0", "1", "M", str(c+1)]+hap1)+"\n")
            out.write(" ".join([str(i), str(i), "0", "0", "1", "P", str(c+1)]+hap2)+"\n")
        #outfam.write(" ".join([str(i), str(i), "0", "0", "1", "-9"])+"\n")
    out.close()
    #outfam.close()
