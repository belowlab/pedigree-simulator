#!/usr/bin/python

#Will print Chrom diagram in bp
import sys
from random import expovariate
from random import uniform
#random.expovariate(lambd)
#Exponential distribution. lambd is 1.0 divided by the desired mean. It should be nonzero. 
def recomb(hap1, hap2, rate=1):  
    #mat and pat are chromosome diagrams
    #[fgl1 stop1 fgl2 stop2 ... fgln stopn]
    #Distance between recombinations = exponential(rate)
    haplos=[hap1, hap2]
    chr_length=float(hap1[-1])
    dp = len(str(hap1[-1]).split("."))-1
    if dp==1:
        digits = len(str(hap1[-1]).split(".")[1])
    else:
        digits=0
    lam = 1/float(rate)
    pos=0; breaks=[]
    while pos < chr_length: #Get recombination break points
        nb = round(pos + expovariate(lam), digits)        
        if dp==0:
            nb =int(nb)
        if nb < chr_length:
            breaks.append(nb)
        pos=nb
    start=uniform(0, 1) #choose random chromosome to start on 
    if start < 0.5:
        hap=0
    else:
        hap=1
    newhap=[]
    while len(breaks) > 0:
        nb = breaks.pop(0)
        keep_going=1
        while keep_going == 1: #Continue down chromosme until next break
            fgl=haplos[hap][0]
            stop=haplos[hap][1]
            if stop < nb:
                newhap=newhap+[fgl, stop]
                haplos[hap]=haplos[hap][2:]
            else:
                keep_going=0
                newhap=newhap+[fgl, nb]
        #switch hap
        hap=(hap+1)%2
        fgl=haplos[hap][0]
        stop=haplos[hap][1]
        while stop < nb: #Find position on other haplotype
            haplos[hap]=haplos[hap][2:]
            fgl=haplos[hap][0]
            stop=haplos[hap][1]
    #Fill out rest of chromosome
    newhap=newhap+haplos[hap]
    return(newhap)

def gethap(cur_cd, family, iid, fgl, chr_length, rate):
    #cur_cd dictionary: iid-> [chrom diagram1, chrom diagram2]
    #family: iid ->parents
    #iid
    for parent in [0, 1]: 
        P = family[iid][parent]
        if not len(cur_cd[iid][parent])==0: #already exists
            continue
        if P == "0" : #founder
            newhap=[fgl, chr_length] 
            fgl=fgl+1
            cur_cd[iid][parent] = newhap
            continue
        elif len(cur_cd[P][parent])==0:
            G=gethap(cur_cd, family, P, fgl, chr_length, rate)    
            cur_cd = G[0]; fgl=G[1]
        newhap=recomb(cur_cd[P][0], cur_cd[P][1], rate)
        cur_cd[iid][parent]=newhap
    return([cur_cd, fgl])

if __name__ == '__main__':
    from optparse import OptionParser
    parser =OptionParser(usage = "%prog [options]\n")
    parser.add_option("--fam", dest="fam", default="", \
            help="Name input fam file. (required)")
    parser.add_option("-o", "--out", dest="out", default="out", \
            help="Name an output file.")
    parser.add_option("--chr-lengths", dest="chr_lengths", default="", \
            help="File containing list of chromosome lengths. Default assumes base pair lengths. (required)")
    parser.add_option("--cM", action="store_true", default=False, dest="cM", \
            help="Indicates chromosome lengths are in cM not basepairs.")
    parser.add_option("--bp-to-cM", default=1000000, type=int, dest="bp_cM", \
            help="Recombination rate per base-pair. Default is 1Mb = 1cM or --bp-to-cM=1000000.")
    parser.add_option("--print-cM", default=False, action="store_true", dest="print_cM", \
            help="Print chromosome diagrams in cM. Not yet implemented.")
    (opts, args) = parser.parse_args()
    if (not opts.fam) or opts.chr_lengths=="":
        sys.exit("Usage error. Please give a fam file with --fam and use the --chr-lengths option.")
    try:
        inp = open(opts.fam, "rt")
        cl=open(opts.chr_lengths, "rt")
    except Exception as e:
        print(e)
        sys.exit("recomb_sim.py died: Error opening input file.")
    c_lengths=[l.strip() for l in cl]
    cl.close()
    n = len(c_lengths)+1
    if opts.cM: #store cM as float
        c_lengths=[float(c) for c in c_lengths]
        rate=100
    else: #Store bp as integers
        c_lengths=[int(c) for c in c_lengths]
        rate=100*opts.bp_cM

    fam = dict([(l.split()[1],[l.split()[2], l.split()[3]]) for l in inp]) #PID MID
    inp.seek(0)
    info=dict([(l.split()[1],[l.split()[0], l.split()[4]]) for l in inp]) #FID SEX
    inp.seek(0)
    iids = [l.split()[1] for l in inp]
    inp.close()
 
    out=open(opts.out, "wt") #FID IID PID MID SEX CHR_ORIGIN CHR CHR_DIAGRAM
    for chr in range(1, n):
        fgl=1 
        clength=c_lengths[chr-1]
        cur_cd = dict([(i,[[], []] ) for i in iids])
        for i in iids:
            G=gethap(cur_cd, fam, i, fgl, clength, rate)    
            cur_cd = G[0]; fgl=G[1]
        for i in iids:
            for X in [0, 1]: #M for maternal chrom. P for paternal chrom
                P = ["P", "M"][X]
                out.write(" ".join([info[i][0], i, fam[i][0], fam[i][1], info[i][1], P, str(chr)]+ [str(c) for c in cur_cd[i][X]])+"\n")
    out.close()
