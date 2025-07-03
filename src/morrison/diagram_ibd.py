#!/usr/bin/python
# Modified 9-22-2013 to remove extra options
# Modified 4-18-2012 to print IBD1 IBD2 components
# Modified 4-24 for new format
#File format from simu (Chris's program modified by Ruth)
#IID Generation MotherID FatherID Sex(0=F) Mat/Pat(0=Mat) Chrom Chrom Diagram
#  0    1           2       3       4       5               6      7...
#File format from recomb_sim.py (By Jean)
#FID IID FatherID MotherID Sex(1=M,2=F) Mat/Pat(M=Mat,P=Pat) Chrom Chrom Diagram
# 0   1    2          3      4             5                   6      7...
if __name__ == '__main__':
    import sys
    import simu_tools
    import get_relationship
    from optparse import OptionParser
    parser =OptionParser(usage = "%prog [options] simulated_data.txt\n")
    parser.add_option("--out", dest="o", default="out.txt", help="Output file name.")
    #parser.add_option("--id-list", default = "", dest="idl", help="ID list to include")
    #parser.add_option("--make-clean", default = -1, type="int", dest="ngen", help = "Number of generations to clean pedigrees to.")
    #parser.add_option("--fam-id", default="0", dest="f", help="Family ID to appear in fam file. This is for the old format. Use --old-format option")
    #parser.add_option("--old-format", default=False, action="store_true", dest="of", help="Indicates old format")
    parser.add_option("--chrom_num", default=22, type="int", dest="cn", help="Number of chromosomes")
    (opts, args) = parser.parse_args()

    if  not len(args) == 1:
        sys.exit("One argument required. Use -h for help.")
    data_file = args[0]
    simu_tools.write_fam_file(data_file, args[0]+".fam")

    try:
        input = open(data_file, "rt")
        output = open(opts.o, "wt")
        output.write("ID1 ID2 IBD_Mb IBD1 IBD2 IBD_proportion N_segs Relationship MRCA1 MRCA2\n")
    except Exception as e:
        print(e)
        sys.exit("Failed to open "+data_file+" or "+opts.o)
    
    names = set([(l.split()[1]) for l in input])
    input.seek(0)
    names = list(names)
    print(("\nCalculating pairwise IBD for "+str(len(names))+" individuals ..."))

    parents = dict([(l.split()[1], l.split()[2:4]) for l in input]); input.seek(0)
    print(("Reading "+data_file))
    #Get Chr sizes
    total_length=0.0
    cur_chr=0
    for i in range(1, opts.cn+1):
        while not cur_chr==i:
            l = input.readline()
            cur_chr=int(l.split()[6])
            chr_size = int(l.split()[-1])
        #print(("Chromosome "+str(cur_chr)+" size: "+"%0.4f"% (chr_size/1000000.0)+" Mb."))
        total_length = total_length+chr_size
    #print(("Total size "+"%0.4f"% (total_length/1000000)+" Mb."))
    input.seek(0)

    chroms = {} # id --> [chromosome data1, chrdata2, ...., chrdata22]
    for n in names: #Set chromosome dictionary
        chroms[n] =  [None for i in range(2*opts.cn)]
    for l in input:
        id = l.split(None, 2)[1] 
        if id not in chroms:
            continue
        chr_origin=l.split()[5]
        if chr_origin=="M" or chr_origin=="0":
            chr_origin=0
        else:   
            chr_origin=1
        p=(int(l.split()[6])-1)+(chr_origin*opts.cn)
        chroms[id][p]=l.split()[7:]
    input.close()
    #print(("Done reading "+data_file))
    
    print("Calculating IBD ...")
    counter = 1
    for i in range(len(names)):
        id1 = names[i]
        self_ibd = [0, 0] # [ibd_bp, n_segs]
        for cnum in range(0, opts.cn):
            sibd = simu_tools.get_ibd(chroms[id1][cnum], chroms[id1][opts.cn+cnum])
            self_ibd = [self_ibd[0]+sibd[0], self_ibd[1]+sibd[1]] 
                        #IBD BP                 #IBD NSEGS
        #output.write("ID1 ID2 IBD_Mb IBD1 IBD2 IBD_proportion N_segs Relationship MRCA1 MRCA2\n")
        output.write(" ".join([id1, id1, "%0.4f" % (self_ibd[0]/1000000.0), "NA", "NA", "%0.4f"%(self_ibd[0]/(2*total_length)), str(self_ibd[1]), "SELF", "0", "0"])+"\n")
        if counter % 100 == 0:
            print((str(counter)))
        counter = counter+1
        for j in range(i+1, len(names)):
            id2 = names[j]
            relationship = get_relationship.get_rel(parents, id1, id2)
            if relationship == "No":
                continue #Skip unrelated pairs
            ibd = [0, 0, 0,  0] # [ibd_bp, ibd1_bp, ibd2_bp, n_segs]
            for cnum in range(0, opts.cn):
                my_ibd = simu_tools.get_ibd_four(chroms[id1][cnum], chroms[id1][opts.cn+cnum], chroms[id2][cnum], chroms[id2][opts.cn+cnum])
                ibd[0] = ibd[0] + (my_ibd[0] + 2.0*(my_ibd[1]))
                ibd[1] = ibd[1] + my_ibd[0]; ibd[2]=ibd[2]+my_ibd[1]
                ibd[-1] = ibd[-1] + my_ibd[2]
            relationship = get_relationship.get_rel(parents, id1, id2)
            mrca = get_relationship.mrca(parents, id1, id2)
            output.write(" ".join([id1, id2, "%0.4f" % (ibd[0]/1000000.0), "%0.4f"%(ibd[1]/(total_length)), "%0.4f"%(ibd[2]/(total_length)), "%0.4f"%(ibd[0]/(2.0*total_length)), str(ibd[-1]), relationship, str(mrca[0]), str(mrca[1])])+"\n")
    output.close()
    print("Done!")
    
