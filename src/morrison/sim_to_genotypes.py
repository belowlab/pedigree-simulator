#!/usr/bin/env python

#Modified April 2012 to work with new format from recomb_simulation.py
#Modified May 17 to output names of founders used for each FGL

import gzip
import sys
import re

from typing import TextIO
from enum import Enum

class FileFormat(Enum):
    VCF = 1
    BGL = 2

def translate_simu(fgl_dict, cm_list, chrom):
    #fgl_dict: fgl -> [g1, g2, ...]
    #cmlist = [pos1, pos2, ...]
    #Stops in chrom assumed to be in cM
    my_chrom = chrom[:]
    my_cm_list = [float(x) for x in cm_list]
    simu_genos = []
    next_stop = 0
    cur_fgl = "N"
    for i in range(len(my_cm_list)):
        while next_stop < my_cm_list[i]: #time to get the next fgl
            if len(my_chrom) == 0:
                sys.exit("Error in translate_simu:\nChromosome length must be greater than largest snp position")
            cur_fgl = int(my_chrom.pop(0))
            next_stop = my_chrom.pop(0)
            if cur_fgl not in fgl_dict:
                sys.exit("Founder genome label "+cur_fgl+" not found.")
        simu_genos.append(fgl_dict[cur_fgl][i])
    return(simu_genos)

def get_fnames_vcf(f: TextIO):
    f.seek(0)
    for line in f:
        if line.startswith("##"):
            # Skip all meta information lines
            continue
        if line.startswith("#"):
            # Header line. VCF returns 1 record per diplotype
            out_list = []
            for i in line.split()[9:]:
                out_list.append(i)
                out_list.append(i)

            return out_list # We need to double so we get haplotypes

def get_fnames_bgl(f: TextIO):
    f.seek(0)
    return f.readline().split()[2:] #List of haplotype names in this file

def get_snp_names(f: TextIO, file_format: FileFormat) -> list[str]:
    """Returns a list of SNPs from a BGL file or list of positions from VCF"""
    if file_format == FileFormat.VCF:
        snp_names = []
        for line in f:
            # Skip all header lines
            if line.startswith("#"):
                continue
            cols = line.split(None, 9)[:9]
            #cols[7] = "." # Remove INFO, its no longer relevant
            snp_names.append("\t".join(cols))
        f.seek(0)
        for line in f:
            # Return the cursor to end of header
            if line.startswith("#"):
                continue
        return snp_names
    elif file_format == FileFormat.BGL:
        snp_names = [l.split(None, 2)[1] for l in f]
        f.seek(0)
        next(f)
        return snp_names
    else:
        raise Exception(f"Unknown file format for reference population file")


##bgl file format (same as phased):
#I rsid Id1 Id1 Id2 Id2 Id3 Id3 ...
#M rsnum g1 g1 g2 g2 g3 g3 ...
#File format from simu (Chris's program modified by Ruth)
#IID Generation MotherID FatherID Sex(0=F) Mat/Pat(0=Mat) Chrom Chrom Diagram
#  0    1           2       3       4       5               6      7...
#File format from recomb_simulation.py (By Jean)
#FID IID FatherID MotherID Sex(1=M,2=F) Mat/Pat(M=Mat,P=Pat) Chrom Chrom Diagram
# 0   1    2          3      4             5                   6      7...
#Output of admix_sim.py is the same as above but all the fgls correspond to pop. ids.
#Pop. ids are numbered from 1

if __name__ == '__main__':

    from random import shuffle
    from optparse import OptionParser

    parser =OptionParser(usage = "%prog [options] -c N -s simulated_data.txt -m data.map fg1.bgl fg2.bgl ...\n")
    parser.add_option("--mode", type="choice",  choices=("IBD", "admix"), default="IBD", help="To use with output of recomb_sim.py use 'IBD' mode. To use with output of admixture_sim.py use 'admix' mode. Admix mode requires a file with population ID for each individual in the genotype files.") 
    parser.add_option("-c", "--chr", dest="c", default = "", help="Chromosome number. (required)")
    parser.add_option("-s", "--simulation", default = "", dest="s", help="Chromosome diagram output from recomb_sim.py or admix_sim.py. (required)")
    parser.add_option("-m", "--map", default="", dest="m", help="File giving marker positions in cM or bp if using bp-cm-map option. File must contain exactly the snps in the bgl files. (required)")
    parser.add_option("--out", dest="out", default="out.txt", help="Output file name.")
    parser.add_option("--chunk", default = 100, type = "int", dest="chunk", help="Number of SNPs per chunk to process.")
    parser.add_option("--gzip", action="store_true", default = False, dest="gz", help="Indicates genotype files are gzipped")
    parser.add_option("--bp-cm-sim", default=1000000, dest="bptm", help="If stop points are in bp give the conversion. Defaulta assumes bp and uses 1Mb=1cM. For files in cM use --bp-cm-sim=0 ", type="int")
    parser.add_option("--bp-cm-map", type="int", dest="bptm_map", help="If map file is in bp give conversion to cM. Default assumes map file is in bp and uses 1 Mb = 1cM. For map file in cM use --bp-cm-map=0.", default=1000000)
    parser.add_option("--founder-ids", default = "", dest="fids", help="Text file containing list of founder ids to use.")
    parser.add_option("--pop-ids", default="", dest="pid", help="A file with population id for all ids in bgl files. Use this option with 'admix' mode.")
    parser.add_option("--write-names", default="" , dest="write", help="Write out haplotypes used for each FGL")
    parser.add_option("--read-names", default="" , dest="read", help="Provide haplotype assignments from a previous chromosome.")
    parser.add_option("--bind-haplos", default=False, action="store_true", dest="bind", help="Bind founder haplotypes")
    (opts, args) = parser.parse_args()

    #Check options
    if not (opts.c and opts.s and opts.m):
        sys.exit("Improper usage. Use -h for help.")
    if opts.mode == "admix":
        if not opts.pid:
            sys.exit("Improper usage. Use -h for help.")
        if opts.fids:
            sys.exit("--founder-ids used only for IBD mode")
        if opts.read:
            sys.exit("--read-names used only for IBD mode")
        #if opts.write:
        #    sys.exit("--write-names used only for IBD mode")
    if len(args) == 0:
        sys.exit("Improper usage. At least one founder genome file is required. Use -h for help.")
    if opts.fids and opts.read:
        sys.exit("Please don't use both --read-names and --founder-ids")
    
    # Determine BGL or VCF by looking at extension of first file.
    if ".vcf" in args[0].lower():
        file_mode = FileFormat.VCF
    elif ".bgl" in args[0].lower():
        file_mode = FileFormat.BGL
    else:
        sys.exit(f"""Failed to determine file format for file {args[0]}.
Expected vcf or bgl. Ensure that .vcf or .bgl is in file name""")

    if file_mode == FileFormat.BGL:
        # Read map
        # BGL format needs a map.
        try:
            inp = gzip.open(opts.m, "rt")
        except Exception as e:
            print(e)
            sys.exit("Failed to open "+opts.m )
        cm_list = [l.split()[1] for l in inp]
        if opts.bptm_map:
            #print(("Converting map bp positions to cM: "+str(opts.bptm_map)+"=1cM"))
            cm_list = [float(x)/opts.bptm_map for x in cm_list]
        inp.close()
        #print(("Last SNP position: "+str(cm_list[-1])))

        #Read simulated chromosome diagrams 
        try:
            inp = open(opts.s, "rt")
        except Exception as e:
            print(e)
            sys.exit("Failed to open "+opts.s )
        simu_data1 = {} #Chroms for first haplotype, keyed by IID
        simu_data2 = {} #For second haplotype
        chrom_length="N"
        fgls_used = set([])
        for l in inp:
            l = l.split()
            if not l[6] == opts.c:
                continue
            n=l[1] #IID
            chromid=l[5] #M or P
            l=l[7:] #Chromosome diagram
            fgls_used = fgls_used.union(set([int(l[x]) for x in range(0, len(l), 2)]))
            if opts.bptm:
                for i in range(1, len(l), 2):
                    l[i] = float(l[i])/opts.bptm
            if chrom_length == "N":
                chrom_length = float(l[-1])
            elif not chrom_length == float(l[-1]):
                sys.exit("Chromosome lengths are not the same")
            if chromid=="M":
                simu_data1[n]= l[:]
            elif chromid=="P":
                simu_data2[n] = l[:]
            else:
                sys.exit("Error on "+l[0]+" chrmoid "+chromid)
        inp.close()
        #print(("Chromosome length: "+str(chrom_length)))
        names = list(simu_data1.keys())
        if not len(list(simu_data2.keys())) == len(names):
            sys.exit("Something wonky...")
        #print(("I found "+str(len(names))+" simulated individuals."))
        
        #Read provided haplotype correspondences
        if opts.read:
            try:
                inp=open(opts.read, "rt")
            except Exception as e:
                print(e)
                sys.exit("Failed to open"+opts.read)
            set_haplos=dict([l.split() for l in inp])
            inp.close()
            

        #Open data files and read header
        files = [None for _ in args]
        founder_names = [] #List of working names of founder haplotypea from data
                            #If haplotypes are not named uniquely they will be coerced into unique names
        tagending=0 
        coerced_names = {} #Old name -> [newname1, newname2, ...]
        hap_name_set = set([]) #set of unique haplotype names
        for i in range(len(files)):
            if ".vcf" in args[i].lower():
                raise Exception(f"{args[0]} was BGL but {args[i]} detected as VCF") 
            elif ".bgl" in args[i].lower():
                pass
            else:
                sys.exit(f"Failed to determine file format for file {args[i]}. Expected vcf or bgl.")

            try:
                files[i] = gzip.open(args[i], "rt")
            except Exception as e:
                print(e)
                sys.exit("Failed to open and read file "+args[i])

            #List of haplotype names in this file
            fnames = get_fnames_vcf(files[i]) if file_mode == FileFormat.VCF else get_fnames_bgl(files[i])
            #print(("File "+str(i+1)+" has "+str(len(fnames))+ " haplotypes."))
            for j in range(len(fnames)):
                if fnames[j] in hap_name_set: #Name is not unique
                    if fnames[j] in coerced_names:
                        coerced_names[fnames[j]].append(fnames[j]+"_"+str(tagending))
                    else:
                        coerced_names[fnames[j]]=[fnames[j]+"_"+str(tagending)]
                    fnames[j] = fnames[j]+"_"+str(tagending)
                    tagending=tagending+1
                hap_name_set.add(fnames[j])
            founder_names.extend(fnames)
        #print(("I found "+str(len(founder_names))+" founder haplotypes total."))

        # Limit founder names to maximum of 2x names
        if len(founder_names) >= len(names) * 2:
            founder_names = founder_names[:(len(names) * 2)]
            #print(f"Only using {len(names) * 2} to simulate {len(names)} individuals")

        snp_names = get_snp_names(files[0], file_mode)

        if not len(snp_names) == len(cm_list):
            sys.exit("Map file must contain the same snps as bgl files.")
        nchunks = len(snp_names)/opts.chunk + 1
        #print(("There will be "+str(nchunks)+" chunks."))
        #founder_names contains full set of working haplotype names

        #Read population dictionary
        if opts.mode == "admix":
            try:
                inp = open(opts.pid, "rt")
                #Format: Name PopID
            except Exception as e:
                print(e)
                sys.exit("Failed to open "+opts.pid)
            my_fns = set(founder_names)
            pop_ids = {} #population id --> [list of founder names]
            fids=set([])
            for l in inp:
                [n, p] = l.split()
                p = int(p)
                if not n in my_fns:
                    continue
                if p in pop_ids:
                    pop_ids[p].append(n)
                else:
                    pop_ids[p] = [n]
                my_fns.remove(n)
                fids.add(n)
                if n in coerced_names:
                    for new_name in coerced_names[n]:
                        pop_ids[p].append(new_name) 
                        my_fns.remove(new_name)
                        fids.add(new_name)
            if not len(fgls_used - set(pop_ids.keys()))== 0 :
                print((list(pop_ids.keys())))
                print(fgls_used)
                sys.exit("Some of the population ids in the simulated data have no associated founders.")
            #print(("I found populaion information for " +str((len(fids))) + " ids."))
            print((str(len(my_fns)) +" haplotypes have no population information and will be discarded."))
            for population in list(pop_ids.keys()):
                shuffle(pop_ids[population])

        #Read ids to include from --founder-ids option
        if opts.fids:
            try:
                inp = open(opts.fids, "rt")
            except Exception as e:
                print(e)
                sys.exit("Failed to open "+ opts.fids)
            fids = set([])
            for l in inp:
                l=l.split()
                fids = fids.union(set(l))
                for x in l:
                    if x in coerced_names:
                        fids=fids.union(set(coerced_names[x]))
            inp.close()
            print((str(len(fids))+" read from "+opts.fids))
        elif opts.mode=="IBD":
            fids=set(founder_names)
        
        #founder_names contains full set of working haplotype names
        #fids is subset of working haplotype names that will be used
        founder_range = list(range(len(founder_names))) 
        my_fns = fids.copy()
        for i in range(len(founder_names)):
            if not founder_names[i] in fids:
                founder_range.remove(i)
            else:
                my_fns.remove(founder_names[i])

        if opts.read:
            for fgl in list(set_haplos.keys()):
                founder_range[int(fgl)-1]=founder_names.index(set_haplos[fgl])
        elif opts.bind:
            even_range = [founder_range[x] for x in range(0, len(founder_range), 2)]
            even_index = dict([(founder_range[x], x) for x in range(0, len(founder_range), 2)])
            shuffle(even_range)
            newrange =[]
            for k in range(len(even_range)):
                newrange.append(even_range[k])
                indx=even_index[even_range[k]]+1
                newrange.append(founder_range[indx])
            founder_range=newrange[:]
        else:
            shuffle(founder_range)
        founder_names = [founder_names[x] for x in founder_range]
        if opts.write :
            #write out
            #fgl founder_name(working)
            try:
                write_out = open(opts.write, "wt")
            except:
                sys.exit("Failed to open "+opts.write)
            for i in range(len(founder_range)):
                write_out.write(" ".join([str(i+1), founder_names[i]])+"\n")
            write_out.close()
        if not len(my_fns) == 0:
            print((str(len(my_fns))+" Ids not found\n"))
        print((str(len(founder_names))+" founder haplotypes remain."))
        #founder_names now contains only working names to use
        #their file indices are in founder_range
        #index in founder_range will be fgl
        
        #Recode ids in chromosome diagrams for admix mode
        if opts.mode == "admix":
            if opts.write:
                #write out
                #name p:newfgl1,newfgl2 p:newfgl1,newfgl2, ...
                try:
                    write_out = open(opts.write+".admix", "wt")
                except:
                    sys.exit("Failed to open "+opts.write+".admix")
            for n in names:
                if opts.write:
                    write_out.write(n+" ") 
                for population in list(pop_ids.keys()):
                    try:
                        myname1 = pop_ids[population].pop(0)
                        myname2 = pop_ids[population].pop(0)
                    except:
                        sys.exit("Insufficient number of individuals for population "+str(population))
                    new_fgl = str(founder_names.index(myname1)+1)
                    if opts.write:
                        write_out.write(str(population)+":"+new_fgl+",")
                    for j in range(0, len(simu_data1[n]), 2):
                        if simu_data1[n][j] == str(population):
                            simu_data1[n][j] = new_fgl
                    new_fgl = str(founder_names.index(myname2)+1)
                    if opts.write:
                        write_out.write(new_fgl+" ")
                    for j in range(0, len(simu_data2[n]), 2):
                        if simu_data2[n][j] == str(population):
                            simu_data2[n][j] = new_fgl
                if opts.write:
                    write_out.write("\n")
            if opts.write:
                write_out.close()          
        #Open output file and write header
        try:
            output = open(opts.out, "wt")
            output.write("I rsid "+" ".join([n+" "+n for n in names])+"\n")
        except: 
            sys.exit("Error opening "+opts.out)
        
        founder_genomes = dict([(i , []) for i in range(1, len(founder_names)+1)]) #FGLs are integers
        #fgl --> [genotypes]
                    
        counter = 0 #How many snps have I read so far?
        chunk_start = 0
        for snp in snp_names:
            genos = []
            for i in range(len(files)):
                genos.extend(files[i].readline().split()[2:])
            genos = [genos[x] for x in founder_range]
            for i in range(1, len(genos)+1):
                founder_genomes[i].append(genos[i-1])
            counter = counter+1

            if counter % opts.chunk == 0: #Time to print
                #print("Writing chunk "+str(counter/opts.chunk -1 )+" of "+str(nchunks))
                my_genos1 = {}; my_genos2 = {}
                for n in names:
                    my_genos1[n] = translate_simu(founder_genomes, cm_list[chunk_start:counter], simu_data1[n])
                    my_genos2[n] = translate_simu(founder_genomes, cm_list[chunk_start:counter], simu_data2[n])
                for i in range(opts.chunk):
                    output.write(" ".join(["M", snp_names[chunk_start+i]]+[my_genos1[n][i]+" "+my_genos2[n][i] for n in names])+"\n")
                founder_genomes = dict([(i , []) for i in range(1, len(founder_names)+1)]) #Reset fg hash at end of chunk
                chunk_start = counter
        my_genos1 = {}; my_genos2 = {} #Print the last chunk
        for n in names:
            my_genos1[n] = translate_simu(founder_genomes, cm_list[chunk_start:], simu_data1[n])
            my_genos2[n] = translate_simu(founder_genomes, cm_list[chunk_start:], simu_data2[n])
        for i in range(counter-chunk_start):
            output.write(" ".join(["M", snp_names[chunk_start+i]]+[my_genos1[n][i]+" "+my_genos2[n][i] for n in names])+"\n")
        for f in files:
            f.close()
        output.close()

    elif file_mode == FileFormat.VCF:
        # VCF doesn't need a map
        # Positions are in col 2
        try:
            inp = gzip.open(args[0], "rt")
        except Exception as e:
            print(e)
            sys.exit("Failed to open "+args[0] )
        cm_list = []
        for l in inp:
            if l.startswith("#"):
                # Skip all header lines
                continue
            cm_list.append(l.split(None, 2)[1])
        if opts.bptm_map:
            #print(("Converting map bp positions to cM: "+str(opts.bptm_map)+"=1cM"))
            cm_list = [float(x)/opts.bptm_map for x in cm_list]
        inp.close()
        #print(("Last SNP position: "+str(cm_list[-1])))

        #Read simulated chromosome diagrams 
        try:
            inp = open(opts.s, "rt")
        except Exception as e:
            print(e)
            sys.exit("Failed to open "+opts.s )
        simu_data1 = {} #Chroms for first haplotype, keyed by IID
        simu_data2 = {} #For second haplotype
        chrom_length="N"
        fgls_used = set([])
        for l in inp:
            l = l.split()
            if not l[6] == opts.c:
                continue
            n=l[1] #IID
            chromid=l[5] #M or P
            l=l[7:] #Chromosome diagram
            fgls_used = fgls_used.union(set([int(l[x]) for x in range(0, len(l), 2)]))
            if opts.bptm:
                for i in range(1, len(l), 2):
                    l[i] = float(l[i])/opts.bptm
            if chrom_length == "N":
                chrom_length = float(l[-1])
            elif not chrom_length == float(l[-1]):
                sys.exit("Chromosome lengths are not the same")
            if chromid=="M":
                simu_data1[n]= l[:]
            elif chromid=="P":
                simu_data2[n] = l[:]
            else:
                sys.exit("Error on "+l[0]+" chrmoid "+chromid)
        inp.close()
        #print(("Chromosome length: "+str(chrom_length)))
        names = list(simu_data1.keys())
        if not len(list(simu_data2.keys())) == len(names):
            sys.exit("Something wonky...")
        #print(("I found "+str(len(names))+" simulated individuals."))
        
        #Read provided haplotype correspondences
        if opts.read:
            try:
                inp=open(opts.read, "rt")
            except Exception as e:
                print(e)
                sys.exit("Failed to open"+opts.read)
            set_haplos=dict([l.split() for l in inp])
            inp.close()
            

        #Open data files and read header
        files = [None for _ in args]
        founder_names = [] #List of working names of founder haplotypea from data
                            #If haplotypes are not named uniquely they will be coerced into unique names
        tagending=0 
        coerced_names = {} #Old name -> [newname1, newname2, ...]
        hap_name_set = set([]) #set of unique haplotype names
        for i in range(len(files)):
            if ".vcf" in args[i].lower():
                pass
            elif ".bgl" in args[i].lower():
                raise Exception(f"All input files must be the same format. {args[0]} was VCF but {args[i]} is BGL.")
            else:
                raise Exception(f"Failed to determine file format for file {args[i]}. Expected vcf or bgl.")

            files[i] = gzip.open(args[i], "rt")

            #List of haplotype names in this file
            fnames = get_fnames_vcf(files[i])
            #print(("File "+str(i+1)+" has "+str(len(fnames))+ " haplotypes."))
            for j in range(len(fnames)):
                if fnames[j] in hap_name_set: #Name is not unique
                    if fnames[j] in coerced_names:
                        coerced_names[fnames[j]].append(fnames[j]+"_"+str(tagending))
                    else:
                        coerced_names[fnames[j]]=[fnames[j]+"_"+str(tagending)]
                    fnames[j] = fnames[j]+"_"+str(tagending)
                    tagending=tagending+1
                hap_name_set.add(fnames[j])
            founder_names.extend(fnames)
        #print(("I found "+str(len(founder_names))+" founder haplotypes total."))

        # Limit founder names to maximum of 2x names
        if len(founder_names) >= len(names) * 2:
            founder_names = founder_names[:(len(names) * 2)]
            #print(f"Only using {len(names) * 2} to simulate {len(names)} individuals")

        snp_names = get_snp_names(files[0], file_mode)

        if not len(snp_names) == len(cm_list):
            raise Exception("Got different count of physical position vs cm position in VCF file.")
        nchunks = len(snp_names)/opts.chunk + 1
        #print(("There will be "+str(nchunks)+" chunks."))
        #founder_names contains full set of working haplotype names

        #Read population dictionary
        if opts.mode == "admix":
            try:
                inp = open(opts.pid, "rt")
                #Format: Name PopID
            except Exception as e:
                print(e)
                sys.exit("Failed to open "+opts.pid)
            my_fns = set(founder_names)
            pop_ids = {} #population id --> [list of founder names]
            fids=set([])
            for l in inp:
                [n, p] = l.split()
                p = int(p)
                if not n in my_fns:
                    continue
                if p in pop_ids:
                    pop_ids[p].append(n)
                else:
                    pop_ids[p] = [n]
                my_fns.remove(n)
                fids.add(n)
                if n in coerced_names:
                    for new_name in coerced_names[n]:
                        pop_ids[p].append(new_name) 
                        my_fns.remove(new_name)
                        fids.add(new_name)
            inp.close()
            if not len(fgls_used - set(pop_ids.keys()))== 0 :
                print((list(pop_ids.keys())))
                print(fgls_used)
                sys.exit("Some of the population ids in the simulated data have no associated founders.")
            #print(("I found populaion information for " +str((len(fids))) + " ids."))
            print((str(len(my_fns)) +" haplotypes have no population information and will be discarded."))
            for population in list(pop_ids.keys()):
                shuffle(pop_ids[population])

        #Read ids to include from --founder-ids option
        if opts.fids:
            try:
                inp = open(opts.fids, "rt")
            except Exception as e:
                print(e)
                sys.exit("Failed to open "+ opts.fids)
            fids = set([])
            for l in inp:
                l=l.split()
                fids = fids.union(set(l))
                for x in l:
                    if x in coerced_names:
                        fids=fids.union(set(coerced_names[x]))
            inp.close()
            #print((str(len(fids))+" read from "+opts.fids))
        elif opts.mode=="IBD":
            fids=set(founder_names)
        
        #founder_names contains full set of working haplotype names
        #fids is subset of working haplotype names that will be used
        founder_range = list(range(len(founder_names))) 
        my_fns = fids.copy()
        for i in range(len(founder_names)):
            if not founder_names[i] in fids:
                founder_range.remove(i)
            else:
                my_fns.remove(founder_names[i])

        if opts.read:
            for fgl in list(set_haplos.keys()):
                founder_range[int(fgl)-1]=founder_names.index(set_haplos[fgl])
        elif opts.bind:
            even_range = [founder_range[x] for x in range(0, len(founder_range), 2)]
            even_index = dict([(founder_range[x], x) for x in range(0, len(founder_range), 2)])
            shuffle(even_range)
            newrange =[]
            for k in range(len(even_range)):
                newrange.append(even_range[k])
                indx=even_index[even_range[k]]+1
                newrange.append(founder_range[indx])
            founder_range=newrange[:]
        else:
            shuffle(founder_range)
        founder_names = [founder_names[x] for x in founder_range]
        if opts.write :
            #write out
            #fgl founder_name(working)
            try:
                write_out = open(opts.write, "wt")
            except:
                sys.exit("Failed to open "+opts.write)
            for i in range(len(founder_range)):
                write_out.write(" ".join([str(i+1), founder_names[i]])+"\n")
            write_out.close()
        if not len(my_fns) == 0:
            print((str(len(my_fns))+" Ids not found\n"))
        #print((str(len(founder_names))+" founder haplotypes remain."))
        #founder_names now contains only working names to use
        #their file indices are in founder_range
        #index in founder_range will be fgl
        
        #Recode ids in chromosome diagrams for admix mode
        if opts.mode == "admix":
            if opts.write:
                #write out
                #name p:newfgl1,newfgl2 p:newfgl1,newfgl2, ...
                try:
                    write_out = open(opts.write+".admix", "wt")
                except:
                    sys.exit("Failed to open "+opts.write+".admix")
            for n in names:
                if opts.write:
                    write_out.write(n+" ") 
                for population in list(pop_ids.keys()):
                    try:
                        myname1 = pop_ids[population].pop(0)
                        myname2 = pop_ids[population].pop(0)
                    except:
                        sys.exit("Insufficient number of individuals for population "+str(population))
                    new_fgl = str(founder_names.index(myname1)+1)
                    if opts.write:
                        write_out.write(str(population)+":"+new_fgl+",")
                    for j in range(0, len(simu_data1[n]), 2):
                        if simu_data1[n][j] == str(population):
                            simu_data1[n][j] = new_fgl
                    new_fgl = str(founder_names.index(myname2)+1)
                    if opts.write:
                        write_out.write(new_fgl+" ")
                    for j in range(0, len(simu_data2[n]), 2):
                        if simu_data2[n][j] == str(population):
                            simu_data2[n][j] = new_fgl
                if opts.write:
                    write_out.write("\n")
            if opts.write:
                write_out.close()          
        #Open output file and write header FOR VCF
        try:
            output = open(opts.out, "wt")
            files[0].seek(0)
            for line in files[0]:
                if line.startswith("##"):
                    output.write(line)
                else:
                    break

            output.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+"\t".join(names)+"\n")
        except: 
            raise Exception("Error opening and writing header to "+opts.out)
        
        founder_genomes = dict([(i , []) for i in range(1, len(founder_names)+1)]) #FGLs are integers
        #fgl --> [genotypes]
                    
        counter = 0 #How many snps have I read so far?
        chunk_start = 0
        for snp in snp_names:
            genos = []
            for i in range(len(files)):
                #genos.extend(re.split(r'[\t\|]', next(files[i]).rstrip(), (len(names) * 2) + 9)[9:])
                genos.extend(re.split(r'[\t\|]', next(files[i]).rstrip(), maxsplit=(len(names) * 2) + 9)[9:])
            genos = [genos[x] for x in founder_range]
            for i in range(1, len(names)+1):
                founder_genomes[i].append(genos[i-1])
            counter = counter+1

            if counter % opts.chunk == 0: #Time to print
                #print("Writing chunk "+str(counter/opts.chunk -1 )+" of "+str(nchunks))
                my_genos1 = {}; my_genos2 = {}
                for n in names:
                    my_genos1[n] = translate_simu(founder_genomes, cm_list[chunk_start:counter], simu_data1[n])
                    my_genos2[n] = translate_simu(founder_genomes, cm_list[chunk_start:counter], simu_data2[n])
                for i in range(opts.chunk):
                    output.write("\t".join([snp_names[chunk_start+i]]+[my_genos1[n][i]+"|"+my_genos2[n][i] for n in names])+"\n")
                founder_genomes = dict([(i , []) for i in range(1, len(founder_names)+1)]) #Reset fg hash at end of chunk
                chunk_start = counter
        my_genos1 = {}; my_genos2 = {} #Print the last chunk
        for n in names:
            my_genos1[n] = translate_simu(founder_genomes, cm_list[chunk_start:], simu_data1[n])
            my_genos2[n] = translate_simu(founder_genomes, cm_list[chunk_start:], simu_data2[n])
        for i in range(counter-chunk_start):
            output.write("\t".join([snp_names[chunk_start+i]]+[my_genos1[n][i]+"|"+my_genos2[n][i] for n in names])+"\n")
        for f in files:
            f.close()
        output.close()
        print(f"\nDone adding genotypes for chromosome {opts.c}")
