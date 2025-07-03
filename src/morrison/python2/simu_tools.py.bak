#!/usr/bin/python
#File format from simu (Chris's program modified by Ruth)
#IID Generation MotherID FatherID Sex(0=F) Mat/Pat(0=Mat) Chrom Chrom Diagram
#  0    1           2       3       4       5               6      7...
#File format from recomb_simulation.py (By Jean)
#FID IID FatherID MotherID Sex(1=M,2=F) Mat/Pat(M=Mat,P=Pat) Chrom Chrom Diagram
# 0   1    2          3      4             5                   6      7...
import sys
def clean_simu(file, output_file, perfect=True, ngen=3, fid_present=True):
    #fid_present indicates new file format
    #This function will not be needed with new file format anyway...
    from get_relationship import mrca 
    try:
        inp = open(file, "rb")
        out = open(output_file, "wb")
        if not fid_present:
            out.write(inp.readline()); home_pos = inp.tell() #First line is seed
        else:
            home_pos=0
    except:
        sys.exit("Failed to open "+file+" or "+output_file+".")
    id_index=0; #Jean note: Temporary while transitioning to new file format.
    if fid_present:
        id_index=1
    parents = dict([(l.split()[id_index], l.split()[2:4]) for l in inp])
    inp.seek(home_pos)
    print("Found "+str(len(parents.keys()))+" individuals in "+file)
    print("Cleaning...")
    for l in inp:
        id = l.split()[id_index]
        if not parents.has_key(id):
            continue
        if "0" in parents[id]:
            out.write(l)
            continue
        if not (parents.has_key(parents[id][0]) and parents.has_key(parents[id][1])):
            del(parents[id])
            continue
        par_mrca = mrca(parents, parents[id][0], parents[id][1])
        if (perfect) and (not par_mrca == ["NA", "NA"]):
            del(parents[id])
            continue
        elif (not perfect) and (not par_mrca == ["NA", "NA"]):
            if par_mrca[0] <= ngen or par_mrca[0] <= ngen:
                del(parents[id])
                continue
        out.write(l)
    print("Wrote "+str(len(parents.keys()))+" to "+output_file)
    out.close(); inp.close()
            
def write_fam_file(file, output_file, header=False, FID="0", fid_present=True):
    #fid_present indicates new file format
    #fid_present supercedes FID argument
    #With new file type fam file is simply the first 5 columns
    #Fam file format
    #FID IID PID MID SEX
    inp = open(file, "rb")
    inp.readline()
    outp = open(output_file, "wb")
    if not fid_present:
        inp.readline()
        famdata = dict([(l.split()[0],[l.split()[3], l.split()[2], l.split()[4], "-9"]) for l in inp])
    else:
        famdata = dict([(l.split()[1],[l.split()[0], l.split()[2], l.split()[3], l.split()[4], l.split()[5]]) for l in inp])
    inp.close()
    if header:
        outp.write("FID IID PID MID SEX\n")
    for key in famdata.keys():
        if not fid_present:
            if famdata[key][2]=="0": 
                famdata[key][2]="2" #Rrecode female as 2 and male as 1
            outp.write(" ".join([FID,key]+famdata[key])+"\n")
        else:
            outp.write(" ".join([famdata[key][0],key]+famdata[key][1:])+"\n")
    outp.close()
    
    
def get_ibd (chrom1, chrom2):
    #chrom1 and chrom2 are in the format as output by simu
    fgl1 = [chrom1[2*i] for i in range((len(chrom1))/2)]
    fgl2 = [chrom2[2*i] for i in range((len(chrom2))/2)]
    stop1 = [int(chrom1[2*i+1]) for i in range((len(chrom1))/2)]
    stop2 = [int(chrom2[2*i+1]) for i in range((len(chrom2))/2)]
    start=0; stop=0; gl1=fgl1[0]; gl2=fgl2[0]
    ibd_bp=0; nsegs = 0
    is_ibd = 0 #Flag to help with accurate segment count
    while stop1:
        start = stop
        gl1=fgl1[0]; gl2=fgl2[0]
        if stop1[0] <= stop2[0]:
            stop = stop1.pop(0); gl1 = fgl1.pop(0)
        else:
            stop = stop2.pop(0); gl2 = fgl2.pop(0)
        if gl1 == gl2:
            ibd_bp = ibd_bp+stop - start
            if not is_ibd:
                nsegs = nsegs+1
                is_ibd=1
        else:
            is_ibd=0
    return [ibd_bp, nsegs]

def get_ibd_four(chrom1, chrom2, chrom3, chrom4):
    #assumes first two chromosomes from one ind. 
        #second two chromosomes from second individual
    #will return [ibd1, ibd2, nsegs]
    #inbred ibd will not be recorded separately 
        #(1,2,3)(4) will count as ibd1
        #(12)(34) is ibd 0, (1234) is ibd2
    fgl = []; stop=[]
    fgl.append([chrom1[2*i] for i in range((len(chrom1))/2)])
    fgl.append([chrom2[2*i] for i in range((len(chrom2))/2)])
    fgl.append([chrom3[2*i] for i in range((len(chrom3))/2)])
    fgl.append([chrom4[2*i] for i in range((len(chrom4))/2)])
    stop.append([int(chrom1[2*i+1]) for i in range((len(chrom1))/2)])
    stop.append([int(chrom2[2*i+1]) for i in range((len(chrom2))/2)])
    stop.append([int(chrom3[2*i+1]) for i in range((len(chrom3))/2)])
    stop.append([int(chrom4[2*i+1]) for i in range((len(chrom4))/2)])
    start=0; my_stop=0
    gl = [fgl[i][0] for i in range(0, 4)]
    ibd_bp1=0; ibd_bp2=0; nsegs = 0
    is_ibd = [0, 0, 0, 0] #Flag indicating previous section is ibd. Should lead to more accurate segment count
            #02 03 12 13
    while stop[0]:
        start = my_stop
        gl = [fgl[i][0] for i in range(0, 4)]
        stops=[stop[i][0] for i in range(0, 4)]; x= stops.index(min(stops))
        my_stop = stop[x].pop(0); fgl[x].pop(0)
        if(gl[0] == gl[2]): 
            is_ibd[1] = 0 ; is_ibd[2] = 0#Cannot match two strands to one i.e. (1, 2, 3)(4)
            if not is_ibd[0]:
                nsegs = nsegs+1
                is_ibd[0] = 1
            if (gl[1] == gl[3]): #IBD2
                ibd_bp2 = ibd_bp2 + my_stop - start
                if not is_ibd[3]:
                    nsegs = nsegs+1
                    is_ibd[3] = 1
            else: #IBD1 
                ibd_bp1 = ibd_bp1 + my_stop - start
                is_ibd[3] = 0
        elif (gl[0] == gl[3]):
            is_ibd[0] = 0 ; is_ibd[3] = 0
            if not is_ibd[1]:
                nsegs = nsegs+1
                is_ibd[1] = 1
            if (gl[1] == gl[2]): #IBD2
                ibd_bp2 = ibd_bp2 + my_stop - start
                if not is_ibd[2]:
                    nsegs = nsegs+1
                    is_ibd[2] = 1
            else: #IBD1 
                ibd_bp1 = ibd_bp1 + my_stop - start
                is_ibd[2] = 0
        elif(gl[1] == gl[2]):
            is_ibd[0] = 0; is_ibd[1] = 0; is_ibd[3] = 0
            ibd_bp1 = ibd_bp1 + my_stop - start
            if not is_ibd[2]:
                nsegs= nsegs+1
                is_ibd[2] = 1
        elif(gl[1] == gl[3]):
            is_ibd[0] = 0; is_ibd[1] = 0; is_ibd[2] = 0
            ibd_bp1 = ibd_bp1 + my_stop - start
            if not is_ibd[3]:
                nsegs= nsegs+1
                is_ibd[3] = 1
        else:
            is_ibd = [0, 0, 0, 0]
    return([ibd_bp1, ibd_bp2, nsegs])
