#!/usr/bin/python
 
def get_rel(parents, id1, id2):
    if not parents.has_key(id1):
        print(id1 +" not found.")
        return(1)
    if not parents.has_key(id2):
        print(id2 +" not found.")
        return(1)
    missing = set(["0",])
    rel = "No"
    p1 = set(parents[id1]).difference(missing)
    p2 = set(parents[id2]).difference(missing)
    if id2 in p1 or id1 in p2:
        return("PC")
    elif len(p1.intersection(p2))== 2:
        return("FSib")
    elif len(p1.intersection(p2))== 1:
        if len(p1) == 2 and len(p2) == 2:
            return("HSib")
        else:
            return("Sib")
    else:
        anc1 = p1.copy(); anc1.add(id1)
        anc2 = p2.copy(); anc2.add(id2)
        while len(p1) > 0:
            name = p1.pop()
            if parents.has_key(name):
                anc1 = anc1.union(set(parents[name]))
                p1 = p1.union(set(parents[name]))
        while len(p2) > 0:
            name = p2.pop()
            if parents.has_key(name):
                anc2 = anc2.union(set(parents[name]))
                p2 = p2.union(set(parents[name]))
        anc1= anc1.difference(missing)
        anc2= anc2.difference(missing)
        if len(anc1.intersection(anc2)) > 0:
            rel = "Yes" 
    return rel

def mrca(parents, id1, id2):
    #From parents dictionary, returns a list [mrca1, mrca2] of integers
    #mrca1 is the number of generations back from id1 you find the most recent common ancestor
        #0 is the generation of ID1
    #If ID1 is the parent of ID2 the function will return [0, 1]
    # If the two Ids are sibs or half sibs it will return [1, 1]
    #NA indicates the two ids are unrelated
    if not parents.has_key(id1):
        print(id1 +" not found.")
        return(["ERROR", "ERROR"])
    if not parents.has_key(id2):
        print(id2 +" not found.")
        return(["ERROR", "ERROR"])
    missing = set(["0",])
    gens = {}
    cur_gen1 = set([id1,]).difference(missing)
    cur_gen2 = set([id2,]).difference(missing)
    if id1 == id2:
        print("Ids are the same.\n"+id1)
        return([0, 0])
    gens[id1] = [0, "NA"]
    gens[id2] = ["NA", 0]
    gen_count = 0; mrca = ""
    while cur_gen1 or cur_gen2:
        gen_count = gen_count+1
        new_gen1 = set([]); new_gen2 = set([])
        for n in cur_gen1:
            if parents.has_key(n):
                new_gen1 = new_gen1.union(set(parents[n]))
        for n in cur_gen2:
            if parents.has_key(n):
                new_gen2 = new_gen2.union(set(parents[n]))
        new_gen1 = new_gen1.difference(missing)
        new_gen2 = new_gen2.difference(missing)
        if not (new_gen1 or new_gen2):
            break
        for n in new_gen1:
            if gens.has_key(n):
                gens[n][0] = gen_count
                mrca = n;break
            else:
                gens[n] = [gen_count, "NA"]
        if mrca:
            break
        for n in new_gen2:
            if gens.has_key(n):
                gens[n][1] = gen_count
                mrca = n; break
            else:
                gens[n] = ["NA", gen_count]
        if mrca:
            break
        cur_gen1 = new_gen1.copy(); cur_gen2 = new_gen2.copy()
    #print("Max generations searched: "+str(gen_count))
    if mrca:
        #print(mrca)
        return(gens[mrca])
    else:
        return(["NA", "NA"])
