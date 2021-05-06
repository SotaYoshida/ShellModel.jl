'''
 python script to make snt file "properly ordered"
 V(a,b,c,d,J) -> a<=b & c<=d & a<=c
'''
import copy
import numpy as np
import re #glob,os,sys,time,subprocess,re,copy
import sys

#from multiprocessing import Pool
#import multiprocessing as multi
#import random

def exnum(a) :
    pattern=r'([+-]?[0-9]+\.?[0-9]*)'
    return re.findall(pattern,a)[0]

def make_mydict(p_sps,n_sps):
    pcheck = sorted(p_sps, key=lambda x:(x[1],-x[2],-x[3]), reverse=True)
    ncheck = sorted(n_sps, key=lambda x:(x[1],-x[2],-x[3]), reverse=True)
    mydict ={ }        
    for i,nljtz in enumerate(pcheck):
        mydict[nljtz[0]] = i+1
        pcheck[i][0] = i+1
    for i,nljtz in enumerate(ncheck):
        mydict[nljtz[0]] = i+len(pcheck)+1
        ncheck[i][0] = i+len(pcheck)+1
    return pcheck,ncheck,mydict 
    
def change(label,mydict):
    return [ mydict[label[0]],mydict[label[1]],mydict[label[2]],mydict[label[3]],label[4]]

def orderedTBME(abcdJ,nljs):
    fac = 1.0
    abcdJ=list(map(int,abcdJ))
    totJ = abcdJ[4]
    oabcdJ=[abcdJ[0],abcdJ[1],abcdJ[2],abcdJ[3],totJ]
    if abcdJ[0] > abcdJ[1]:
        oabcdJ[1] = abcdJ[0]
        oabcdJ[0] = abcdJ[1]
        fac = fac * (-1.0)**( (nljs[abcdJ[0]][3] + nljs[abcdJ[1]][3])/2 + totJ + 1)
    if abcdJ[2] > abcdJ[3]:
        oabcdJ[2] = abcdJ[3]
        oabcdJ[3] = abcdJ[2]
        fac = fac * (-1.0)**( (nljs[abcdJ[2]][3] + nljs[abcdJ[3]][3])/2 + totJ + 1)
    if oabcdJ[0] > oabcdJ[2]:
        return [oabcdJ[2],oabcdJ[3],oabcdJ[0],oabcdJ[1],totJ],fac
    elif oabcdJ[0]== oabcdJ[2] and oabcdJ[1] > oabcdJ[3]:
        return [oabcdJ[2],oabcdJ[3],oabcdJ[0],oabcdJ[1],totJ],fac
    else:
        return [oabcdJ[0],oabcdJ[1],oabcdJ[2],oabcdJ[3],totJ],fac


def rewritesnt(inpf):
    inp=open(inpf,"r");lines=inp.readlines();inp.close()
    hit = 0
    count = 0
    p_sps = [];n_sps=[]
    pSPE = []; nSPE=[]
    labels = []; MEs = []
    for line in lines :
        if line.strip()[0]=="!":
            continue
        if hit ==0:
            lp,ln,cp,cn = list(map(int,line.split()))
            hit += 1
            count = 0
            continue
        if hit ==1:
            count += 1
            i,n,l,j,tz = list(map(int,line.split()))
            if count <= lp:
                p_sps += [ [i,n,l,j,tz] ]
            else:
                n_sps += [ [i,n,l,j,tz] ]
            if count == lp+ln:
                hit += 1
            continue
        if hit == 2:
            nsp,zero = list(map(int,line.split()))
            hit += 1
            count = 0
            continue
        if hit == 3:
            count += 1
            i,j,SPE = line.rstrip().split()
            i = int(i); SPE = float(SPE)
            if count <= lp:
                pSPE += [SPE]
            else:
                nSPE += [SPE]
            if count == lp+ln:
                hit += 1
            continue
        if hit == 4:
            tl = line.rstrip().split()
            nTBME,massop,Aref = list(map(int,tl[:-1]))
            p = float(tl[-1])
            hit += 1
            continue
        tl = line.rstrip().split()
        labels += [ list(map(int,tl[:-1])) ]
        MEs += [float(tl[-1])]
    SPEs = pSPE+nSPE
    p_sps,n_sps,mydict = make_mydict(p_sps,n_sps)
    sps = [ ["dummy"] ]+p_sps + n_sps
    nlabels = []
    for label in labels:
        nlabels += [ change(label,mydict) ] 
    olabels = []; nMEs = []
    for i,ME in enumerate(MEs):
        tmp = orderedTBME(nlabels[i],sps)
        olabels += [ copy.copy(tmp[0]) ]
        nMEs += [ tmp[1]* ME]

    ## write snt file
    oup = open("ordered_"+inpf.split("/")[-1],"w")
    print(lp," ",ln,"  ",cp," ",cn,file=oup)
    for tmp in p_sps:
        dum,n,l,j,tz = tmp
        print(str("%3i" % dum), str("%4i" % n),str("%4i" % l),str("%4i" % j),str("%4i" % tz),file=oup)
    for tmp in n_sps:
        dum,n,l,j,tz = tmp
        print(str("%3i" % dum),str("%4i" % n),str("%4i" % l),str("%4i" % j),str("%4i" % tz),file=oup)
    print("  ",nsp,"  ",zero,file=oup)
    for i in range(1,lp+ln+1):
        print(str("%4i" % i),str("%4i" % i),str("%12.6f" % SPEs[mydict[i]-1]),file=oup)
    print(" ",nTBME,"  ",massop, " ",Aref, "  ",p,file=oup)
    for i in range(nTBME):
        a,b,c,d,J = olabels[i]
        print(str("%4i" % a),str("%4i" % b),str("%4i" % c),str("%4i" % d),
              str("%6i" % J),str("%12.6f" % nMEs[i]), file=oup)
    oup.close()


def systematic(targets):
    for target in targets:
        readsnt(target)
    
if __name__ == '__main__':
    args = sys.argv
    if len(args) <= 1:
        print("Usage: python3 make_ordered_snt.py target.snt ")
    inpf = args[1] 
    rewritesnt(inpf)
