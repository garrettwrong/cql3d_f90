#!/usr/bin/python
###########################################################################
## script for converting CQL3D files into parallel form
## usage: doparallel.py in.f out.f insertions.f
## files: patterns.mpi - strings for cleaning (like "write(*,*)")
##        block.mpi    - list of output modules for blocking (YuP: not used now, but can be if needed)
## YuP[07-2016] removed usage of block.mpi as it was inserting lines at wrong place
##        dec.mpi      - list of fortran declaration insructions
###########################################################################


import string, sys

inf = []  ## input lines
out = []  ## output lines

###########################################################################
## clearing all terminal outputs (inf->out)
###########################################################################

pat = [] ## pattern for cleaning
buf = [] ## buffer for continuation lines
str = '' ## joined line

## clean single line 
def cleans(s):
    if s=='':
        return ''
    if s[0] in 'cC!':
        return ''
    state = 0
    r = s[:6]
    for c in s[6:]:
        if state==0:
            if c=='"':
                state = 1
            elif c=="'":
                state = 2
            elif c=='!':
                break
            else:
                r += c
        elif state==1:
            if c=='"':
                state = 0
        elif state==2:
            if c=="'":
                state = 0
    if string.strip(s)=='':
        return ''
    if len(r)<6:
        r += ' '*6
    return r

## saving buffer data into out list (c - prefix)
def printbuf(c=''):
    global buf, out
    if c!='' and string.strip(buf[0][:6])!='':
        ## it is labeled line and it should be commented, so we use new line
        ## with given label and CONTINUE instruction
        out.append(buf[0][:6]+'CONTINUE !***')
    for b in buf:
        if len(b)>0 and b[0] in 'cC':
            out.append(b)
        else:
            out.append(c+b)

## processing buffer data 
def procbuf():
    if buf==[]:
        return
    s = string.join(string.split(str),'')
    for p in pat:
        if p in s:
            printbuf('!MPI :::')
            return
    printbuf()

## main cleaning module (edit patterns.mpi for a list of lines that will be processed->commented)
def clean():
    global buf, str, pat
    ## reading patterns for cleaning
    P = open('mpi/patterns.mpi','r').readlines()
    for p in P:
        pat.append(string.strip(p))
    ## start processing
    buf = []
    for s in inf:
        ss = cleans(s)
        if ss=='': ## empty or comment line
            buf.append(s)
        elif ss=='&':   #ss[5]!=' ': ## continuation line
            # NEED TO COMMENT THE NEXT LINE
            buf.append(s)
            str += ss[6:]
            continue
        else: ## normal line
            procbuf()
            buf = [s]
            str = ss
    procbuf()

###########################################################################
## blocking all outputs from working processes (out->out)
###########################################################################

## detecting endpoint of declaration (YuP: Need to revise for f90)
def start(s, D):   # Used in block() only, which is not important now.
    if s=='':
        return False
    if s[0] in 'cC!':
        return False
    if len(s)>5 and s[5]!=' ':
        return False
    c = string.strip(s[6:])
    if len(c)==0:
        return False
    pos = 0
    com = ''
    while c[pos] in string.ascii_letters:
        com += c[pos]
        pos += 1
        if pos==len(c):
            break
    if string.lower(com) in D:
        return False
    return True

## main blocking module
## YuP[07-2016] removed usage of block.mpi as it was inserting lines at wrong place.
## YuP: Instead, use appropriate lines from mpins_par.f 
##     (such as !MPIINSERT_IF_RANK_NE_0_RETURN) in subroutines which are
##     supposed to be called by mpirank=0 only.
def block():
    ## blocking lines
    ins = ' '*6+"include 'mpilib.h'\n"+' '*6+"if(mpirank.ne.0) return"
    ## checking that current file is in the list of blocked files
    check = False
    for l in open('mpi/block.mpi','r').readlines():
        if string.strip(l)==sys.argv[1]:
            check = True
            break
    if not check: ## current file shouldn't be blocked
        return
    ## reading declaration instructions for detecting inserting point
    D = open('mpi/dec.mpi','r').readlines()
    for i in range(len(D)):
        D[i] = string.strip(D[i])
    ## searching insertion point and inserting blocking instruction
    state = 0
    for i in range(len(out)):
        if state == 0 and start(out[i],D):
            out.insert(i,ins)
            break

###########################################################################
## inserting mpi instructons (out->inf)
###########################################################################

## main inserting module
def dompi():
    global out, inf
    ins = {}
    inf = []
    ## reading insertions
    L = open(sys.argv[3],'r').readlines()
    #print L
    for l in L:
        if (l[0]=='!'):   
            #In *.f90 files, these lines start with !MPI, example: !MPIINSERT_INCLUDE
            key = string.strip(l)
            ins[key] = []
            #print '1. key=',key 
        else:
            ins[key].append(l[:-1])
            #print '2. key=',key
    ## searching and inserting 
    mode = ''
    for l in out:
        if mode=='replace':
            inf.append('!MPI :::'+l)
            mode = ''
            continue
        ll = string.strip(l)
        if (l[0:4]=='!MPI') and ins.has_key(ll):
            inf.append('!MPI >>>')
            for i in ins[ll]:
                inf.append(i)
            inf.append('!MPI <<<')
            if l[0:11]=='CMPIREPLACE':  # Not present anymore (could skip)
                mode = 'replace'
            continue
        inf.append(l)

###########################################################################
## chief function
###########################################################################

def main():
    global inf, out

    ## reading input file into inf buffer
    f = open(sys.argv[1],'r')
    s = f.readline()
    while s!='':
        inf.append(string.rstrip(s))
        s = f.readline()
    f.close()
    
    ## doing parallelizing    
    clean() ## clear all outputs on terminal (write(*,*) and so on). YuP:Deal with it later
    block() ## block all outputs for workers. YuP: not used anymore.
    dompi() ## inserting mpi instructions into code

    ## writing result
    f = open(sys.argv[2],'w')
    for s in inf:
        print >> f, s
    f.close()

## do all
main()
