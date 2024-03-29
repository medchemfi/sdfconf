# -*- coding: utf-8 -*-

import os
#import sys
import re
import numpy
from collections import OrderedDict as OrDi
from six import string_types

try:
    from future.utils import lmap
except ImportError:
    lmap = map

'''
if sys.version_info[0]==2 and sys.version_info[1]>=7:
    pass
else:
    raise SystemError('Python version must be 2.7. or later, but not 3.x.')
'''
    
def listtostring(tab, fill):
    return ''.join([('{:>'+str(fill)+'}').format(i) for i in tab])

def numify(stri):
    if type(stri) in (int,long):
        return stri
    elif type(stri) == float:
        if float.is_integer(stri):
            return int(stri)
        else:
            return stri
    else:
        try:
            return int(stri)
        except ValueError:
            try:
                return float(stri)
            except ValueError:
                if stri.lower() in ('null','none'):
                    return float('nan')
                else:
                    return stri
        except TypeError:
            if type(stri)==list:
                return lmap(numify,stri)
            else:
                raise TypeError('Wrong type: {}. {}'.format(str(type(stri)),str(stri)))

def numitest(string):
    num = numify(string)
    if type(num)==float:
        return True
    tststr = str(num)
    del(num)
    if tststr != string:
        tststr = tststr.rstrip('0')
        string = string.rstrip('0')
        if string == tststr and re.search(r'\.\d*$', tststr):
            return True
        string = string.strip('0')
        tststr = tststr.strip('0')
        if len(tststr)==0:
            return False
        else:
            return  string == tststr and tststr[0] == '.'
    else:
        return True

def avg(num):
    try:
        return float(sum(num))/len(num)
    except TypeError:
        if type(num) in (int,float,long):
            return num
        else:
            raise TypeError()

def sub(num):
    if len(num)<2:
        return num
    else:
        return num[0]-sum(num[1:])

def div(num):
    try:
        if len(num)<2:
            return 1.0/num
        else:
            return numpy.asscalar(numpy.array(num[0],dtype=float)/numpy.prod(num[1:]))
    except ZeroDivisionError:
        return float('inf')

def divBCK(num):
    try:
        if len(num)<2:
            return 1.0/num
        else:
            return numpy.asscalar(float(num[0])/numpy.prod(num[1:]))
    except ZeroDivisionError:
        return float('inf')

def myprod(num):
    try:
        if len(num)<2:
            return 1.0*num
        else:
            return numpy.asscalar(numpy.array(num[0],dtype=float)*numpy.prod(num[1:]))
    except ZeroDivisionError:
        return float('inf')

def myprodBCK(num):
    try:
        if len(num)<2:
            return 1.0*num
        else:
            return numpy.asscalar(float(num[0])*numpy.prod(num[1:]))
    except ZeroDivisionError:
        return float('inf')

def remainder(num):
    return numpy.asscalar(numpy.remainder(num[0],num[1]))

def mypowTEST(num):
    try:
        if len(num)<2:
            return numpy.asscalar(numpy.array(num, dtype=float))
        else:
            return numpy.asscalar(numpy.array(num[0],dtype=float)**mypow(num[1:]))
    except ZeroDivisionError:
        return float('inf')
    
def mypowBCK(num):
    if type(num) in (tuple, list):
        return numpy.asscalar(numpy.power(*lmap(float,num[:2])))
    else:
        return None
    
def mypow(num):
    try:
        if len(num)<2:
            return numpy.asscalar( num )
        else:
            return numpy.asscalar( numpy.array(num[0],dtype=float)**numpy.array(num[1],dtype=float) ) #numpy.power(num[1:]))
    except ZeroDivisionError:
        return float('inf')
    
'''        
def myprod(num):
    if type(num) in (tuple, list):
        #return numpy.asscalar(numpy.prod(*lmap(float,num[:2])))
        return numpy.asscalar(numpy.prod(*lmap(float,num)))
    else:
        return None
'''
def readCsv(path,sep='\t'):
    f=open(path,'r')
    matrix = csvtomatrix(f.readlines(),sep)
    f.close()
    return matrix

def csvtomatrix(lines,sep):
    chibuti=re.compile(r'[\s\n;]+')
    spli = re.compile(sep)
    return [[numify(chibuti.sub('',cell)) for cell in spli.split(line)] for line in lines]

def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)


def splitter(stringofparams):
    return parentifier([cell.strip() for cell in stringofparams.split(',')],',')

def parentifier(original,separator):
    stapa=re.compile(r'\{|\[|\(') #Starting parentesis
    endpa=re.compile(r'\}|\]|\)') #Ending parenthesis
    
    m=[]
    s=False
    pars=(stapa,endpa)
    for i, item in enumerate(original):
        sta=list(pars[0].finditer(item))
        end=list(pars[1].finditer(item))
        if len(sta)>len(end):
            if s:
                m[-1]=[i]
            else:
                s=True
                m.append([i])
        elif s:
            if len(sta)<len(end):#end:
                m[-1].append(i)
                s=False
    new=[]
    i=0
    for li in m:
        if len(li)==0:
            continue
        new.extend(original[i:li[0]])
        new.extend([separator.join(original[li[0]:li[1]+1])])
        i=li[1]+1
    new.extend(original[i:])
    if len(original)>len(new):
        return parentifier(new,separator)
    else:
        return new

def lister(string):
    string = string.strip()
    if string[0] == '[' and string[-1] == ']':
        li = splitter(string[1:-1])
        for i, cell in enumerate(li):
            li[i] = lister(cell)
        return li
    else:
        return string
        
def listsep(tosplit, delim):
    '''
    Splits a list by delimiter characters
    Return list and list of delimiters
    not in action: If length is 1, return None
    '''
    splitted = re.split('( *['+delim+'] *)', tosplit)
    delim=splitted[1::2]
    splitted = splitted[0::2]        
    return (splitted,delim)
        
def allsame(listordict):
    '''
    Checks if all items in list, tuple or dict are of same type,
    types being float, int and string. If includes int and float, returns float
    If all same, return type, else return None
    '''
    if type(listordict) == OrDi:
        stuff = listordict.values()
    else:
        stuff = listordict
    #types = set(lmap(type, stuff))
    types = set()
    for thing in stuff:
        if isinstance(thing, string_types):
            types.add(str)
        else:
            types.add(type(thing))
    if long in types:
        types.remove(long)
        types.add(int)
    if len(types) == 1:
        return types.pop()
    elif not str in types:
        return float
    else:
        return None

class InputException(Exception):
    pass
