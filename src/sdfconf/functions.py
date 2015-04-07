#!/usr/bin/env python
# -*- coding: latin-1 -*-

import os
import sys
import re
#import math
#import argparse
import numpy
#import time
#import copy
#import operator
from collections import OrderedDict as OrDi
#import warnings
import bisect as bi

try:
    from future.utils import lmap
except ImportError:
    lmap = map

if sys.version_info[0]==2 and sys.version_info[1]>=7:
    pass
else:
    raise SystemError('Python version must be 2.7. or later, but not 3.x.')


import sdf
import mol2

def listtostring(tab, fill):
    return ''.join([('{:>'+str(fill)+'}').format(i) for i in tab])

def numify(stri):
    if type(stri) == int:
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
        if string == tststr and re.search('\.\d*$', tststr):
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
        if type(num) in (int,float):
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
            return numpy.asscalar(float(num[0])/numpy.prod(num[1:]))
    except ZeroDivisionError:
        return float('inf')

def remainder(num):
    return numpy.asscalar(numpy.remainder(num[0],num[1]))
    
def mypow(num):
    if type(num) in (tuple, list):
        return numpy.asscalar(numpy.power(*lmap(float,num[:2])))
    else:
        return None

def readcsv(path,sep='\t'):
    f=open(path,'r')
    matrix = csvtomatrix(f.readlines(),sep)
    f.close()
    return matrix

def csvtomatrix(lines,sep):
    chibuti=re.compile('[\s\n;]+')
    spli = re.compile(sep)
    return [[numify(chibuti.sub('',cell)) for cell in spli.split(line)] for line in lines]

def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)


def splitter(stringofparams):
    return parentifier([cell.strip() for cell in stringofparams.split(',')],',')

def parentifier(original,separator):
    stapa=re.compile('\{|\[|\(') #Starting parentesis
    endpa=re.compile('\}|\]|\)') #Ending parenthesis
    
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
    types = set(lmap(type, stuff))
    if len(types) == 1:
        return types.pop()
    elif not str in types:
        return float
    else:
        return None


class InputException(Exception):
    pass

class Findable(object):
    
    
    def __init__(self, iterable=None, **kwargs):
        self._maindict = None #ordered dict jonka alkiot sisältävät koordinaatit listana, indeksi avaimena
        self._orders = [] #Lista joka tulee sisältämää dimensioiden määrän listapareja (tuplessa), jotka sisältävät kunkin koordinaattipisteen indeksin ja yhden dimension koordinaatin. Listat järjestetty kunkin dimension mukaan kasvavaan järjestykseen.
        self._thelen = 0 #Dimensioiden lukumäärä
        if iterable:
            self.set_iterable(iterable, kwargs.get('ignores', ['H'])) #Luodaan datarakenne
                
        else:
            self.set_iterable([]) #Tyhjä rakenne
            
    def set_iterable(self, iterable, ignores=[]):
        if type(iterable) in (list, tuple):
            pass
        elif type(iterable) in (sdf.Sdfmole,mol2.Mol2Mol):
            iterable = [list(atom[1]) for atom in iterable.atomsGenerator(ignores=ignores)]
        else:
            raise TypeError('Input must be of type list, tuple, Sdfmole, Mol2Mol.')
        
        if type(iterable[0]) not in (list, tuple):
            iterable = (iterable,) 
            
        self._thelen = len(iterable[0]) #ks init
        self._maindict = OrDi(enumerate(iterable)) #ks init
        self._orders = [] #ks init
        for i in range(self._thelen): #joka dimensiolle
            self._orders.append(([],[])) #lisätään listapari
            for key in sorted(self._maindict, key=lambda x: self._maindict[x][i]): #jokaiselle pisteelle (koordinaattien lukujärjestyksessä)
                self._orders[i][0].append(  key  ) #listaan 0 lisätään indeksi
                self._orders[i][1].append(  self._maindict[key][i]  ) #listaan 0 lisätään koordinaatti
        
    def findinrange(self, coord, radius, indexes=None , level = 0): #kutsuttaessa anna hakukoordinaatti ja säde, muut auttavat rekursiivisessa dimensioiden käsittelyssä
        '''
        tehdään ensin haut joilla haetaan säteen sisällä olevat pisteen yksittäisten dimensioiden mukaan ja vasta siten löytyneille pisteille tehdään todellinen etäisyyshaku.
        Eli ensin etsitään säteen mukainen kuutio ja vasta sen sisältä pallo :P (kun kolme ulottuvuutta)
        '''
        
        if not indexes:
            indexes = self._maindict.keys() #haetaan alkuperäiset indeksit
        
        ind_left  = bi.bisect_right( self._orders[level][1] , coord[level]-radius ) #puolitushaku piste - radius
        ind_right =  bi.bisect_left( self._orders[level][1] , coord[level]+radius ) #puolitushaku piste + radius
        
        if indexes is not self._maindict.keys(): #mikäli ei ensimmäinen taso
            new_indexes = list( set(indexes) & set(self._orders[level][0][ind_left:ind_right]) ) #ulosmenevistä indekseistä poistetaan ne joita ei ollut edellisessä haussa
        else:
            new_indexes = self._orders[level][0][ind_left:ind_right] #indeksit puolitushakujen välistä löytyvät
        
        if len(new_indexes)==0: #mikäli ei löytynyt mitään, palaa
            return new_indexes
        elif level<self._thelen-1: #mikäli ei olla viimeisessa dimensiossa
            return self.findinrange(coord, radius, new_indexes, level +1) #haetaan myös seuraavasta dimensiosta, annetaan myös täällä löytyneet
        else: #mikäli ollaan viimeisessä dimensiossa
            final = []
            co1 = numpy.array(coord)
            for i in new_indexes: #käydään läpi löydetyt pisteet
                if  (sum((co1-numpy.array(self._maindict[i]))**2))**0.5 <= radius: #etäisyystesti
                    final.append(i)
            return final

