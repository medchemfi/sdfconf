#!/usr/bin/env python
# -*- coding: latin-1 -*-

import sys
import re
import math
import numpy
from collections import OrderedDict as OrDi

try:
    import functions
except ImportError:
    from sdfconf import functions

#from sdfconf import functions
#from . import functions

try:
    from future.utils import lmap
except ImportError:
    lmap = map

if sys.version_info[0]==2 and sys.version_info[1]>=7:
    pass
else:
    raise SystemError('Python version must be 2.7. or later, but not 3.x.')


class Mol2File(object):
    def __init__(self, path = None):
        self._molecules = []
        if path:
            self.readfile(path)
    
    def __getitem__(self, ind):
        if isinstance( ind, slice ) or isinstance( ind, int ):
            return self._molecules[ind]
        else:
            raise TypeError("Invalid argument type.")
    
    def readfile(self, path):
        with open(path,'r') as f:
            curflag = None
            for line in f:
                if re.match('@<TRIPOS>',line.strip()):
                    curflag = line[9:].strip()
                    if curflag == 'MOLECULE':
                        self._molecules.append( Mol2Mol() )
                    self._molecules[-1][curflag]=[]
                elif curflag:
                    self._molecules[-1][curflag].append(line.strip('\n'))
    
    
    def writefile(self, path):
        with open(path, 'w') as f:
            for mol in self:
                for field in mol:
                    f.write('@<TRIPOS>'+field+'\n')
                    for line in mol[field]:
                        f.write(line+'\n')
    
    

#End of Mol2File

class Mol2Mol(OrDi):
    def __init__(self, *args, **kwargs):
        super(Mol2Mol, self).__init__(*args, **kwargs)
        self._delims = None
        #self._array = None
        
    def pickatomdata(self, column):
        newdic = OrDi()
        for line in self['ATOM']:
            things = re.split( '\s+', line.strip())
            if len(things)>1:
                if len(things[0])==0:
                    offset = 1
                else:
                    offset = 0
                newdic[functions.numify(things[0+offset])] = things[column+offset]
        return newdic
    
    #def atomsGenerator(self, ignores=['H']):
    def atomsGenerator(self, **kwargs):
        
        gettab = [lambda: int(things[offset+0]), lambda: numpy.array( lmap(functions.numify,things[offset+2:offset+5])) ]
        
        ignores = kwargs.get('ignores',['H'])
        at = kwargs.get('types',False)
        
        if at:
            gettab.append(lambda: things[offset+1].strip() )
        
        if kwargs.get('tabs',False):
            tablist = kwargs.get('tabs')
            if not isinstance(tablist, (list, tuple)):
                tablist = (tablist, )
            gettab[1:]=[lambda x: functions.numify(things[offset+i-1]) for i in tablist]
        
        for line in self['ATOM']:
            things = re.split( '\s+', line.strip())
            if len(things)>1:
                if len(things[0])==0:
                    offset = 1
                else:
                    offset = 0
            else:
                break
            if things[offset+1].strip() in ignores:
                continue
            yield tuple(gette() for gette in gettab)
    
    def injectatomdata(self, data, column, defaultValue=0.0, prec=4):
        
        if data._datastruct == OrDi:
            offset = 0
        else:
            offset = -1
        
        def findindex(linetab, cols):
            realcols = []
            j=-1
            for i, item in enumerate(linetab):
                if len(item.strip())!=0:
                    j+=1
                    if j in cols:
                        realcols.append(i)
                        if len(cols)==len(realcols):
                            return realcols
        
        form = '{:.' + str(prec) +'f}'
        injectinfo = []
        for line in self['ATOM']:
            things = re.split( '(\s+)', line)
            mycol = None
            try:
                ind, mycol = findindex(things, (0, column))
            except ValueError:
                print(things)
                continue
            if mycol:
                try:
                    inject = data._data[functions.numify(things[ind])+offset]
                except KeyError :
                    inject = defaultValue
                if math.isnan(inject) or math.isinf(inject):
                    inject = defaultValue
                stachar = sum( lmap(len, things[:mycol]) )
                endchar = stachar+len(things[mycol])
                injectinfo.append( (form.format(inject), stachar, endchar))
            else:
                continue
        
        sta = min( [item[1] for item in injectinfo] )
        end = max( [item[2] for item in injectinfo] )
        maxlen = max( [ len(item[0]) for item in injectinfo] )
        
        for i, line in enumerate(self['ATOM']):
            self['ATOM'][i] = line[:sta] + (maxlen-len(injectinfo[i][0]))*' ' +  injectinfo[i][0] + line[end+1:]
        
#End of Mol2Mol
