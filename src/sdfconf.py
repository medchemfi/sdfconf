#!/usr/bin/env python
# -*- coding: latin-1 -*-

#imports
import os
import sys
import re
import math
import argparse
import numpy
import time
#import pylab
#import matplotlib.pyplot as pylab
import matplotlib.pyplot as plt
import matplotlib as mpl
import copy
import operator
from collections import OrderedDict as OrDi
import warnings
from types import NoneType

#Common regular expressions used.  
confchop= re.compile('\{\[.+\]\}') #re.compile('\{\[\d+\]\}') #gets conformation number #changed from number to everything

#metachop = re.compile('>\s+<') #Detects beginning of metafield name
stapa=re.compile('\{|\[|\(') #Starting parentesis
endpa=re.compile('\}|\]|\)') #Ending parenthesis
goodchop = re.compile('\s*,{0,1}\s*') #CSV-separator

parre = re.compile('[\{\[\(\)\]\}\"\']') #Matches all parentheses

ckey = "confnum"
atomcut=(0, 10, 20, 30, 31, 34, 36, 39, 42, 45, 48, 51, 54, 57, 60, 63, 66, 69)

class Sdffile(object):
    '''
    Class that represents one .sdf-file. Includes all molecules in nested dictionaries, first by molecule name, then by conformation number.
    Order of molecules (in the .sdf-file) is represented by a list containing cells with molecule name and conformation number.
    '''
    
    molchop = re.compile('^\${4}') #Separates molecules in file
    
    def __init__(self, path=''):
        #Initianilizes an empty file.
        self._dictomoles = dict() # {'aspirin':{1:aspi1, 2:aspi2, ...}, 'bentzene':{1:benz1, 2:benz2, ...}}
        self._orderlist = list()
        if path != '':
            self.xreadself(path)
    
    def __copy__(self):
        #Shallow copy function used by copy.copy()
        new = Sdffile()
        #new._dictomoles = copy.copy(self._dictomoles)
        new._dictomoles = dict()
        for key in self._dictomoles:
            new._dictomoles[key] = copy.copy( self._dictomoles[key] )
        new._orderlist = copy.copy(self._orderlist)
        return new
    
    def __deepcopy__(self,memo):
        #Deep copy function used by copy.deepcopy()
        new = Sdffile()
        new._dictomoles = copy.deepcopy(self._dictomoles,memo)
        new._orderlist = copy.deepcopy(self._orderlist,memo)
        return new
    
    def __iter__(self):
        '''
        New iterator method with an actual generator
        '''
        for molkey in self._dictomoles:
            for confkey in self._dictomoles[molkey]:
                yield self._dictomoles[molkey][confkey]
    
    
    def __len__(self):
        #Return total number of conformations in file.
        return len(self._orderlist)
    
    def __str__(self):
        #Return a single string representing the .sdf-file
        tab = []
        for molinfo in self._orderlist:
            for line in self._dictomoles[molinfo[0]][molinfo[1]].selftolist():
                tab.append(line)
        return ''.join(tab)
        
    def __getitem__(self, ind):
        '''
        return conformation(s) as from a list
        molecules in the same order as in represented file
        '''
        if isinstance( ind, slice ) :
            mols = self._orderlist[ind]
            return [self._dictomoles[inf[0]][inf[1]] for inf in mols]
        elif isinstance( ind, int ) :
            mols = self._orderlist[ind]
            return self._dictomoles[mols[0]][mols[1]]
        else:
            raise TypeError, "Invalid argument type."
        
        
    @staticmethod
    def sdfseparator(strings):
        '''
        Separates a list of strings (as in .sdf-file) into list of lists of strings.
        lists of strings represent single conformations
        return list of lists of strings
        '''
        #Hajoittaa .sdf tiedoston listaksi listoja, joista yksi sisaltaa yhden tiedoston molekyyleista
        _moles = []
        _lines = []
        for line in strings:
            if Sdffile.molchop.match(line):
                _moles.append(_lines)
                _lines = []
            else:
                _lines.append(line)
        del(_lines)
        return _moles
    
    @staticmethod
    def xsdfseparator(xfile):
        '''
        Separates a list of strings (as in .sdf-file) into list of lists of strings.
        lists of strings represent single conformations
        return list of lists of strings
        '''
        #Hajoittaa .sdf tiedoston listaksi listoja, joista yksi sisaltaa yhden tiedoston molekyyleista
        #_moles = []
        _lines = []
        for line in xfile:
            if Sdffile.molchop.match(line):
                #_moles.append(_lines)
                yield _lines
                _lines = []
            else:
                _lines.append(line)
        #del(_lines)
        #return _moles
    
    def add(self, stringsofone):
        '''
        Adds a molecule to datastructure. 
        stringsofone is a list of strings representing a single conformation
        '''
        new = Sdfmole(stringsofone)
        name = new.getname()
        
        if not name in self._dictomoles:
            self._dictomoles[name]=dict()
        n = self.tempconfn(new, name)
            
        self._dictomoles[name][n] = new
        self._orderlist.append([name,n])
    
    def remove(self, name, confn):
        #Remove a molecule by name and conformation number.
        del(self._dictomoles[name][str(confn)])
        self._orderlist.remove([name,confn])
    
    def sdfmetacombi(self, other, bolist=[True, True], overwrite=False):
        '''
        Combine metadata to current file from another sdf-file.
        bolist = [require same name, require same confnumber]
        '''
        for omol in other:
            melist = []
            if bolist[0]:
                try:
                    medict = self._dictomoles[omol.getname()]
                    if bolist[1]:
                        melist = [medict[omol.getconfn()]]
                    else:
                        melist = [medict[key] for key in medict]
                except KeyError:
                    continue
            else:
                if bolist[1]:
                    melist = []
                    for key in self._dictomoles:
                        try:
                            melist.append(self._dictomoles[key][omol.getconfn()])
                        except KeyError:
                            continue
                else: melist = self
            for mol in melist:
                mol.metacombine(omol, overwrite)
                        
    def sdflistremove(self, other, sameconf=True):
        '''
        Cut operator.
        Remove conformations from current file that are present in the other file.
        ''' 
        for omol in other:
            oname = omol.getname()
            if sameconf:
                try:
                    self.remove(oname,omol.getconfn())
                except KeyError:
                    pass
            else:
                if oname in self._dictomoles:
                    del(self._dictomoles[oname])
                    newlist=[]
                    for info in self._orderlist:
                        if info[0]!=oname:
                            newlist.append(info)
                    self._orderlist = newlist
        self.dictmaint()
    
    def stripbutmeta(self, metalogic):
        '''
        removes all atoms from the file, except for those defined in the given logical meta statement.
        reads numbers from metastatement and keeps same atomnumbers in conformations.
        Indexing starts from 1
        '''
        #levels = leveler(metalogic)
        metas = self.getmollogic(metalogic)
        #for mol in self:
        #    mol.stripbutmeta(levels)
        for molinfo in self._orderlist:
            self._dictomoles[molinfo[0]][molinfo[1]].stripbutmeta(metas[molinfo[0]][molinfo[1]])
    
    #conformation numbers
    
    def tempconfn(self, mol, molname=None):
        #Returns a temporary conformation number for given molecule.
        n = mol.getconfn()
        if not n:
            try:
                dictomol = self._dictomoles[molname]
            except KeyError:
                return -1
            i = -len(dictomol)
            while True:
                j=str(i)
                if not j in dictomol:
                    return j
                else:
                    i -= 1
            n = str(i)
        return n
    
    def uniconfn(self, mol, n='-1'):
        #Returns a unique conformation number for given molecule.
        n=abs(int(n))
        try:
            moli = self._dictomoles[mol.getname()]
            while True:
                if str(n) in moli:
                    n += 1
                else:
                    return str(n)
        except KeyError:
            return '1'
    
    def addconfs(self, bolist = [True, False]):
        '''
        add unique conformation numbers to conformations
        bolist=(toName,toMeta confnum)
        ''' 
        for i, molinfo in enumerate(self._orderlist):
            mol = self._dictomoles[molinfo[0]][molinfo[1]]
            c = mol.getconfn()
            if not c:
                c = self.uniconfn(mol, 1)
            mol.addconf(c, bolist)
            del(self._dictomoles[molinfo[0]][molinfo[1]])
            self._dictomoles[molinfo[0]][c]=mol
            self._orderlist[i]=[molinfo[0],c]
            
    def remconfs(self, bolist = [True, True]):
        '''
        remove conformation numbers
        bolist=(fromName,fromMeta confnum)
        '''
        for mol in self:#self._dictomoles:
            mol.remconf(bolist)
    
    #misc
    
    def metatoname(self, meta, joiner='_'):
        '''
        Changes the name of molecules to whatever found in given metafield.
        Applies only for metafields of type str
        '''
        for molord in self._orderlist:
            olname = molord[0]
            n = molord[1]
            mol = self._dictomoles[olname][n]
            
            name = joiner.join(mol.getmeta(meta).getmetastrings())
            mol._name = name
            nn = self.uniconfn(mol)
            self._orderlist[i]=[name,nn]
            del(self._dictomoles[olname][n])
            if not name in self._dictomoles:
                self._dictomoles[name] = dict()
            self._dictomoles[name][nn] = mol
        self.dictmaint()
    
    def metamerger(self, string):
        '''
        Method to merge multiple metafields into one. Comes with a few operations. Frontend to mergenewmeta.
        All but join can also be made with makenewmeta()
        '''
        [newmeta, task] = re.split('\s*=\s*',string)
        funki = task.find('(')
        funk = task[:funki]
        metas = goodchop.split(task[funki:].strip('{[()]}'))
        
        operators = {'sum':sum, 'max':max, 'min':min, 'avg':avg, 'join':metajoiner, 'prod':numpy.prod, 'sub':sub, 'div':div, 'pow':mypow}
        if funk in operators:
            self.mergenewmeta(newmeta, metas, operators[funk.lower().strip()])
        else:
            if args.verbose:
                print 'Command {} not found'.format(funk)
    
    
    def mergenewmeta(self, newmetaname, metas, oper):
        '''Merges multiple metavalues into one'''
        for mol in self:
            newmeta = Sdfmeta()
            newmeta.setname(newmetaname)
            tab= []
            for meta in metas:
                meta = meta.strip()
                if meta in mol._meta:
                    tab.append(mol.getmeta(meta)) #add meta instead of number
                else:
                    try:
                        num =  mol.logicgetmeta(leveler(meta)) #numify(meta)
                    except ValueError:
                        #meta isn't a statement
                        continue
                    if type(num) == Sdfmeta:
                        tab.append(num)
            newmet = Sdfmeta.metaoper(oper, tab)
            if newmet:
                newmet.cleandelim(True)
                mol.addmeta(newmetaname, newmet)
    
    def makenewmetastr(self, string):
        #Frontend for makenewmeta
        string = string.strip()
        comps = ('>=', '<=', '<', '>', '==', '!=', '=' )
        eqind=string.find('=')
        name=string[:eqind].strip()
        ostri = string[eqind+1:].strip()
        for comp in comps:
            i = ostri.find(comp)
            if i>-1:
                self.makenewmeta(name, ostri[:i].strip(), comp, ostri[i+len(comp):].strip() )
                return
        self.makenewmeta(name, ostri.strip(), None, None )
    
    def makenewmeta(self, name, metastatement, logicchar = None, value = None):
        '''
        make a new metafield by metastatement and pick only those fullfilling given logical statement
        '''
        comps = {'>=':operator.ge, '<=':operator.le, '<':operator.lt, '>':operator.gt, '==':operator.eq, '=':operator.eq, '!=':operator.ne }
        newmetas = self.getmollogic(metastatement)
        if logicchar:
            pickvalues = self.getmollogic(value)
        count = 0
        for molname in newmetas:

            for confn in newmetas[molname]:
                if logicchar:
                    try:
                        newmetas[molname][confn].pickvalues( pickvalues[molname][confn], comps[logicchar] )
                    except KeyError:
                        continue
                if len(newmetas[molname][confn])>0:
                    count += 1
                    self._dictomoles[molname][confn].addmeta(name, newmetas[molname][confn], overwrite=True)
        if count<len(self):
            warnings.warn('Not all new metas were generated',UserWarning)
        
    def nametometa(self, meta):
        #Adds a metafield including the molecule name
        for mol in self:
            mol.addmeta(meta, mol.getname(), literal=True)
            
    def changemetaname(self, oldname, newname):
        #Change the name of a metafield
        for mol in self:
            mol.changemetaname(oldname, newname)
            
    def changemeta(self, who, towhat):
        #Change the content of given metafield by logical statement
        state = leveler(towhat)
        for mol in self:
            if who in mol._meta:
                mol._meta[who] = mol.logicgetmeta(state)
    
    def dictmaint(self): 
        #Maintenance of Sdffile dictionaries
        delkeys=[]
        for key in self._dictomoles:
            if len(self._dictomoles[key])==0:
                delkeys.append(key)
        for key in delkeys:
            del(self._dictomoles[key])
            
    def copymetalist(self):
        '''
        Runs metanewdict for all molecules.
        '''
        for mol in self:
            mol.metanewdict()
    
    def sametest(self, bolist, samelist):
        #Test wether molecules are "the same" as described in bolist.
        return (bolist[0] <= samelist[0]) and (bolist[1] <= samelist[1])
        
    def makecsv(self,stringofmetas,separator='\t'):
        #Make a csv-list containing all molecules and given metafields as columns. '?' gives all fields
        listofmeta = [met.strip() for met in re.split('\s*,|;\s*',stringofmetas)]
        if '?' in listofmeta:
            metalist = self.listmetas()
            for meta in listofmeta[1:]:
                metalist.remove(meta)
            listofmeta = metalist
        
        csv = [separator.join(listofmeta)]
        for info in self._orderlist:
            mol = self._dictomoles[info[0]][info[1]]
            line = []#[mol._name]
            for meta in listofmeta:
                if meta in mol._meta:
                    memeta = mol.getmeta(meta) #.getmetastr()
                    if memeta._datastruct != 'single' or  memeta.dtype() == str:
                        line.append('"'+memeta.getmetastr()+'"')
                    else:
                        line.append(mol.getmeta(meta).getmetastr())
                else:
                    line.append('""')
            csv.append(separator.join(line))
        return csv
        
    
    def makeatomiccsv(self,stringofmetas,separator='\t'):
        #Make a csv-list containing all molecules and given metafields as columns. '?' gives all fields
        listofmeta = [met.strip() for met in re.split('\s*,|;\s*',stringofmetas)]
        if '?' in listofmeta:
            metalist = self.listmetas()
            for meta in listofmeta[1:]:
                metalist.remove(meta)
            listofmeta = metalist
        csv = [separator.join(['atom_number'] + listofmeta)]
        for info in self._orderlist:
            mol = self._dictomoles[info[0]][info[1]]
            
            keys=set()
            for meta in listofmeta:
                if not meta in mol._meta:
                    continue
                elif mol.getmeta( meta )._datastruct == OrDi:
                    
                    for key in mol.getmeta( meta )._data.keys():
                        keys.add(key)
            for key in keys:
                newline = [str(key)]
                
                for meta in listofmeta:
                    if not meta in mol._meta:
                        newline.append("NA")
                    elif mol.getmeta(meta)._datastruct == OrDi:
                        newline.append(str(mol.getmeta(meta)._data.get(key,"NA")))
                    else:
                        newline.append(mol.getmeta(meta).getmetastr())
                csv.append(separator.join(newline))
        return csv
        
    def listmetas(self):
        #Make a list of all metafields present in the file.
        metas = []
        for mol in self:
            for meta in mol._metakeys:
                if not meta in metas:
                    metas.append(meta)
        return metas
    
    def counts(self):
        #Get counts of different conformations and molecules
        counts = []
        for mol in self._dictomoles:
            #counts.append('{}\t{}'.format(mol,len(self._dictomoles[mol])))
            counts.append((mol,len(self._dictomoles[mol])))
        return counts
    
    def closestStr(self, string):
        argus = splitter(string)
        if len(argus)==1:
            self.closest(argus[0])
        elif len(argus)==2:
            self.closest(argus[0], name=argus[1])
    
    #should work
    def closest(self, point, **varargdict): #name, intrests, num
        if not 'name' in varargdict:
            name = ''
        else:
            name = varargdict['name']+'_'
        for mol in self:
            #is coord1 meta or not? If it was a metafield with atomnumber, it would give that one as closest, so no.
            alist = mol.dists(point)
            alist=sorted(alist, key=lambda item: item[0])
            #newlist = []
            moles=[]
            if 'intrests' in varargdict:
                switch = {list:varargdict['intrests'],int:[varargdict['intrests']],str:mol.logicgetmeta(leveler(varargdict['intrests']))._data}
                mylist = switch.get( type(varargdict['intrests']), [])
                
                for atom in alist:
                    if atom[1] in mylist:
                        moles.append(atom)
            else:
                moles = alist
            del(alist)
            if 'num' in varargdict:
                moles = moles[:varargdict['num']]
            od = OrDi()
            for item in moles:
                od[item[1]]=item[0]
            mol.addmeta(name+'Closest_atoms',od)
            
    #should work
    def closer(self, point, meta):
        log = leveler(meta)
        for mol in self:
            alist = sorted(mol.dists(point), key  = lambda x: x[0] )
            molmet = mol.logicgetmeta(log)
            for i, atom in enumerate(alist):
                if atom[1] in molmet._data:
                    mol.addmeta('Closer_atoms_than_'+meta.strip(), i)
                    mol.addmeta('Closest_atom_from_'+meta.strip(), atom[1])
                    break

    
    def closestatoms(self, point, metafield):
        dot = Sdfmole.coorder(point)#self.coordormeta(point)
        if type(dot) != list:
            #if not a good point, use argument 'point' as metafield
            dot = False
        else:
            dot = numpy.array(dot)
        for mol in self:
            try:
                if list(dot):
                    mins, mind = mol.pointlistdists(dot,metafield)
                else:
                    mins, mind = mol.atomlistdistances(point,metafield)
            except TypeError:
                continue
            mol.addmeta(metafield+'_closest_distances',' '.join([str(cell) for cell in mind]))
            ind = mind.index(min(mind))
            mol.addmeta(metafield+'_closest_atom',str(mins[ind]))
            mol.addmeta(metafield+'_closest_distance',str(mind[ind]))

    def closestbybonds(self,fromwfield,towfield,newfield):
        for mol in self:
            if towfield in mol._meta and fromwfield in mol._meta:
                home = mol.getmeta(fromwfield)
                target = mol.getmeta(towfield)
                targets = [i.strip() for i in re.split('\s',target)]
                distances = mol.bonddists(home, targets)
                dists = []
                atomos = []
                for tar in targets:
                    dists.append(distances[int(tar)-1])
                    atomos.append(int(tar))
                mind = min(dists)
                mina = atomos[dists.index(mind)]
                mol.addmeta(newfield+'_closest_bond_atom', str(mina))
                mol.addmeta(newfield+'_closest_bond_distance', str(mind))
        
    def coordormeta(self, canditate):
        cords = map(numify,re.split('\s*',canditate))
        if len(cords)==3:
            if  map(type,cords)==[float,float,float]:
                return cords
        else:
            return None
    
    #File operations
    
    def readself(self, path):
        f = open(path, 'r')
        data = f.readlines()
        f.close()
        moles = Sdffile.sdfseparator(data)
        del (data)
        for mole in moles:
            self.add(mole)
            
    def xreadself(self, path):
        with open(path, 'r') as f:
            #xdata = f.xreadlines()
            for mole in Sdffile.xsdfseparator(f):
                self.add(mole)
        
        
    def selftostring(self, output, **kwargs):
        if output=='getcsv':
            return '\n'.join(self.makecsv(kwargs['getcsv']))+'\n'
        elif output=='getatomcsv':
            return '\n'.join(self.makeatomiccsv(kwargs['getatomcsv']))+'\n'
        elif output=='metalist':
            return '\n'.join(self.listmetas())+'\n'
        elif output=='counts':
            towrite=[]
            if kwargs['counts']==0 or kwargs['counts']==2:
                towrite = ['Total {}\t{}'.format(len(self._dictomoles),len(self))]
            if kwargs['counts']==1 or kwargs['counts']==2:
                towrite.extend(('\t'.join(map(str, count)) for count in self.counts()))
            return '\n'.join(towrite)+'\n'
        elif output=='sdf':
            return str(self)
        elif output == 'donotprint':
            return None
        else:
            return None

    def histogrammer(self, Xname, Yname=None, **kwargs):
        if 'ex' in kwargs:
            sdf = copy.copy(self)
            sdf.mollogicparse(kwargs['ex'])
            del(kwargs['ex'])
        else:
            sdf = self
        #kwargs: title, Xtitle, Ytitle, bins
        newargs = dict()
        titargs = ['Xtitle','Ytitle','title']
        for key in kwargs:
            if key in titargs:
                newargs[key]=kwargs[key]
        for key in newargs:
            del(kwargs[key])
        
        datas = [[]]
        almetas = [sdf.getmollogic(Xname)]
        if Yname:
            datas.append([])
            almetas.append(sdf.getmollogic(Yname))
        
        for mol in sdf:
            molname = mol.getname()
            confn = mol.getconfn()
            
            try:
                metas = tuple( almeta[molname][confn] for almeta in almetas )
            except KeyError:
                continue
            if 0 in map(len, metas): #skips 0 len metas
                continue
            try:
                minlen = min( {len(meta) for meta in metas}-{1} ) # shortest list length, ignore singles
            except ValueError:
                minlen = 1
            try:
                structs = set([meta._datastruct for meta in metas])-{'single'}
                #ostru = iter(structs).next() #needed?
            except StopIteration:
                #ostru = list #needed?
                pass
            
            
            if OrDi in structs:
                keys = None
                keyorder = None
                for meta in metas:
                    if meta._datastruct == OrDi:
                        if not keys:
                            keys = set(meta._data.keys())
                        else:
                            keys.intersection_update(set(meta._data.keys()))
                        if not keyorder:
                            keyorder = meta._data.keys()
                        if len(keys) == 0:
                            break
                keyorder = [key for key in keyorder if key in keys][:minlen]
            for i, meta in enumerate(metas): 
                if meta._datastruct == 'single':
                    datas[i].extend(meta._data[:1] * minlen)
                elif meta._datastruct == list:
                    datas[i].extend(meta._data[:minlen])
                elif meta._datastruct == OrDi:
                    try:
                        datas[i].extend([meta._data[key] for key in keyorder])
                    except KeyError:
                        print mol._name
                        print meta._data
                        print keyorder
                        raise KeyError(str(key))
            
        plt.figure()
        if not Yname:
            larg=[datas[0]] #[getcolumn(data,Xname)]
            if 'bins' in kwargs:
                larg.append(kwargs['bins'])
                del(kwargs['bins'])
            #(n, bins, patches)=
            plt.hist(*larg,**kwargs) #kwargs?
        else:
            cmap = mpl.cm.jet
            X=plt.hist2d(datas[0],datas[1],**kwargs)[0]
            
            ticks=list(numpy.arange(numpy.max(X)+1))
            skip = int(round(float(len(ticks))/20))
            if skip>1:
                ticks = numpy.array(ticks[:1] + ticks[skip:-skip:skip] + ticks[-1:])
            norm = mpl.colors.BoundaryNorm(numpy.arange(-0.5,numpy.max(X)+1.5), cmap.N)
            cbar = plt.colorbar()
            cbar = mpl.colorbar.ColorbarBase(cbar.ax , cmap=cmap,norm=norm,ticks=ticks,spacing='proportional')
        
        if 'Xtitle' in newargs:
            plt.xlabel(newargs['Xtitle'])
        else:
            plt.xlabel(Xname)
        if 'Ytitle' in newargs:
            plt.ylabel(newargs['Ytitle'])
        elif Yname:
            plt.ylabel(Yname)
        if 'title' in newargs:
            plt.title(newargs['title'])
    
    def show(self):
        plt.show()
    
    def histogramFromList(self,plots):
        showflag=True
        #plots = args.histogram.split('|')
        for plot in plots:
            path=None
            params = splitter(plot)
            larg = []
            darg = dict()
            for para in params:
                ind=para.find('=')
                if ind!=-1:
                    darg[para[:ind]]=para[ind+1:]
                else:
                    larg.append(para)
            for key in darg:
                darg[key]=lister(darg[key])
                darg[key]=numify(darg[key])
            if 'save' in darg:
                path=darg['save']
                del(darg['save'])
            if 'noplot' in larg:
                showflag=False
                larg.remove('noplot')
            self.histogrammer(*larg,**darg)
            if path:
                plt.savefig(path, bbox_inches='tight')
        if showflag:
            self.show()
        #times.append(time.time())
                
    
    
    def mollogicparse(self, string):
        '''
        remove conformations not fullfilling given meta comparison
        '''
        self.dealer(*self.mollogic(string))
    
    def propor(self, string):
        '''
        return ratios of molecules and conformations fulfilling given metacomparison
        '''
        (trues, falses) = self.mollogic(string)
        porconf = float(len(trues))/len(trues+falses)
        pormols = float(len(set((info[0] for info in trues))))/len(set((info[0] for info in trues+falses)))
        return (pormols, porconf)
    
    def propomes(self, propor=None):
        if propor:
            (molp, confp) = self.propor(propor)
            return '{:.1f}% of molecules and {:.1f}% of conformation fulfill the statement {}.'.format(molp*100, confp*100, propor)
        else:
            return ''
        
    def dealer(self, picks, drops):
        '''
        picks and drops with same format as _orderlist
        remove molecules in drops and keep those in picks
        '''
        #if pick:
        for info in drops:
            del(self._dictomoles[info[0]][info[1]])
            #self._orderlist.remove(info)
        neworder=[]
        for info in self._orderlist:
            if info in picks:
                neworder.append(info)
        self._orderlist = neworder
        #self._orderlist = picks
        self.dictmaint()
        if len(picks)!=len(self):
            warnings.warn('Number of picked ones doesn\'t match number of remaining ones.')
    
    def mollogic(self, string):
        '''
        return lists of trues and falses for meta comparison in given string
        '''
        def compar( opesplit ):
            '''
            basic logical comparison meta1 > meta2, etc. true if statement is true
            '''
            opera = opesplit[1]
            #tocompare = [ leveler(opesplit[0]) , leveler(opesplit[2]) ]
            trues = []
            falses = []
            
            compares = ( self.getmollogic(opesplit[0]), self.getmollogic(opesplit[2]) )
            
            for info in self._orderlist:
                #mole = self._dictomoles[info[0]][info[1]]
                
                try: # added so comparisons to nonexistent metas won't crash
                    #if comps[opera]( mole.logicgetmeta(tocompare[0]), mole.logicgetmeta(tocompare[1]) ):
                    if comps[opera]( compares[0][info[0]][info[1]], compares[1][info[0]][info[1]] ):
                        trues.append(info)
                    else:
                        falses.append(info)
                except (ValueError, KeyError):
                    falses.append(info)
            return (trues, falses)
        
        def mima( string ):
            '''
            pick number, or percentage of conformations for all molecules.
            smallest or greatest of meta for conformations
            '''
            #string = opesplit[0]
            tear = leveler(string)
            if len(tear) != 2:
                raise ValueError('Weird logic.')
            #funk = tear[0]
            matheus = tabjoin(tear[1][1])
            rcomindex = matheus.rfind(',')
            if rcomindex == -1:
                raise TypeError('Weird logic. No comma.')
            metatab = matheus[:rcomindex]
            numstring = matheus[rcomindex+1:]
            perindex = numstring.rfind('%')
            if perindex > 0:
                num = numify(numstring[:perindex])
                per = True
            else:
                num = numify(numstring)
                per = False
            trues = []
            falses = []
            
            if tear[0]=='max':
                reverse = True
                #bymole = False
            elif tear[0]=='min':
                reverse = False
                #bymole = False
            values = self.getmollogic(metatab)
            #for molec in values:
            for mol in self._dictomoles:
                try:
                    molec = values[mol]
                except KeyError:
                    falses.extend([[mol, moln] for moln in self._dictomoles[mol]])
                    warnings.warn(mol + ' marked false for not having necesary metafield.')
                    continue
                if len(set(map(len,[mole[1] for mole in molec.iteritems()])) - {0,1})>0:
                    warnings.warn('Comparing lists, might lead to weird results.',UserWarning)
                moles = OrDi(sorted(molec.iteritems(), key= lambda xx: xx[1], reverse = reverse))
                if per:
                    grab = int(math.ceil(len(moles)*num/100.0))
                else:
                    grab = int(num)
                trues.extend( [[mol, item] for item in moles.keys()[:grab]] )
                falses.extend( [[mol, item] for item in moles.keys()[grab:]] )
            return (trues, falses)
        
        def uniqu( string ):
            '''
            not implemented
            will remove duplicates of conformations by given meta expression
            '''
            string = string[1:]
            pass
        #end of internal functions
        
        string = string.strip()
        
        if string[:2] == '+ ':
            pick = True
            string = string[2:]
        elif string[:2] == '- ':
            pick = False
            string = string[2:]
        else:
            pick = True
        
        opesplit=[item.strip() for item in re.split('(>=|<=|<|>|==|!=|=)',string)]
        comps = {'>=':operator.ge, '<=':operator.le, '<':operator.lt, '>':operator.gt, '==':operator.eq, '=':operator.eq, '!=':operator.ne }
        
        if len(opesplit)==3:
            trues, falses = compar( opesplit )
        elif len(opesplit)==1:
            if opesplit[0][0] == '!':
                trues, falses = uniqu( opesplit[0] )
            else:
                trues, falses = mima( opesplit[0] )
        
        if pick:
            return (trues, falses)
        else:
            return (falses, trues)
    
    def getmollogic(self, string):
        '''
        return single meta for all conformations inside a nested dict like _dictomoles
        meta is made from meta expression in given string
        '''
        molfunx = { 'max':lambda mets: max([met[0] for met in mets]), 'min': lambda mets: min([met[0] for met in mets]), 'avg': lambda mets: avg([met[0] for met in mets]) }
        #fuf = lambda mets: max([met[0] for met in mets])
        sortfunx = { 'asc':True,'des':False }
        
        metafunx = {'mlen':lambda meta: len(meta), 'mavg':lambda meta: avg((thing for thing in meta))}
        
        #def sorto(dire):
        #    pass
        
        maths = {'+':sum,'-':sub,'*':numpy.prod,'/':div,'**':mypow}
        pars = ('(','[','{','"',"'")
        
        mylevel = levels(string)
        metadic = dict()
        precalc = dict()
        
        def tabiter(conf, tab, par = None):
            '''
            collapse given levels tab generated from some meta expression into a single Sdfmeta or None
            '''
            molname = conf.getname()
            
            def listope(conf,tab,par=None):
                #Work with lists in given structures
                if len(tab)==1: #evaluate
                    return tabiter(conf, tab[0], par)
                elif len(tab)==2: #it a slice / asc/max, or something like that
                    if type(tab[0])==str:
                        if tab[0] in molfunx:
                            #do precalc, etc
                            if not tab[0] in precalc or not molname in precalc[tab[0]]:
                                metas = []
                                for confi in self._dictomoles[molname]:
                                    metas.append(tabiter(self._dictomoles[molname][confi],tab[1]))
                                if not tab[0] in precalc:
                                    precalc[tab[0]] = dict()
                                precalc[tab[0]][molname] = Sdfmeta.construct( molfunx[tab[0]](metas) ) #molfunx functions not implemented
                                del(metas)
                            if tab[0] in precalc and molname in precalc[tab[0]]:
                                return precalc[tab[0]][molname] #start and end missing
                        elif tab[0] in sortfunx:
                            #sort next tuple, etc.
                            meta = tabiter(conf, tab[1] )
                            if meta:
                                meta.sortme(sortfunx[tab[0]])
                                return meta
                            else:
                                return None
                        elif tab[0] in metafunx:
                            return Sdfmeta.construct( metafunx[tab[0]](tabiter(conf, tab[1] )) )
                    #slice
                    meta = tabiter(conf, tab[0])
                    if meta:
                        sli = tabiter(conf, tab[1]) #assumes tuple
                        return meta.slicer(sli, tab[1][0])
                    else: return None
                else: #len(tab)>2 evaluate first against the last
                    metaus = tabiter(conf, tab[0])
                    for item in tab[1:]:
                        metaus = tabiter(conf, [metaus,item])
                    return metaus
                
            def tuplope(conf,tab,par=None):
                #Work with tuples in given structure
                if tab[0] in maths:
                    #do math
                    return Sdfmeta.metaoper(maths[tab[0]],[tabiter(conf, tab[1]),tabiter(conf, tab[2])])
                    #pass
                elif tab[0] in pars:
                    #do slicing, etc.
                    return tabiter(conf, tab[1], tab[0]) 
                #return something. or not
                
            def striope(conf,tab,par=None):
                #Work with strings in given structure
                if par in ('"',"'"):
                    return Sdfmeta.construct(tab)
                elif tab in sortfunx or tab in molfunx:
                    return tab
                elif conf.hasmeta(tab):
                    return copy.copy(conf.getmeta(tab)) #FIXME
                elif par in ('(','[','{'):
                    trytab = [numify(i) for i in re.split('\s*[ ,]{0,1}\s*', tab )]
                    if not str in map(type, trytab):
                        return Sdfmeta.construct( trytab )
                    try:
                        return slice(*[{True: lambda n: None, False: int}[x == ''](x) for x in (tab.split(':') + ['', '', ''])[:3]])
                    except ValueError:
                        #return Sdfmeta.construct( tab )
                        raise ValueError('Your logic makes no sense. '+str(tab))
                else:
                    trytab = [numify(i) for i in re.split('\s*[ ,]{0,1}\s*', tab )]
                    if not str in map(type, trytab):
                        return Sdfmeta.construct( trytab )
                    #else: return None
                    
            def raiser(conf, tab, par=None):
                #function for raising error
                raise TypeError('Bad levels: '+str(tab))
            
            testdi = {list: listope, tuple: tuplope, str: striope, Sdfmeta: lambda conf, me, par=None : me if len(me)>0 else None, NoneType : lambda conf, me, par=None : None }
            return testdi.get(type(tab), raiser )(conf, tab, par)
        
        def confloop(mol, tab):
            #Loop through conformations of single type  molecule
            if not mol in metadic:
                metadic[mol] = dict()
            for confnum in self._dictomoles[mol]:
                try:
                    newmet = tabiter(self._dictomoles[mol][confnum], tab)
                except ValueError:
                    newmet = None
                if newmet and len(newmet)>0:
                    metadic[mol][confnum] = newmet
            if len(metadic[mol])==0:
                del(metadic[mol])
        
        def molloop():
            #Loop through types of molecules
            for mol in self._dictomoles:
                confloop(mol, mylevel)
        
        molloop()
        return metadic

    def sorter(self,sortstring):
        '''
        if datastructure defining the sort has more than 1 entry, picks the 1st.
        '''
        sortstring=sortstring.strip()
        sortstring=sortstring.strip('\"|\'')
        if sortstring[0]=='>':
            rever=True
        elif sortstring[0]=='<':
            rever=False
        else:
            return self._orderlist
        return sorted(self._orderlist, key=lambda avain: numify(self._dictomoles[avain[0]][avain[1]].getmeta(sortstring[1:])[0]) if self._dictomoles[avain[0]][avain[1]].hasmeta(sortstring[1:]) else None, reverse=rever)
        
    def sortme(self,sortstring):
        '''
        frontend for sorter
        '''
        self._orderlist = self.sorter(sortstring)
    
    def sortmetas(self, sorty):
        '''
        Sorts cells of a metafield in ascending or descending order
        '''
        sorty=sorty.strip()
        sorty=sorty.strip('\"|\'')
        if sorty[0] == '<':
            ascending=True
            tosort = sorty[1:]
        elif sorty[0] == '>':
            ascending=False
            tosort = sorty[1:]
        else:
            ascending=True
            tosort = sorty
        if sorty[-1] == '+':
            sorty = sorty[:-1]
            byValue = False
        else:
            byValue = True
        for mol in self:
            if tosort in mol._meta:
                mol.getmeta(tosort).sortme(ascending,byValue)
        #return sorted(self._orderlist, key=lambda avain: numify(self._dictomoles[avain[0]][avain[1]]._meta[sorty[1:]][0]),reverse=rever)
    
    def addcsvmeta(self, path, verbose=False):
        '''
        Add metadata from given csv-file path
        '''
        f = open(path)
        chop = re.compile('\s*[;,\t]\s*')
        csvdata = [[cell.strip('\n ') for cell in chop.split(line)] for line in f.readlines()] #lukee tiedoston, splittaa pilkuista ja poistaa alkioista rivinvaihdot
        f.close()
        for j, line in enumerate(csvdata):
            line = parentifier(line, ',')
            s=-1
            newtab=[]
            for i, cell in enumerate(line):
                if s<0 and len(cell)==0:
                    newtab.append('')
                elif s>=0 and len(cell)==0:
                    continue
                elif s < 0 and cell[0] == '"' and cell[-1] != '"':
                    s = i
                    #pass
                elif s >= 0 and cell[0] != '"' and cell[-1] == '"':
                    newtab.append(','.join(line[s:i+1]).strip('"'))
                    s = -1
                elif s>=0:
                    continue
                else:
                    newtab.append(cell.strip('"'))
            csvdata[j] = newtab
        header = csvdata[0]
        for csvmol in csvdata[1:]:
            m = confchop.search(csvmol[0])
            if m:
                n=m.group(0)[2:-2]
                seeker = csvmol[0][:-(len(n)+4)]
            else:
                seeker = csvmol[0]
            if not seeker in self._dictomoles:
                continue
            if m:
                keys = [n]
            else:
                if seeker in self._dictomoles:
                    keys = self._dictomoles[seeker].keys()
                else:
                    continue
            for key in keys:
                if key in self._dictomoles[seeker]:
                    mol = self._dictomoles[seeker][key]
                else:
                    continue
                for i, newmeta in enumerate(csvmol[1:]):
                    if len(newmeta)!=0:
                        mol.addmeta(header[i+1], newmeta, overwrite=True)
                        #mol._meta[header[i+1]].cleandelim(True)
                        mol.getmeta(header[i+1]).cleandelim(True)
        del(csvdata)
    
    def removemeta(self, metaliststring, pick = False):
        #Remove metafields given in list of metanames
        remoli = re.split('\s*,|;\s*',metaliststring)
        if pick:    #keep selected metafields
            for mol in self:
                deles = []
                for meta in mol._meta:
                    if not meta in remoli:
                        deles.append(meta)
                for dele in deles:
                    del(mol._meta[dele])
                    mol._metakeys.remove(dele)
        elif not pick: #remove selected metafields
            for mol in self:
                for remo in remoli:
                    if remo in mol._meta:
                        del(mol._meta[remo])
                        mol._metakeys.remove(remo)
                        
    def splitme(self, n):
        #positive n means how many parts, negative means how many molecules per file
        nofm = len(self)
        
        if n<0:
            n=abs(n)
        elif n>0:
            noff=n
            n=int(math.ceil(float(nofm)/noff))
        files = [Sdffile()]
        count = 0
        for i, mol in enumerate(self):
            files[-1].add(mol.selftolist()[:-1])
            count +=1
            if count >= n:
                files.append(Sdffile())
                count = 0
        if len(files[-1])==0:
            del(files[-1])
        return files
        
        
    def writer(self, writetype, **kwargs):#path=None, split=False, makefolder=False): #new way of writing things
        '''
        Chaotic writer function that puts wanted result in a file/files or standard output
        '''
        newargs = dict(kwargs)
        if writetype=='none':
            return
        if not 'split' in kwargs:
            if 'path' in kwargs: #output or overwrite
                if 'makefolder' in kwargs:
                        del(newargs['makefolder'])
                        #onepath = kwargs['path'][:dotplace]+'_'+str(i)+'/'+kwargs['path']
                        onepath = kwargs['makefolder']+'/'+kwargs['path']
                        ensure_dir(onepath)
                else:
                    onepath = kwargs['path']
                del(newargs['path'])
                f = open(onepath,'w')
                f.write(self.selftostring(writetype, **newargs)) #TODO
                f.close()
                
            else: #stdout
                f = sys.stdout
                f.write(self.selftostring(writetype, **newargs)) #TODO
                f.flush()
        else: #means split
            del(newargs['split'])
            if not 'path' in kwargs: #recursive fix, not possible to make multioutput for stdout
                self.writer(writetype, **newargs) #TODO
            else:
                del(newargs['path'])
                dotplace=kwargs['path'].rfind('.')
                if dotplace < 0:
                    dotplace = len(kwargs['path'])
                files = self.splitme(kwargs['split'])
                for i, onesdf in enumerate(files):
                    onepath = kwargs['path'][:dotplace]+'_'+str(i)+kwargs['path'][dotplace:]
                    if 'makefolder' in kwargs:
                        del(newargs['makefolder'])
                        onepath = kwargs['makefolder']+'_'+str(i)+'/'+kwargs['path']
                        ensure_dir(onepath)
                    newargs['path'] = onepath
                    onesdf.writer(writetype, **newargs) #TODO
    
    def getMol2DataStr(self, string):
        try:
            path, column, metaname = re.split('\s*,\s*', string.strip())
        except ValueError:
            raise ValueError('Give: path, column, metaname')
        self.getMol2Data(metaname, int(column), path)
    
    def getMol2Data(self, metaname, column, path):
        '''
        Get atomwise information from wanted column in given .mol2-file and add it as metadata to current .sdf-file
        If name in mol2 includes sdfconf-compatible confum, data is added conformation wise. Otherwise by complete molecule names.
        '''
        mol2 = Mol2File(path)
        for mol in mol2:
            newdic = mol.pickatomdata(column)
            name = confchop.sub('',  mol['MOLECULE'][0])
            if name != mol['MOLECULE'][0] :
                conf = confchop.search(mol['MOLECULE'][0]).group().strip('{[]}')
            else:
                conf = None
            mols = []
            if name in self._dictomoles:
                if conf and conf in self._dictomoles[name]:
                    mols.append(self._dictomoles[name][conf])
                elif not conf:
                    for key in self._dictomoles[name] :
                        mols.append(self._dictomoles[name][key])
            for sdfmol in mols:
                sdfmol.addmeta(metaname, newdic)
            #self[i].addmeta(metaname, newdic)
        del(mol2)
        
    def injectMol2Data(self, metaname, column, path, defaultValue=0.0, precision=4, outpath = None):
        '''
        Inject atomwise information from current .sdf-file metadata to given wanted column in .mol2-file and make wanted output file.
        '''
        mol2 = Mol2File(path)
        #metaname = metaname.strip() #No logic meta
        level = leveler(metaname)
        for i, mol in enumerate(mol2):
            #mol.injectatomdata(self[i]._meta[metaname], column, defaultValue, precision) #No logic meta
            mol.injectatomdata(self[i].logicgetmeta(level), column, defaultValue, precision)
        if not outpath:
            outpath=path
        mol2.writefile(outpath)
        
    def injectMol2DataStr(self, string):
        try:
            inputf, output, column, metaname, default, precision = re.split('\s*,\s*', string)
        except ValueError:
            raise ValueError('Give: input, output, column, metaname, default, precision')
        self.injectMol2Data(metaname, int(column), inputf, numify(default), int(precision), output)
    
#end of Sdffile

class Sdfmole(object):
    '''
    Class containing all information concerning a single conformation
    '''
    mes='Conformation number mismatch! {} vs {}'
    
    def __init__(self, stringsofone=None):
        '''
        Initialize an empty molecule with no data at all
        Adds data if list of strings containing lines of sdf-file describing single molecule is given
        '''
        self._name = '' #nimi, ei ID
        self._meta = dict() #sisältää metadatan
        self._metakeys = list() #sisältää metadatan avaimet ja pitää niille järjestyksen
        self._other = list() #sisältää kaiken nimen ja ekan metadatan välissä #myöhemmin ainoastaan headerin nimen ja atomien välissä
        self._numeric = False
        self._ignoretype = []
        
        if stringsofone:
            self.initialize(stringsofone)
            
    def __str__(self):
        '''
        Turns the object into a single string as in .sdf-file 
        '''
        return ''.join(self.selftolist())
    
    def __copy__(self):
        '''
        Shallow copy method
        '''
        new = Sdfmole()
        new._name           = self.name #copy.copy(self._name)
        new._meta           = dict(self._meta) #copy.copy(dicself._meta)
        new._metakeys       = list(self._metakeys) #copy.copy(self._metakeys)
        new._numeric        = self._numeric #copy.copy(self._numeric)
        new._ignoretype     = list(self._ignoretype) #copy.copy(self._ignoretype)
        new._comment        = list(self._comment) #copy.copy(self._comment)
        if new.numeric:
            new._counts     = list(self._counts) #copy.copy(self._counts)
            new._atoms      = list(self._atoms) #copy.copy(self._atoms)
            new._bonds      = list(self._bonds) #copy.copy(self._bonds)
            new._properties = list(self._properties) #copy.copy(self._properties)
        else:
            new._other      = list(self._other) #copy.copy(liself._other)
        return new
        
    def __deepcopy__(self,memo):
        '''
        Deep copy method
        '''
        new = Sdfmole()
        new._name           = self._name #copy.deepcopy(self._name,memo)
        new._meta           = copy.deepcopy(self._meta,memo)
        new._metakeys       = list(self._metakeys) #copy.deepcopy(self._metakeys,memo)
        new._numeric        = self._numeric #copy.deepcopy(self._numeric,memo)
        new._ignoretype     = list(self._ignoretype) #copy.deepcopy(self._ignoretype,memo)
        new._comment        = list(self._comment) #copy.deepcopy(self._comment)
        if new._numeric:
            new._counts     = list(self._counts) #copy.deepcopy(self._counts,memo)
            new._atoms      = copy.deepcopy(self._atoms,memo)
            new._bonds      = copy.deepcopy(self._bonds,memo)
            new._properties = list(self._properties) #copy.deepcopy(self._properties,memo)
        else:
            new._other      = list(self._other) #copy.deepcopy(self._other,memo)
        return new
    
    def initialize(self, strings):
        '''
        Add data to molecule from list of strings as in .sdf-file
        '''
        self._name = strings[0].strip()
        self._comment = [line.strip() for line in strings[1:3]] #[strings[i].strip() for i in [1,2]]
        first = -1
        
        for i, line in enumerate(strings[1:]): #notice index shift: i=0; line=strings[1]
            if Sdfmeta.metaname.match(line):#if line[:3] == '> <':
                first = i+1
                break
            else:
                self._other.append(line.rstrip('\n'))
        add = 0
        for i, line in enumerate(strings[first:]): #some optimazition required in the order of statements
            if line == '\n': #emptyline
                if i-add < 2:
                    add = i
                    continue
                newmeta = Sdfmeta(strings[first+add:first+i+1])
                add = i+1
                #self._meta[newmeta.getname()] = newmeta
                self.addmeta(newmeta.getname(), newmeta)
                #self._metakeys.append(newmeta.getname())
                
    def dists(self, point1, ignores=['H']):
        '''
        return list of distances and atom numbers. single line describes distance of an atom in molecule to given point.
        Atom types in ignores are omitted.
        '''
        ignores = [item.upper() for item in ignores]
        self.numerize()
        coord1 = Sdfmole.coorder(point1)
        if coord1 == None:
            coord1 =  self.getatomloc(self.logicgetmeta(leveler(point1+'[0]'))._data[0])
        alist=[]
        for i in range(len(self._atoms)):
            #this responsibility will be taken elsewhere
            if self.gettype(i+1).strip() in ignores:
                continue
            
            #coord2 = numpy.array(self.getcoord(i+1)) #numpy.array([float(c) for c in atom[0:3]])
            coord2 = self.getatomloc(i+1)
            dist = sum((coord1-coord2)**2)**0.5
            alist.append([float(dist),i+1])
        return alist

    def numerize(self):
        '''
        Turn initial textual format of the atomblock to listlike and numerical format.
        '''
        if not self._numeric:
            self._counts = [self._other[2][i:i+3].strip() for i in range(33)[::3]]+[self._other[2][33:39].strip()]
            blockdiv = [3,3+int(self._counts[0]),3+int(self._counts[0])+int(self._counts[1])]
            self._atoms = [[atom[atomcut[i]:atomcut[i+1]] for i in range(len(atomcut)-1)] for atom in self._other[3:blockdiv[1]]]
            self._bonds = [[bond[i:i+3].strip() for i in range(7*3)[::3]] for bond in self._other[blockdiv[1]:blockdiv[2]]]
            self._properties = self._other[blockdiv[2]:]
            del(self._other) #DANGER!!!
            self._numeric = True
        
    def getconf(self):
        '''
        return conformation number in the name and meta for the conformation
        '''
        ans = [None, None]
        m = confchop.search(self._name)
        if m:
            ans[0] = m.group(0)[2:-2]
        if ckey in self._meta:
            ans[1] = self.getmeta(ckey)[0]
        if ans[0] and ans[1]:
            #if ans[0]!=str(ans[1]):
            if str(ans[0])!=str(ans[1]):
                print(self.mes.format(ans[0],ans[1]))
        return ans
        
    def getconfn(self):
        '''
        Return single conformationnumber
        '''
        a = self.getconf()
        return a[0] or a[1]
        
    def addconf(self, conf, bolist = [True, True]):
        '''
        add given conformation number
        bolist =(toName, toMeta confnum)
        '''
        conf = str(conf)
        already =  self.getconf() 
        if bolist[1]:
            if already[1]:
                if conf != str(already[1]):
                    print(self.mes.format(already[0],already[1]))
            else:
                self.addmeta(ckey, conf, literal = True)
        if bolist[0]:
            if already[0]:
                if already[0] != conf:
                    print(mes.format(ans[0],ans[1]))
            else:
                self._name = self._name + '{[' + conf + ']}'
                
    def remconf(self, bolist = [True, True]):
        '''
        remove conformation number
        bolist =(fromName, fromMeta confnum)
        '''
        already = self.getconf()
        if bolist[0] and already[0]:
            self._name = confchop.sub('',self._name)
        if bolist[1] and already[1]:
            self._metakeys.remove(ckey)
            self._meta.pop(ckey)
            
    def getname(self):
        '''
        return name of molecule without possible conformation number
        '''
        return confchop.sub('', self._name)
    
    def gettype(self, N):
        '''
        return atom type of given atom
        Indexing from 1
        '''
        return self._atoms[N-1][4]
    
    def hasmeta(self, metaname):
        '''
        return boolean value if conformation has given metafield
        '''
        return metaname in self._meta
    
    def getmeta(self, metaname, **kwargs):
        '''
        return Sdfmeta with given metaname.
        wont numerize meta if argument "dummy" is given
        '''
        if self.hasmeta(metaname):
            if 'dummy' in kwargs and  kwargs['dummy']:
                pass
            else:
                self._meta[metaname].numerize()
            return self._meta[metaname]
        else:
            return None
    
    def metacombine(self, othersdfmol,overwrite=False):
        '''
        Add all Sdfmeta from another Sdfmole
        '''
        for key in othersdfmol._metakeys:
            if not overwrite and key in self._meta:
                continue
            self.addmeta(key,othersdfmol.getmeta(key,dummy=True),overwrite=True)
                
    def issame(self, othersdfmol):
        '''
        return boolean values if molecule names and conformation numbers are same
        '''
        same = [confchop.sub('', self._name) == confchop.sub('', othersdfmol._name), self.getconfn() == othersdfmol.getconfn()]
        return same
    
    def make_other(self):
        '''
        return reconstructed atomblock from numerized molecule
        '''
        other = []
        for line in self._comment:
            other.append(line + '\n')
        other.append(listtostring(self._counts[:-1],3)+'{:>6}'.format(self._counts[-1])+'\n') # atomheader!
        for atom in self._atoms:
            other.append(Sdfmole.atomtostring(atom)+'\n')
        for bond in self._bonds:
            other.append(listtostring(bond, 3)+'\n')
        for prop in self._properties:
            other.append(prop+'\n')
        return other
    
    def selftolist(self):
        '''
        return constructed list of strings describing the conformation as in .sdf-file
        '''
        me=[]
        me.append(self._name + '\n')
        
        if self._numeric:
            me.extend(self.make_other())
        else:
            for line in self._other:
                me.append(line + '\n')
        
        for key in self._metakeys:
            mystrings = [line+'\n' for line in self.getmeta(key,dummy=True).selftolistofstrings()]
            me.extend(mystrings)
        me.append('$$$$\n')
        return me
        
    #def addmeta(self,metafield, value, overwrite = True):
    
    def addmeta(self,metafield, value, **dictarg): #overwrite, literal
        '''
        Create a new Sdfmeta with given name and value and add it to the conformation
        '''
        if 'overwrite' in dictarg and dictarg['overwrite']:
            overwrite = True
            del(dictarg['overwrite'])
        else:
            overwrite = False
        
        if metafield in self._meta and not overwrite:
            return
        if isinstance(value, Sdfmeta): #type(value) == Sdfmeta:
            insert = value
        else:
            insert = Sdfmeta.construct(value, name = metafield, **dictarg)
        self._meta[metafield] = insert
        if not self.getmeta(metafield,dummy=True).isdumb():
            self.getmeta(metafield).cleandelim(True)
        self.getmeta(metafield,dummy=True).setname(metafield)
        
        if not metafield in self._metakeys:
            self._metakeys.append(metafield)
        
    def metasort(self, reverse=False):
        '''
        Sort the order of metafields appearing in .sdf-file
        '''
        self._metakeys.sort(reverse, key=lambda meta: meta.lower())
        
    def changemetaname(self, oldname, newname):
        '''
        Change name of given Sdfmeta
        '''
        if oldname in self._meta:
            i=self._metakeys.index(oldname)
            self._metakeys[i]=newname
            self.addmeta(newname, self.getmeta(oldname, dummy=True))
            #self._meta[newname]=self._meta[oldname]
            #self._meta[newname].setname(newname)
            del(self._meta[oldname])
            
    def metanewdict(self):
        '''
        Makes a copy of dict holding the metadata. Assists in case you have multiple sdffiles and want to make new metas with same names.
        '''
        self._meta = dict(self._meta)
        self._metakeys = list(self._metakeys)
            
    def atomlistdistances(self,metaatom,metalist):
        '''
        used to give distances from one atom in metafield to many atoms in metafield
        in future will:
        add new meta with 
        '''
        self.numerize()
        try:
            atom=numify(self.getmeta(metaatom))
            manyatom=self.getmeta(metalist)
        except KeyError:
            return None
        atoms = [numify(one) for one in re.split('\s+',manyatom.strip())]
        dists = []
        for oneatom in atoms:
            dists.append(self.atomdist(atom,oneatom))
        return (atoms,dists)
    
    @staticmethod
    def atomtostring(tab):
        d=[]
        for i in range(len(atomcut)-1):
            n=atomcut[i+1]-atomcut[i]
            if n==10:
                d.append('{:>10.4f}'.format(float(tab[i])))
            elif n==1:
                d.append(' ')
            elif i==4:
                d.append('{:<3}'.format(tab[i]))
            else:
                d.append(('{:>'+str(n)+'}').format(tab[i]))
        return ''.join(d)
        
    def getatomloc(self,n):
        '''
        return coordinates of given atom
        '''
        self.numerize()
        return numpy.array([numify(coord) for coord in self._atoms[n-1][0:3]])
    
    def atomdist(self, atom1, atom2):
        '''
        return distance between 2 atoms
        '''
        self.numerize()
        return sum((self.getatomloc(atom1)-self.getatomloc(atom2))**2)**0.5
    
    def pointlistdists(self, point, metalist):
        '''
        return lists of atoms and distances from given point to atoms in given metalist
        '''
        alldists = self.dists(point)
        try:
            atoms = self.getmeta(metalist)
        except KeyError:
            return None
        dists = []
        for oneatom in atoms:
            dists.append(alldists[oneatom-1])
        return (atoms,dists)
    
    
    def bonddists(self, atom1, intrests=None):
        '''
        Doesn't work right now
        #dijikstra algorithm for finding path length to all atoms, or to list of intresting atoms
        '''
        self.numerize()
        n=len(self._atoms)
        visited = list(numpy.zeros(n,dtype=bool))
        distances = list(numpy.ones(n,dtype=int)*numpy.inf)
        atom=numify(atom1)
        distances[atom-1]=0
        if intrests:
            intrests = map(numify, intrests)
        while True:
            links = []
            for bond in self._bonds:
                if numify(bond[0])==atom:
                    links.append(numify(bond[1]))
                if numify(bond[1])==atom:
                    links.append(numify(bond[0]))
            for target in links:
                if distances[target-1]>distances[atom-1]+1:
                    distances[target-1]=distances[atom-1]+1
            visited[atom-1]=True
            
            if intrests:
                if atom in intrests:
                    intrests.remove(atom)
                    if len(intrests)==0:
                        for i,vis in enumerate(visited):
                            if not vis:
                                distances[i]=numpy.inf
                        break
            
            d = numpy.inf
            nextone = -1
            for i,nex in enumerate(distances):
                if nex<d and not visited[i]:
                    d=nex
                    nextone = i+1
            if nextone==-1:
                break
            else:
                atom=nextone
        return distances
    
    def pickmetavalues(self, meta, value, operator,copyme=False):
        '''
        Removes entries in single Sdfmeta that won't fullfill given comparison
        '''
        if self.hasmeta(meta):
            if type(value)==str:
                comps = {'>=':operator.ge, '<=':operator.le, '<':operator.lt, '>':operator.gt, '==':operator.eq, '=':operator.eq, '!=':operator.ne }
                operator = comps.get(operator.strip())
            newmeta = self.getmeta(meta)
            if copyme:
                newmeta = copy.copy(newmeta)
            newmeta.pickvalues(value, operator)
            self._meta[meta] = newmeta
    
    def logicgetmetastr(self, string):
        #return meta defined in given meta expression
        return self.logicgetmeta(leveler(string))
    
    def logicgetmeta(self, partab):
        #return meta defined in given list structure
        return self.collapser(self.levopemap(partab))
        
    def levopemap(self, tab, par=None, length=None):
        '''
        partly deprecated
        Accepts lists made by leveler method. Maps this list for mathematical operators and metafield names. Also slices.
        '''
        if isinstance(tab, tuple):
            #if par and par in ('"',"'"):
            if tab[0] in ('"',"'"):
                return Sdfmeta.construct(levjoin(tab[1]))
                #return levjoin(tab[1])
            else:
                return (tab[0], self.levopemap(tab[1], tab[0]))
        elif isinstance(tab, list):
            for j, thing in enumerate(tab):
                if not isinstance(thing, str):
                    continue
                thing = thing.strip()
                #new
                sear = re.search('[+-]', thing)
                if sear:
                    i = sear.start()
                    ie = sear.end()
                if not sear:
                    sear = re.search('(^|[^\*])[\*/]([^\*]|$)', thing)
                    if sear:
                        helpse = re.search('[\*/]', thing[sear.start():sear.end()])
                        i = sear.start()+helpse.start()
                        ie = sear.start()+helpse.end()
                        del(helpse)
                if not sear:
                    sear = re.search('[\*]{2}', thing)
                    if sear:
                        i = sear.start()
                        ie = sear.end()
                if sear:
                    tab1=[]
                    tab1.extend(tab[:j])
                    if i > 0:
                        tab1.append(thing[:i])
                    if ie < len(thing):
                        tab2 = [thing[ie:]]
                    else:
                        tab2 = []
                    tab2.extend(tab[j+1:])
                    return (thing[i:ie],self.levopemap(tab1),self.levopemap(tab2))
                #new ends
            return [self.levopemap(item, par,len(tab)) for item in tab]
        else:
            
            tab = tab.strip()
            
            if self.hasmeta(tab):
                return self.getmeta(tab)
            
            if length==1 and par!=None:
                return str(tab)
            else:
                tab = numify(tab)
                if not isinstance(tab, str):
                    return Sdfmeta.construct(tab)
                else:
                    raise ValueError('Your logic makes no sense: '+tab)
                    #return Sdfmeta.construct(tab)
    
    @staticmethod
    def coorder(point):
        coord = numpy.array( map(numify, goodchop.split(point.strip('{[()]}'))) )
        if type(coord) == numpy.ndarray and len(coord) > 1:
            return coord
        else:
            return None
    
    def collapser(self, tab, para=None):
        '''
        deprecated
        collapses given structure to single Sdfmeta
        '''
        maths = {'+':sum,'-':sub,'*':numpy.prod,'/':div,'**':mypow}
        pars = ('(','[','{')
        if isinstance(tab,list):
            if len(tab)==2:
                #it's a "slice"
                meta = self.collapser(tab[0])
                sli = self.collapser(tab[1][1],tab[1][0])
                return meta.slicer(sli, tab[1][0]) #THIS IS IT!
                
            elif len(tab)==1:
                #evaluate
                return self.collapser(tab[0],para)
            elif len(tab)>2:
                metaus = self.collapser(tab[0])
                for item in tab[1:]:
                    metaus = self.collapser([metaus,item])
                return metaus
        elif isinstance(tab, tuple):
            if tab[0] in maths:
                return Sdfmeta.metaoper(maths[tab[0]],[self.collapser(tab[1]),self.collapser(tab[2])]) #OR THIS
            elif tab[0] in pars:
                return self.collapser(tab[1], tab[0])
        elif isinstance(tab, str):
            if para in ('(','[','{'):
                trytab = [numify(i) for i in re.split('\s*[ ,]{0,1}\s*', tab )]
                if not str in map(type, trytab):
                    return Sdfmeta.construct( trytab )
            try:
                return slice(*[{True: lambda n: None, False: int}[x == ''](x) for x in (tab.split(':') + ['', '', ''])[:3]])
            except ValueError:
                #return Sdfmeta.construct( tab )
                raise ValueError('Your logic makes no sense. '+str(tab))
        elif isinstance(tab, Sdfmeta):
            return tab #OR THIS
        else:
            raise TypeError('Unknown logic')
            
    
    #def stripbutmeta(self, levemeta):
    def stripbutmeta(self, metaobj):
        '''
        Remove atoms except those defined in given liststructure
        '''
        self.numerize()
        #atoms = self.logicgetmeta(levemeta)
        atoms = metaobj
        newatoms = []
        for n in atoms._data:
            newatoms.append(self._atoms[n-1])
        self._atoms = newatoms
        self._bonds = []
        self._counts[0] = len(newatoms)
        self._counts[1] = 0

#end of Sdfmole

class Sdfmeta(object):
    '''
    Class including information in single metafield
    '''
    
    metaname = re.compile('\>(.*)\<(.+)\>') #Match gets metafield name
    ematch =   re.compile('[ ]{0,2}')
    
    def __init__(self, listofstrings=None):
        '''
        Initializes an empty metafield
        If listofstrings is given (a single metafield in .sdf-file) adds that data
        '''
        self._name = (None,None) #Metafield name
        self._datatype = None #int, float, str
        self._datastruct = None #list, dict, single
        self._data = None #The actual data
        self._delims = []
        self._dumb = True
        self._dumbcontent = []
        if listofstrings:
            self.initialize(listofstrings)
        
    def __getitem__(self, ind):
        '''
        get a slice from list or dict Sdfmeta
        '''
        if self._datastruct == OrDi:
            mylist = self._data.values()
        else:
            mylist = self._data
        if isinstance( ind, (slice, int) ):
            return mylist[ind]
        else:
            raise TypeError, "Invalid argument type."
            
    def __iter__(self):
        '''
        initializes an iterator to loop data in Sdfmeta
        returns the iterator
        '''
        if type(self._data) == OrDi:
            return iter(self._data.values())
        else:
            return iter(self._data)
            
    def __len__(self):
        #number of entries in Sdfmeta
        if self._data:
            return len(self._data)
        else:
            return 0
            
    #comparisons
    def __lt__(self, other):
        #less than
        return self._compare(other,operator.lt)
        
    def __le__(self, other):
        #less or equal
        return self._compare(other,operator.le)
        
    def __eq__(self, other):
        #equal
        return self._eqcompare(other,operator.eq)
        
    def __ge__(self, other):
        #greater or equal
        return self._compare(other,operator.ge)
        
    def __gt__(self, other):
        #greater than
        return self._compare(other,operator.gt)
        
    def __ne__(self, other):
        #not equal
        return not self._eqcompare(other,operator.eq)
    
    def _compare(self, other, oper):
        '''
        method to actually do given comparisons from example __gt__
        '''
        if self.isdumb():
            self.numerize()
        if type(other) != Sdfmeta:
            othermeta = Sdfmeta.construct( numify(other))
        else:
            othermeta = other
        if type(self._data) == OrDi:
            li1 = self._data.values()
        else:
            li1 = list(self._data)
        
        if type(othermeta._data) == OrDi:
            li2 = othermeta._data.values()
        else:
            li2 = list(othermeta._data)
        
        for item2 in li2:
            fulfils = True
            for item1 in li1:
                if not oper(item1, item2):
                    fulfils = False
                    break
            if fulfils:
                return True
        return False
    
    def _eqcompare(self, other, oper):
        '''
        same as _compare but for eq and ne
        '''
        if self.isdumb():
            self.numerize()
        if type(other) != Sdfmeta:
            othermeta = Sdfmeta.construct( numify(other))
        else:
            othermeta = other
        if type(self._data) == OrDi:
            li1 = self._data.values()
        else:
            li1 = list(self._data)
        
        if type(othermeta._data) == OrDi:
            li2 = othermeta._data.values()
        else:
            li2 = list(othermeta._data)
        
        for item2 in li2:
            for item1 in li1:
                if oper(item1, item2):
                    return True
        return False
        
    def initialize(self, listofstrings):
        '''
        parse the metadata from list of strings
        '''
        #get name of metafield
        if listofstrings[0][0] != '>' :
            self._name = ['  ','']
            fi = 0
        else:
            self._name = Sdfmeta.metaname.match(listofstrings[0]).groups()
            fi = 1
            if Sdfmeta.ematch.match(self._name[0]):
                self._name = ('  ', self._name[1])
        #Remove the last empty line and linechanges
        if not len(listofstrings[-1].strip())==0:
            #Raise error?
            return
        else: #             rstrip?
            self._data = [listofstrings[i].strip('\n ') for i in range(fi,len(listofstrings)-1) ]
    
    def numerize(self):
        '''
        Metas are initially stored in dumb-form (text). This method does the actual parsing/numerizing
        '''
        if self._dumb:
            
            #instead of commented section, merge lines with '' and do whattype
            #if no delimiter or ' ' in the end of line, add ' '
            mylines = self._data
            self._data = None
            
            newlines = []
            for line in mylines[:-1]:
                newlines.append(line)
                if not re.match('[ ,;\t]',line[-1]):
                    newlines.append(' ')
            newlines.append(mylines[-1])
            (data,dtype,delims) = Sdfmeta.whattype(''.join(newlines))
    
            #if string, it's special
            if dtype == str and type(data) != OrDi:
                self._datatype = str
                self._data = [line.strip('\n') for line in mylines]
                self._datastruct = list
                self._delims = ['' for line in mylines[:-1]]
                return
            #not just string
            
            self._data = data
            self._datatype = dtype
            if type(data) == list and len(data)==1:
                self._datastruct = 'single'
            else:
                self._datastruct = type(data)
                self._delims = delims
            
            self._dumb = False
        else:
            return
        
    @staticmethod
    def construct(data, **dictarg): #name, delims, literal
        '''
        return a new Sdfmeta from actual data (list, Ordered dict, string, int, float, ....
        '''
        new = Sdfmeta()
        new._dumb = False
        if 'name' in dictarg:
            name = dictarg['name']
            new.setname(name)
        else:
            name = ''
        if type(data) == str:
            
            if 'literal' in dictarg and dictarg['literal']:
                new._datastruct = 'single'
            else:
                (newdata, newtype, newdelims) = Sdfmeta.whattype(data)
                if newtype != str:
                    return Sdfmeta.construct(newdata, delims = newdelims, name = name)
            data = [data]
        elif type(data) in (OrDi, dict, list):
            if type(data)==dict:
                data = OrDi(data)
                new._datastruct = OrDi
            elif type(data)==list and len(data)==1:
                new._datastruct = 'single'
            else:
                new._datastruct = type(data)
        elif type(data) in (int, float):
            new._datastruct = 'single'
            data = [data]
        elif isinstance(data, Sdfmeta): #type(data) == Sdfmeta:
            out = copy.deepcopy(data)
            out.setname(name)
            if 'delims' in dictarg:
                out._delims = dictarg['delims']
            return out
        else:
            print type(data)
            raise TypeError('Datastructure type not list, dict, OrderedDict, str, int or float.')
        if not ('literal' in dictarg and dictarg['literal']):
            if type(data) == list:
                data = map(numify, data)
                
            elif type(data) == OrDi:
                for key in data:
                    newkey = numify(key)
                    data[newkey] = numify(data[key])
                    if key != newkey:
                        del(data[key])
                
        if type(data) == list:
            types = set(map(type, data))
        elif type(data) == OrDi:
            types = {type(data[key]) for key in data}
        else:
            types = {type(data)}
        
        new._data = data
        if not all({typ in (int, float, str) for typ in types}):
            raise TypeError('Data type not str, int or float.')
        if len(types)==1:
            new._datatype = iter(types).next()
        elif len(types)>1 and str in types:
            raise TypeError('Mixed datatypes. Strings and numbers.')
        else:
            new._datatype = float
            if type(new._data) == list:
                new._data = map(float, new._data)
            elif type(new._data) == OrDi:
                for key in new._data:
                    new._data[key] = float(new._data[key])
            else:
                raise TypeError('DaFaq?')
        if type(new._data) in (list, OrDi) and not 'delims' in dictarg:
            new._delims =  [ ', ' ]*(len(new._data)-1) 
        return new
            
            
    def __copy__(self):
        '''
        Shallow copy method
        '''
        new = Sdfmeta()
        new._name = tuple( self._name ) #new._name = copy.copy( self._name )
        new._datatype = self._datatype #new._datatype = copy.copy( self._datatype )
        new._datastruct = self._datastruct
        new._data = type(self._data)(( self._data )) #new._data = copy.copy( self._data )
        new._delims = list( self._delims ) #new._delims = copy.copy( self._delims )
        new._dumb = self._dumb 
        return new
        
    def __deepcopy__(self, memo):
        '''
        Deep copy method is the same as __copy__
        '''
        new = Sdfmeta()
        new._name = tuple( self._name ) #new._name = copy.copy( self._name )
        new._datatype = self._datatype #new._datatype = copy.copy( self._datatype )
        new._datastruct = self._datastruct
        new._data = type(self._data)(( self._data )) #new._data = copy.copy( self._data )
        new._delims = list( self._delims ) #new._delims = copy.copy( self._delims )
        new._dumb = self._dumb
        return new
    
    def getname(self):
        '''
        return name of Sdfmeta
        '''
        return self._name[-1]
    
    def isdumb(self):
        '''
        Has Sdfmeta been numerized yet
        '''
        return self._dumb
    
    def setname(self, newname):
        '''
        Change the name of Sdfmeta
        '''
        if type(newname) in (list, tuple):
            self._name = tuple(newname[:2])
        elif type(newname) == str:
            #self._name[1] = newname
            if not self._name[0]:
                #self._name[0] = '  '
                self._name = ('  ', newname)
            else:
                self._name = (self._name[0], newname)
        elif newname == None:
            self.setname('')
        else:
            raise TypeError('Name defined incorrectly! '+str(newname))
        
    def dtype(self):
        return self._datatype
        
    def extend(self, other):
        '''
        Used to merge more data to a Sdfmeta
        '''
        
        #If new Sdfmeta, do this
        if not self._datastruct:
            self._data =       copy.deepcopy( other._data )
            self._datastruct = copy.deepcopy( other._datastruct )
            self._datatype =   copy.deepcopy( other._datatype )
            self._delims =     copy.deepcopy( other._delims )
            self._dumb =       copy.deepcopy( other._dumb )
            self._dumbcontent =copy.deepcopy( other._dumbcontent )
            return
            
        floatflag = False
        if self._datatype != other._datatype:
            if self._datatype == str or other._datatype == str:
                raise TypeError('Mixed datatypes')
            else:
                floatflag = True
        
        elif self._datastruct == other._datastruct:
            if self._datastruct == list:
                self._data.extend(other._data)
            elif self._datastruct == OrDi:
                self._data.update(other._data)
            elif self._datastruct == 'single':
                self._datastruct = list
                self._data.extend(other._data)
                
        
        elif self._datastruct == OrDi or other._datastruct == OrDi:
            raise TypeError('Mixed datastructures')
        
        else:
            self._datastruct = list
            self._data.extend(other._data)
            
        try:
            self._delims.append(self._delims[-1])
            self._delims.extend(other._delims)
        except IndexError:
            if self._data:
                self._delims = [', ']*(len(self._data))
                self._delims.extend(other._delims)
                
        if floatflag:
            if self._datastruct == list:
                self._data = map(float, self._data)
            else:
                self._data = {k: float(v) for k, v in self._data.items()}
        
        
    def getmetastrings(self, length=float('inf')):
        def floattosting(flo):
            try1 = str(round(flo,4))
            if try1 == '0.0' and flo != 0:
                return '{:.4e}'.format(flo)
            else:
                return try1
        
        #stringers = {int:str, float:lambda x: str(round(x,4)), str:str}
        stringers = {int:str, float:floattosting, str:str}
        strifu = stringers[self._datatype]
        
        dictflag = self._datastruct == OrDi
        
        if dictflag:
            key = self._data.keys()[0]
            #tmp=[str(key)+':'+str(self._data[key])]
            tmp=[str(key)+':'+strifu(self._data[key])]
            del(key)
            itera = enumerate(self._data.keys()[1:])
            #itera = iter(self._data.keys()[1:])
        else:
            tmp=[strifu(self._data[0])]
            itera = enumerate(self._data[1:])
            #itera = iter(self._data[1:])
        l = len(tmp[-1])
        
        strings=[]
        for i, item in itera:
        #for item in itera:
            if dictflag:
                stuff = str(item)+':'+strifu(self._data[item])
            else:
                stuff = strifu(item)
            #add delimiter
            try:
                lde = len(self._delims[i])
            except IndexError:
                lde = 0
                print self._name
                print self._data
                print self._delims
            if l +lde > length:
                strings.append(''.join(tmp))
                tmp=[]
                l=0
            l += lde
            tmp.append(self._delims[i])
            #add data
            lda = len(stuff)
            if l +lda > length and lda < length: #FIX_extraline_strings_1
                strings.append(''.join(tmp))
                tmp=[]
                l=0
            l += lda
            tmp.append(stuff)
        if len(tmp)>0:
            strings.append(''.join(tmp))
        return strings
        
    def getmetastr(self):
        return self.getmetastrings()[0]
    
    def selftolistofstrings(self):
        '''return in .sdf-file format. No linechanges''' #TODO
        #Make list
        strings = []
        #Append name
        strings.append('>'+self._name[0]+'<'+self._name[1]+'>')
        
        if self.isdumb():
            #Add data
            strings.extend(self._data)
        else:
            #Do other things
            strings.extend(self.getmetastrings(74))
            
        #Apply the last blanck line
        strings.append('')
        return strings
        
    def selftostring(self):
        '''return in .sdf-file format. With linechanges'''
        return '\n'.join(self.selftolistofstrings)+'\n'
        
    def cleandelim(self, unify=False):
        '''Clean excess whitespaces from delimiters'''
        if len(self._delims) == 0:
            if self._data:
                self._delims = [', ']*(len(self._data)-1)
            return
        
        for i in range(len(self._delims)):
            self._delims[i] = self._delims[i].strip() + ' '
        if unify:
            ones = list(set(self._delims))
            nums = [0] * len(ones)
            for item in self._delims:
                nums[ones.index(item)] += 1
            newdelim = ones[nums.index(max(nums))]
            #for i in range(len(self._delims)):
                #self._delims[i] = newdelim
            self._delims = [newdelim] * (len(self._data)-1)
                
        
    @staticmethod
    def metaoper(oper, metas):
        '''
        method to make operations for Sdfmeta objects
        return resulting metavalue of operation.
        resulting Sdfmeta is nameless
        '''
        structs = [meta._datastruct for meta in metas]
        types   = [meta._datatype   for meta in metas]
        mathopers = (sum, sub, numpy.prod, max, min, avg, div, mypow)
        workmetas = [copy.copy(meta) for meta in metas] 
        if oper in mathopers:
            if str in types:
                #cannot calculate strings, doh
                return None
            
            if len(set(structs))>2:
                raise TypeError('Mixed datastructs')
            
            
            #if len(set(structs))==2 and 'single' in structs:
            #singles and (list or OrDi) in structs, extend singles to lists
            
            try:
                minlen = min( {len(meta._data) for meta in workmetas}-{1} ) # shortest list length, ignore singles
            except ValueError:
                minlen = 1
                
            try:
                ostru = iter(set(structs)-{'single'}).next() #the other structuretype
                singles = True
            except StopIteration:
                ostru = list
                singles = False
                
                
            
            if ostru == list:
                #make lists have the same length
                for i, meta in enumerate(workmetas): 
                    if meta._datastruct == 'single':
                        workmetas[i]._data = meta._data * minlen
                    else:
                        workmetas[i]._data = meta._data[:minlen]
                return Sdfmeta.construct( Sdfmeta.listoper(oper, workmetas, singles) ) #Nameless meta
            elif ostru == OrDi:
                keys = None
                keyorder = None
                for meta in workmetas:
                    if meta._datastruct == OrDi:
                        if not keys:
                            keys = set(meta._data)
                        else:
                            keys.intersection_update(set(meta._data))
                        if not keyorder:
                            keyorder = meta._data.keys()
                
                for i,meta in enumerate(workmetas):
                    if meta._datastruct == 'single':
                        newdict = OrDi()
                        for key in keyorder:
                            if key in keys:
                                newdict[key]=meta._data[0]
                        workmetas[i]=Sdfmeta.construct(newdict, name = meta.getname())
                #'''
                return Sdfmeta.construct(Sdfmeta.OrDioper(oper,workmetas))
        
        elif oper == metajoiner: #Therefore, it's joiner
            newmeta = metajoiner(metas)
            newmeta.cleandelim(True)
            return newmeta
        
        else:
            raise ValueError('Non-existent operator ')
        
    @staticmethod
    def listoper(oper, metas, singles=True):
        tab = []
        i=0
        while True:
            caltab = []
            for meta in metas:
                try:
                    caltab.append(meta._data[i])
                except IndexError:
                    if meta._datastruct == 'single' and singles: #if singles == False, chop the meta even if length == 1
                        caltab.append(meta._data[0])
                    else:
                        return tab
            tab.append(float(oper(caltab)))
            i+=1
        return tab
        
    @staticmethod
    def OrDioper(oper, metas):
        keys = set.intersection(*[set(meta._data) for meta in metas if meta._datastruct != 'single'])
        newkeys=[]
        if len(metas)>0:
            for key in metas[0]._data:
                if key in keys:
                    newkeys.append(key)
        keys = newkeys
        del(newkeys)
        dic = OrDi()
        for key in keys:
            caltab = []
            for meta in metas:
                if meta._datastruct == 'single':
                    caltab.append(meta._data[0])
                else:
                    caltab.append(meta._data[key])
            dic[key] = float(oper(caltab))
        return dic #OrDi(dic)
        
    @staticmethod
    def singleoper(oper, metas):
        if len({meta._datastruct for meta in metas} - {'single'})>0:
            return None
        else:
            return [oper([meta._data[0] for meta in metas])]
        
    def pickvalues(self, value, operator):
        if type(value)==Sdfmeta:
            pickvalue=value[0]
        else:
            pickvalue=numify(value)
        #MAY CREATE METAS WITH LENGT OF 0
        if type( self._data ) == OrDi:
            newdata = OrDi()
            for key in self._data:
                if operator(self._data[key], pickvalue):
                    newdata[key] = self._data[key]
        else:
            newdata = []
            for thing in self._data:
                if operator(thing, pickvalue):
                    newdata.append(thing)
        self._data = newdata
        self._delims = self._delims[:len(self._data)-1]
        
    def slicer(self, sliceorindex, paren):
        
        def itsslice( toget, slic ):
            #print 'SLICE!'
            if paren == '(':
                return Sdfmeta.construct(  OrDi( [ (i, toget[i]) for i in self._data.keys()[slic] ] ) )
            else:
                return Sdfmeta.construct(  toget[slic] )
        
        def itsmeta( toget, meta ):
            if type(meta._data) == OrDi:
                indexes = meta._data.values()
            else:
                indexes = meta._data
            if paren == '(':
                if type(self._data) != OrDi:
                    raise TypeError('"(" not applicaple for lists')
                #return Sdfmeta.construct(  OrDi([(i, toget[i]) for i in indexes]) )
                return Sdfmeta.construct(  OrDi([(i, toget[i]) for i in indexes if i in toget]) )
            le = len(toget)
            return Sdfmeta.construct(  [toget[i] for i in indexes if i < le ] )
        
        if paren == '{':
            toget = self._data.keys()
        elif paren == '[' and type(self._data) == OrDi:
            toget = self._data.values()
        else:
            toget = self._data
        
        if not type(toget) in (list, OrDi):
            raise TypeError('There is something wrong with our toget going to slicer') 
        slices = {Sdfmeta:itsmeta, slice:itsslice, NoneType: lambda toget, sliceorindex: None }
        return slices[type(sliceorindex)](toget, sliceorindex)
        
    def sortme(self, ascending=True, byValue = True):
        '''
        Sorts the cells to ascending or descending order, by value or key
        '''
        if byValue and type(self._data) == OrDi:
            self._data = OrDi( sorted( self._data.iteritems(), key=lambda x: x[1], reverse = not ascending ) )
        elif not byValue and type(self._data) == OrDi:
            self._data = OrDi( sorted( self._data.iteritems(), key=lambda x: x[0], reverse = not ascending ) )
        else:
            self._data = sorted( self._data, reverse = not ascending )
    
    @staticmethod
    def whattype(onestring):
        '''
        Tries to find out what type your string is
        Return (parsed,parsedtype,delim)
        '''
        if type(onestring) == str:
            if numitest(onestring):
                numoutof = numify(onestring)
            else:
                numoutof = onestring.strip()
        else:
            numoutof = numify(onestring)
        
        ty = type(numoutof)
        if ty in (int, float):
            #Number, return
            return ([numoutof],ty,[])
        #not just number, go on.
        del(numoutof,ty)
        
        
        def listtypetest(onestring, delimiter):
            #Is it a list delimited by given delimiter ?
            splitted = listsep(onestring,delimiter)
            #is dict?
            dicti = OrDi() #OrderedDict
            notdict = False
            for cell in splitted[0]:
                dicted = cell.split(':')
                if len(dicted)==2:
                    dicti[numify(dicted[0].strip())] = numify( dicted[1].strip() )
                    #good
                else:
                    #abort
                    del(dicti, dicted)
                    notdict = True
                    #print 'ding! not dict!'
                    break
            if notdict:
                if len(splitted[0])==1:
                    #therefore, no list of strings is generated
                    return None
                lister = [ numify( cell.strip() ) for cell in splitted[0] ]
                mytype = allsame(lister)
                if mytype: #all same
                    if mytype == str:
                        return None
                    elif mytype == float:
                        lister = map(float,lister)
                    return (lister,mytype,splitted[1])
                else:
                    pass
            else: #is dict
                mytype = allsame(dicti)
                if mytype:
                    return (dicti, mytype, splitted[1])
            return None
    
        chars = [';',',','\t',' ']
        for char in chars:
            news = listtypetest(onestring,char)
            #print news
            if news:
                return news
        return ([onestring],str,[])
    
#End of Sdfmeta

class Mol2File(object):
    def __init__(self, path = None):
        self._molecules = []
        if path:
            self.readfile(path)
    
    def __getitem__(self, ind):
        if isinstance( ind, slice ) or isinstance( ind, int ):
            return self._molecules[ind]
        else:
            raise TypeError, "Invalid argument type."
    
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
        
    def pickatomdata(self, column):
        newdic = OrDi()
        for line in self['ATOM']:
            things = re.split( '\s+', line.strip())
            if len(things)>1:
                if len(things[0])==0:
                    offset = 1
                else:
                    offset = 0
                newdic[numify(things[0+offset])] = things[column+offset]
        return newdic
    
    def injectatomdata(self, data, column, defaultValue=0.0, prec=4):
        
        if data._datastruct == list:
            offset = -1
        else:
            offset = 0
        
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
            #return []
        
        form = '{:.' + str(prec) +'f}'
        injectinfo = []
        for line in self['ATOM']:
            things = re.split( '(\s+)', line)
            mycol = None
            try:
                ind, mycol = findindex(things, (0, column))
            except ValueError:
                print things
                continue
            if mycol:
                try:
                    inject = data._data[numify(things[ind])+offset]
                except KeyError :
                    inject = defaultValue
                if math.isnan(inject) or math.isinf(inject):
                    inject = defaultValue
                stachar = sum( map(len, things[:mycol]) )
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

class Runner(object):
    order = ('input', 'tofield', 'toname', 'nametometa', 'removeconfname', 'removeconfmeta',  'cut', 
             'allcut', 'combine', 'allcombine', 'addcsv', 'getmol2', 'config', 
             'closestatom', 'closeratoms', 'changemeta', 'mergemeta', 
             'makenewmeta', 'sortmeta', 'stripbutmeta', 'extract', 
             'removemeta', 'pickmeta', 'putmol2', 'histogram', 
             'overwrite', 'output', 'getcsv', 'getatomcsv', 'metalist', 
             'counts', 'donotprint', 'split', 'makefolder')
    
    simplelist =    ('tofield', 'toname', 'removeconfname', 'removeconfmeta', 
                     'nametometa', 'removemeta', 'pickmeta', 'input', 
                     'histogram', 'proportion'
                     )
    
    simpleloops =   ('getmol2', 'closestatom', 'closeratoms', 'changemeta', 
                     'mergemeta', 'sortorder', 'sortmeta', 
                     'stripbutmeta', 'extract','makenewmeta','config',
                     'cut', 'allcut', 'combine', 'allcombine', 'addcsv', 
                     )
    
    #writers =       {'overwrite':lambda x : self.inpath, 'output': lambda x : x, 'stdout': lambda x: None} #default stdout
    writers =       ('overwrite', 'output', 'stdout') #default stdout
    writetypes =    {'getcsv':True,'getatomcsv':True,'metalist':True,'counts':True,'donotprint':True,'sdf':True,'split':False,'makefolder':False} #default 'sdf'
    #graphers =      ('histogram')
    
    
    #also verbose, propor
    
    #def __init__(self, inputfile, options=None):
    def __init__(self, options=dict()):
        self.options = options if options else dict()
        self.times = [time.time()]
        self.inpath = None
        #self.outpath = None
        self.wriarg = dict()
        self.writetype = 'sdf'
        self.writer = None
        self.sdf = Sdffile()
        '''
        self.path = inputfile
        if type(self.path) == str:
            self.sdf = Sdffile(self.path)
        else:
            self.sdf = Sdffile()
        self.times.append(time.time())
        '''
        
        #pop vs get
        #self.config = options.pop('config', None)
        self.verbose = options.pop('verbose', False)
        #self.propor = options.pop('proportion', False)
        self.setPropor(options.pop('proportion', False))
        
    def setPropor(self,parameter=None):
        if not parameter or parameter is True:
            self.propor = False
        else:
            self.propor = parameter.strip()
    
    def runOptions(self):
        
        #after config, do the rest anyway.
        for option in Runner.order:
            if option in self.options:
                self.times.append(time.time())
                self.funcselector(option, self.options[option])
        
        if self.verbose:
            print 'Run complete, it took {} seconds.'.format(self.times[-1]-self.times[0])
        
        #plt.show()
    
    def runConfig(self,confpath):
        
        def parsecon(oneline):
            
            params = [para.strip() for para in oneline.partition('#')[0].partition('::')]
            if params[0] != '' and params[1] == '::' and params[2] != '':
                #print params
                return (params[0], [numify(cell.strip()) for cell in params[2].split(';;')])
            elif params[0] != '':
                return (params[0], (True,))
            else:
                return None
        #use run-file
        #if self.config:
        with open(confpath, 'r') as confile:
            configures = (parsecon(confline)  for confline in confile )
            
            for option in configures:
                if option is None:
                    continue
                if option[0] in Runner.order: #self.options:
                    self.times.append(time.time())
                    self.funcselector(option[0], option[1])
        
    
    def taskLib(self,task,option,param=None,**kwargs):
        steps = kwargs.get('steps',1)
        timedif = lambda : self.times[-1]-self.times[-(1+steps)]
        #old       'task'           :(function,                 (parameters), (loopformat, loopdata), (initialformat, initialdata), (finalformat, finaldata))
        #new       'task'           :(function(*parameters),     (loopformat, loopdata), (initialformat, initialdata), (finalformat, finaldata))

        NoLamb = lambda : None
        tasks =   {
                    'tofield'        :lambda i : (self.sdf.addconfs, lambda : ((False,True),),      NoLamb,NoLamb,lambda : ('Conformation numbers added to metadata. It took {} seconds.', (timedif(),)))[i], 
                    'toname'         :lambda i : (self.sdf.addconfs, lambda : ((True,False),),      NoLamb,NoLamb,lambda : ('Conformation numbers added to names. It took {} seconds.', (timedif(),)))[i], 
                    'removeconfname' :lambda i : (self.sdf.remconfs, lambda : ((True,False),),      NoLamb,NoLamb,lambda : ('Conformation numbers removed from names. It took {} seconds.',(timedif(),)))[i], 
                    'removeconfmeta' :lambda i : (self.sdf.remconfs, lambda : ((True,False),),      NoLamb,NoLamb,lambda : ('Conformation numbers removed from metafield \'confnum\'. It took {} seconds.',(timedif(),)))[i], 
                    'nametometa'     :lambda i : (self.sdf.nametometa, lambda : (param,),           NoLamb,NoLamb,lambda : ('Name written to metafieldfield {}. It took {} seconds.',(param,timedif())))[i], 
                    'removemeta'     :lambda i : (self.sdf.removemeta, lambda : (param,False),      lambda : ('Metafieldfield {} Removed. It took {} seconds.',(param,timedif())),NoLamb,NoLamb)[i], 
                    'pickmeta'       :lambda i : (self.sdf.removemeta, lambda : (param,True),        NoLamb,NoLamb,NoLamb)[i], 
                    'getmol2'        :lambda i : (self.sdf.getMol2DataStr, lambda : (param,),       NoLamb,NoLamb,NoLamb)[i], 
                    'closestatom'    :lambda i : (self.sdf.closest, lambda : (param,),                  NoLamb,NoLamb,NoLamb)[i], 
                    'closeratoms'    :lambda i : (lambda x: self.sdf.closer(splitter(x)), lambda : (param,),     NoLamb,NoLamb,NoLamb)[i], 
                    'changemeta'     :lambda i : (lambda x: self.sdf.changemetaname([item.strip() for item in x.split('>')]), lambda : (param,), NoLamb, NoLamb, NoLamb)[i], 
                    'mergemeta'      :lambda i : (self.sdf.metamerger, lambda : (param,),           lambda : ('New metafield merged',None),NoLamb,NoLamb)[i], 
                    'sortorder'      :lambda i : (self.sdf.sortme, lambda : (param,),               lambda : ('Sort {} done.',(param,)),NoLamb,NoLamb)[i], 
                    'sortmeta'       :lambda i : (self.sdf.sortmetas, lambda : (param,),            lambda : ('Sort {} done.',(param,)),NoLamb,NoLamb)[i], 
                    'stripbutmeta'   :lambda i : (self.sdf.stripbutmeta, lambda : (param,),         lambda : ('All atoms, execpt for those in statement {} removed!',(param,)),NoLamb,NoLamb)[i], 
                    'extract'        :lambda i : (self.sdf.mollogicparse, lambda : (param,),        lambda : ('\n'.join(('After logical chop {}, sdf-file has {} molecules and {} conformations left.',self.sdf.propomes(self.propor))).strip() ,(param,len(self.sdf._dictomoles),len(self.sdf))), lambda : ('\n'.join(('Initially sdf-file has {} molecules and {} conformations.',self.sdf.propomes(self.propor))).strip(),(len(self.sdf._dictomoles),len(self.sdf))), lambda : ('Chopping complete, it took {} seconds.',(timedif(),)))[i], 
                    'makenewmeta'    :lambda i : (self.sdf.makenewmetastr, lambda : (param,),       lambda : ('New metafield {} made.',(param,)), NoLamb, lambda : ('Making new metafields done. It took {} seconds.',(timedif(),)))[i], 
                    'addcsv'         :lambda i : (self.sdf.addcsvmeta, lambda : (param, self.verbose,), lambda : ('Metadata from csv-file {} added. It took {} seconds.',(param,timedif,)),NoLamb,NoLamb)[i], 
                    'input'          :lambda i : (self.sdf.xreadself, lambda : (param,),            NoLamb, lambda : ('Starting to read file {}',(param,)), lambda : ('Reading file done. It took {} seconds.', (timedif(),)))[i], 
                    'config'         :lambda i : (self.runConfig, lambda : (param,),                lambda : ('Config-file {} done.',(param,)), lambda : ('Run config-files.',()), lambda : ('Running config-files done',()))[i], 
                    'histogram'      :lambda i : (self.sdf.histogramFromList, lambda : (param,),    NoLamb, lambda : ('Start plotting',()), lambda : ('Plotting done',()))[i], 
                    'cut'           :lambda i : (lambda x :self.sdf.sdflistremove(Sdffile(x.strip()),True), lambda : (param,), lambda : ('Removing all molecules (matching name) present in {} complete. It took {} seconds.', (param, timedif()), NoLamb, NoLamb))[i], 
                    'allcut'        :lambda i : (lambda x :self.sdf.sdflistremove(Sdffile(x.strip()),False), lambda : (param,), lambda : ('Removing all molecules (matching confnum) present in {} complete. It took {} seconds.', (param, timedif()), NoLamb, NoLamb))[i], 
                    'combine'       :lambda i : (lambda x :self.sdf.sdflistremove(Sdffile(x.strip()),(True,True)), lambda : (param,), lambda : ('Combining metadata from {} (matching name) complete. It took {} seconds.', (param, timedif()), NoLamb, NoLamb))[i], 
                    'allcombine'    :lambda i : (lambda x :self.sdf.sdflistremove(Sdffile(x.strip()),(True,False)), lambda : (param,), lambda : ('Combining metadata from {} (matching confnum) complete. It took {} seconds.', (param, timedif()), NoLamb, NoLamb))[i], 
                    'proportion'    :lambda i : (self.setPropor, lambda : (param,) ,NoLamb, NoLamb, lambda : ('Propor set to {}.',(param,)))[i], 
                   }
        selector = {'func':0,'loop':2,'initial':3,'final':4}
        try:
            if task == 'func':
                mytask=tasks.get(option)
                return mytask(0)(*mytask(1)())
            else:
                return tasks.get(option)(selector.get(task))()
        except KeyError:
            raise KeyError('Wrong task or option.')
    
    def funcselector(self,option,params,):
        
        def messenger(task, mypara=params,**kwargs):
            if task in ('final'):
                self.times.append(time.time())
            if self.verbose:
                mes = self.taskLib(task, option, mypara, **kwargs)
                if mes:
                    print mes[0].format(*mes[1])
        
        if option == 'input':
            messenger('initial')
            self.inpath = params
            self.taskLib('func', option, params)
            messenger('final')
            
        
        elif option in Runner.simpleloops or option in Runner.simplelist:
            if option in Runner.simplelist:
                params = (params,)
                
            messenger('initial')
            for oneparam in params:
                self.taskLib('func', option, oneparam)
                messenger('loop',oneparam)
                
            messenger('final', steps=len(params))
        elif option in Runner.writetypes:
            if Runner.writetypes[option]:
                self.writetype = option
            self.wriarg[option]=params
            
        elif option in Runner.writers:
            writers = {'overwrite':lambda x : self.inpath, 'output': lambda x : x, 'stdout': lambda x: None} #default stdout
            try:
                self.wriarg['path'] = writers.get(option, None)(params)
            except TypeError:
                pass
            self.writer = option
            self.sdf.writer(self.writetype,**self.wriarg)
    
    #look http://parezcoydigo.wordpress.com/2012/08/04/from-argparse-to-dictionary-in-python-2-7/
#End of Runner



    
def listtostring(tab, fill):
    return ''.join([('{:>'+str(fill)+'}').format(i) for i in tab])
    
'''
def nestfind(alist, subindex, tofind):
    return next((i for i, sublist in  enumerate(alist) if sublist[subindex]==tofind),-1)
'''
    
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
                return map(numify,stri)
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
            raise TypeError
    
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
            return numpy.asscalar((num[0])/numpy.prod(num[1:]))
    except ZeroDivisionError:
        return float('inf')
        
def mypow(num):
    if type(num) in (tuple, list):
        return numpy.asscalar(numpy.power(*num[:2]))
    else:
        return None

def metajoiner(metas, **params):
    '''
    Joins multiple metavalues into one
    params : name, delims
    '''
    newmeta = Sdfmeta()
    if 'name' in params:
        newmeta.setname(params['name'])
    else:
        newmeta.setname('')
    for meta in metas:
        newmeta.extend(meta)
    newmeta.cleandelim(True)
    return newmeta
    
def getcolumn(tab,head):
    ind = tab[0].index(head.strip())
    column = []
    for row in tab[1:]:
        column.append(row[ind])
    return column

def readcsv(path):
    f=open(path,'r')
    matrix = csvtomatrix(f.readlines())
    f.close()
    return matrix

def csvtomatrix(lines):
    chibuti=re.compile('[\s\n;]+')
    return [[numify(chibuti.sub('',cell)) for cell in line.split('\t')] for line in lines]

def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)


def splitter(stringofparams):
    return parentifier([cell.strip() for cell in stringofparams.split(',')],',')

def parentifier(original,separator):
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
    #if string.find('[')==0 and (string.rfind(']')-len(string))==-1:
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
    types = set(map(type, stuff))
    if len(types) == 1:
        return types.pop()
    elif not str in types:
        return float
    else:
        return None

    
def leveler(string):
    pair = {'"':'"',"'":"'",'(':')','[':']','{':'}'}
    pars = [(item.group(), item.start()) for item in parre.finditer(string)]
    
    tab = []
    
    if len(pars)<2:
        return [string]
    
    curpar = None
    backpar = None
    
    if pars[0][1]>0:
        stastring = string[:pars[0][1]]
        #if len(stastring)>0:
        tab.append(stastring)
        
    for item in pars:
        
        if not curpar:
            curpar = item
            level = 1
            if backpar:
                addy = string[backpar[1]+1:curpar[1]]
                if len(addy) > 0:
                    tab.append(addy)
                backpar = None
        elif item[0] == pair[curpar[0]]:
            level -= 1
            if level == 0:
                tab.append((curpar[0],leveler(string[curpar[1]+1:item[1]]))) #curpar[1],item[1],
                curpar = None
                backpar = item
        elif item[0] == curpar[0]:
            level += 1
    endstring = string[pars[-1][1]+1:]
    if len(endstring)>0:
        tab.append(endstring)
    return tab

def dumb_levopemap(tab, par=None, length=None):
    '''
    Accepts lists made by leveler method. Maps this list for mathematical operators and metafield names. Doesn't understand local metafields. a dumb version.
    '''
    if isinstance(tab, tuple):
        #if par and par in ('"',"'"):
        if tab[0] in ('"',"'"):
            #return Sdfmeta.construct(levjoin(tab[1]))
            return (tab[0], levjoin(tab[1]))
        else:
            return (tab[0], dumb_levopemap(tab[1], tab[0]))
    elif isinstance(tab, list):
        for j, thing in enumerate(tab):
            if not isinstance(thing, str):
                continue
            thing = thing.strip()
            #new
            sear = re.search('[+-]', thing)
            if sear:
                i = sear.start()
                ie = sear.end()
            if not sear:
                sear = re.search('(^|[^\*])[\*/]([^\*]|$)', thing)
                if sear:
                    helpse = re.search('[\*/]', thing[sear.start():sear.end()])
                    i = sear.start()+helpse.start()
                    ie = sear.start()+helpse.end()
                    del(helpse)
            if not sear:
                sear = re.search('[\*]{2}', thing)
                if sear:
                    i = sear.start()
                    ie = sear.end()
            if sear:
                tab1=[]
                tab1.extend(tab[:j])
                if i > 0:
                    tab1.append(thing[:i])
                if ie < len(thing):
                    tab2 = [thing[ie:]]
                else:
                    tab2 = []
                tab2.extend(tab[j+1:])
                return (thing[i:ie],dumb_levopemap(tab1),dumb_levopemap(tab2))
            #new ends
        return [dumb_levopemap(item, par,len(tab)) for item in tab]
    else:
        return tab.strip()

def levels(string):
    return dumb_levopemap(leveler(string))

def tabjoin(taber):
    pair = {'(':')','[':']','{':'}'}
    string = []
    for item in taber:
        if isinstance(item, str):
            string.append(item)
        elif isinstance(item, (list,tuple)):
            string.extend([item[0],tabjoin(item[1]),pair[item[0]]])
    return ''.join(string)
    
def levjoin(taber):
    #print taber
    pair = {'(':')','[':']','{':'}',"'":"'",'"':'"'}
    string = []
    for item in taber:
        if isinstance(item, str):
            string.append(item)
        elif isinstance(item, list):
            string.append(levjoin(item))
        elif isinstance(item, tuple):
            if item[0] in pair:
                string.extend([ item[0], levjoin(item[1]), pair[item[0]] ])
            else:
                string.extend([levjoin(item[1][0]), item[0], levjoin(item[1][1]) ])
    return ''.join(string)

    
    

if __name__ == "__main__":
    
    #main should only collect arguments...
    
    arger = argparse.ArgumentParser(description= 'Some bad-ass manipulation of SDF-files. Notice that documentation is not up to date.' )
    
    arger.add_argument("input", metavar = 'path', nargs='+', type = str,               help="Specify the input file")
    choicewrite = arger.add_mutually_exclusive_group()
    choicewrite.add_argument("-out", "--output", type = str, default=None,  help = "Specify output file.")
    choicewrite.add_argument("-o", "--overwrite", action='store_true',      help = "overwrite. you don't need to specify output file")
    
    arger.add_argument("-con","--config", metavar = 'path', nargs='+', type = str, help="Specify the config file")
    
    arger.add_argument("-cf", "--tofield", action = "store_true",           help = "add conformation number to metadata from name. if number in name doesn't exist, make a new one.")
    
    arger.add_argument("-cn", "--toname",  action = "store_true",           help = "add conformation number to name from metadata. if number in metadata doesn't exist, make a new one.")
    arger.add_argument("-mtn", "--metatoname",  type = str,                 help = "Change the name of molecule to the data in given metafield")
    arger.add_argument("-ntm", "--nametometa",  type = str,                 help = "Copy the name of molecule into given metafield")
    arger.add_argument("-rcn", "--removeconfname",  action="store_true",    help = "remove conformation number from name.")
    arger.add_argument("-rcm", "--removeconfmeta",  action="store_true",    help = "remove conformation number from metafield 'confnum'.")
    
    choicecombi = arger.add_mutually_exclusive_group()
    choicecombi.add_argument("-co", "--combine", type = str, nargs='+',     help = "Combine metadata from specified file to the data of original file. Confromation numbers must match.")
    choicecombi.add_argument("-aco", "--allcombine", type = str, nargs='+', help = "Combine metadata from specified file to the data of original file. Names must match.")
    
    choicecut = arger.add_mutually_exclusive_group()
    choicecut.add_argument("-cu", "--cut", type = str, nargs='+',           help = "Remove molecules in specified file from original file. Confromations must match.")
    choicecut.add_argument("-acu", "--allcut", type = str, nargs='+',       help = "Remove molecules in specified file from original file. Names must match. Not tested")
    
    arger.add_argument("-csv", "--addcsv", type = str, nargs='+',           help = "Add metadata from csv-file. File must have a 1-line header, it gives names to metafields. Names of molecules must be leftmost. If name includes confnumber, meta is only added molecules with same confnumber.")
    arger.add_argument("-ex", "--extract", type = str, nargs='+',           help = "Pick or remove molecules from file by metafield info. Either with logical comparison or fraction of molecules with same name. Closest_atoms{:5}==soms, 2.5>Closest_atoms(soms)[], Closest_atoms[:3]<5.5, ID='benzene'. Takes multiple statements separated with | ")
    arger.add_argument("-pro", "--proportion", type = str,                  help = "Takes one exctract-like metastatement and prints proportion of molecules and conformations fulfilling it after every chop, if you are on verbose mode.")
    
    choiceremo = arger.add_mutually_exclusive_group()
    choiceremo.add_argument("-rm", "--removemeta", type = str,              help = "Remove metadata from molecules. Takes multiple values, separaterd by comma(,) or semicolon(;). If first is '?', means 'all but'")
    choiceremo.add_argument("-pm", "--pickmeta", type = str,                help = "Remove all nonspecified metadata from molecules. Takes multiple values, separaterd by comma(,) or semicolon(;). If first is '?', means 'all but'")
    
    arger.add_argument("-s", "--split", type = int,                         help = "Split the file into even pieces. Positive means number of files while negative means number of molecules per file. 0 doesn't apply.")
    arger.add_argument("-mf", "--makefolder", action = "store_true",        help = "Put outputfile(s) into folder(s).")
    
    outputtype = arger.add_mutually_exclusive_group()
    outputtype.add_argument("-gc", "--getcsv", type = str ,                 help = "Writes a .csv-file istead of .sdf-file. Specify which fields you'll need, separated by ','.")
    outputtype.add_argument("-gac", "--getatomcsv", type = str ,            help = "Writes a .csv-file istead of .sdf-file. Specify which fields you'll need, separated by ','. Writes dictionaries to separate lines")
    outputtype.add_argument("-ml", "--metalist", action = "store_true",     help = "Writes a list of metafields.")
    outputtype.add_argument("-nm", "--counts", type = int,                  help = "Number of different molecules and different conformations. 0=just sums, 1=by molecule name, 2=both")
    outputtype.add_argument("-dnp", "--donotprint", action = "store_true",  help = "No output")
    
    arger.add_argument("-ca", "--closestatom", type = str, nargs='+',       help = "Calculates the closest atoms (distances by atom number) to given point. Adds 'Closest_atoms' metafield with optional prefix. Needs either coordinates separated by ',' or or single atom number")
    arger.add_argument("-cla", "--closeratoms", type = str, nargs='+',      help = "Calculates number of atoms closer to the given point, than the ones given adds metafields 'Closest_atom_from_{meta}' and 'Closer_atoms_than_{meta}'. Needs point and metafield name separated by ',', point first. Takes multiple parameters separated by '|'")
    arger.add_argument("-v", "--verbose", action = "store_true" ,           help = "More info on your run.")
    arger.add_argument("-mm", "--mergemeta", type = str, nargs='+',         help = "Makes a new metafield based on old ones. newmeta=sum(meta1,meta2). operator are sum, max, min, avg, prod, div, power and join. Takes�multiple arguments separated by |")
    arger.add_argument("-mnm", "--makenewmeta", type = str, nargs='+',      help = "Makes a new metafield based on logical statement and value picking inside the metafield. newmeta = meta1 + meta2 < 50. Takes�multiple arguments separated by |")
    arger.add_argument("-cm", "--changemeta", type = str, nargs='+',        help = "Changes names of metafields. [olname1>newname1|oldname2>newname2]. ")
    arger.add_argument("-sm", "--sortmeta", type = str, nargs='+',          help = "Sorts cells of a metafield in ascending [<metaname] or descending [>metaname] order. Additional '+' as the last character sorts by key in dictionary type metas. Takes multiple values separated by |")
    arger.add_argument("-so", "--sortorder", type = str, nargs='+',         help = "Sorts molecules of a file in order of metafield. <MolecularWeight|>Id Sorts molecules first ascending by weight, then descenting by name")
    arger.add_argument("-hg", "--histogram", type = str, nargs='+',         help = "Plots a 1D or 2D histogram, multiple plots with '|'. 'Xname,Yname,Xtitle=x-akseli,Ytitle=y-akseli,bins=[30,30]'")
    
    arger.add_argument("-gm2", "--getmol2", type = str, nargs='+',          help = "Reads atom block column data from mol2-file and adds it to sdf-file as metadata. pathto.mol2, column, metaname.")#meta column path
    arger.add_argument("-pm2", "--putmol2", type = str, nargs='+',          help = "Injects meta-data from sdf-file and adds it to mol2-file as atom block column data. input.mol2,output.mol2, column, metaname, default, precision.")#metaname, column, path, defaultValue, precision, outpath
    
    arger.add_argument("-sbm", "--stripbutmeta", type = str, nargs='+',     help = "Removes all atoms from molecules, except for those in given logical statement. Takes multiple parameters separated by '|'")
    
    args = arger.parse_args()
    manyfiles = args.input
    
    
    def main(inputfile, options):
        #variab = vars(args)
        options['input'] = inputfile
        onerun = Runner(options)
        onerun.runOptions()
        
        if onerun.writer==None:
            onerun.sdf.writer(onerun.writetype,**onerun.wriarg)
        
    
    for onefile in manyfiles:
        options = vars(args)
        for key in options.keys():
            if not options[key] and type(options[key])!=int:
                del(options[key])
        main(onefile, options)
    
    plt.show()

