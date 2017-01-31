# -*- coding: utf-8 -*-

import sys
import re
import math
import numpy
import copy
import operator
from collections import OrderedDict as OrDi
import warnings
from six import string_types

try:
    from future.utils import lmap
except ImportError:
    lmap = map

try:
    import functions, mol2
except ImportError:
    from sdfconf import functions, mol2
    
try:
    iterange = xrange
except NameError:
    iterange = range

class Sdffile(object):
    '''
    Class that represents one .sdf-file. Includes all molecules in nested dictionaries, first by molecule name, then by conformation number.
    Order of molecules (in the .sdf-file) is represented by a list containing cells with molecule name and conformation number.
    '''
    
    molchop = re.compile('^\${4}') #Separates molecules in file
    
    
    def __init__(self, path=None, **kwargs ):
        '''
        Initianilizes an empty file. If path is specified, add data from it.
        '''
        self._dictomoles = dict() # {'aspirin':{1:aspi1, 2:aspi2, ...}, 'bentzene':{1:benz1, 2:benz2, ...}, ...}
        self._orderlist = list()
        #self._ignores = kwargs.get('ignores', ['H'])
        self.setIgnores(kwargs.get('ignores', ['H']))
        if path is not None:
            self.xreadself(path)
    
    def __copy__(self):
        '''
        Shallow copy function used by copy.copy()
        '''
        new = Sdffile()
        new._dictomoles = dict()
        for key in self._dictomoles:
            new._dictomoles[key] = copy.copy( self._dictomoles[key] )
        new._orderlist = copy.copy(self._orderlist)
        new._ignores = copy.copy(self._ignores)
        return new
    
    def __deepcopy__(self,memo):
        '''
        Deep copy function used by copy.deepcopy()
        '''
        new = Sdffile()
        new._dictomoles = copy.deepcopy(self._dictomoles,memo)
        new._orderlist = copy.deepcopy(self._orderlist,memo)
        new._ignores = copy.copy(self._ignores)
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
            raise TypeError("Invalid argument type.")
        
    def setIgnores(self, ignores=['H']):
        self._ignores = list(ignores)
        #print self._ignores
        for mol in self:
            mol.setIgnores(self._ignores)
        
        
    @staticmethod
    def sdfseparator(strings):
        '''
        Separates a list of strings (as in .sdf-file) into list of lists of strings.
        lists of strings represent single conformations
        return list of lists of strings
        '''
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
        _lines = []
        for line in xfile:
            if Sdffile.molchop.match(line):
                #_moles.append(_lines)
                yield _lines
                _lines = []
            else:
                _lines.append(line)
    
    def add(self, stringsofone):
        '''
        Adds a molecule to datastructure. 
        stringsofone is a list of strings representing a single conformation
        '''
        new = Sdfmole(stringsofone, ignores = self._ignores)
        name = new.getname()
        
        if not name in self._dictomoles:
            self._dictomoles[name]=dict()
        n = self.tempconfn(new, name)
            
        i=0
        while n in self._dictomoles[name]:
            i+=1
            n=str( max(lmap(int,self._dictomoles[name].keys())) + i)
            
        self._dictomoles[name][n] = new
        self._orderlist.append([name,n])
        
            
    def remove(self, name, confn):
        #Remove a molecule by name and conformation number.
        try:
            del(self._dictomoles[name][str(confn)])
            self._orderlist.remove([name,confn])
        except ValueError:
            warnings.warn('{}{{[{}]}} not in list.'.format(name, confn))
        
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
                
    def listremove(self, path, sameconf=True):
        '''Frontend to sdflistremove and csvlistremove'''
        other = Sdffile(path)
        if len(other)==0:
            del(other)
            liomo = functions.csvtomatrix(path,'[,;\t ]')
            self.csvlistremove(liomo, sameconf)
        else:
            self.sdflistremove(other, sameconf)
            
    def csvlistremove(self, others, sameconf=True):
        '''
        Cut operator.
        Remove conformations from current file that are present in the other csvfile.
        ''' 
        for moli in others:
            oname = Sdfmole.confchop.sub('',moli[0]).strip('\"\'')
            if sameconf:
                try:
                    confn = Sdfmole.confchop.search( moli[0] ).group().strip('{[]}\"\'')
                    self.remove(oname,confn)
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
    
    def sdflistremove(self, other, sameconf=True):
        '''
        Cut operator.
        Remove conformations from current file that are present in the other sdffile.
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
        metas = self.getmollogic(metalogic)
        for molinfo in self._orderlist:
            self._dictomoles[molinfo[0]][molinfo[1]].stripbutmeta(metas[molinfo[0]][molinfo[1]])
    
    #conformation numbers
    
    def tempconfn(self, mol, molname=None):
        '''
        Returns a temporary conformation number for given molecule.
        '''
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
        '''
        Returns a unique conformation number for given molecule.
        '''
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
        for i, molord in enumerate(self._orderlist):
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
    
    def makenewmetastr(self, string): #no more accepts new = old < 5; now you must write new = old(<5); you can also new = old2(old1(<5){})(>3)
        '''Frontend for makenewmeta'''
        things = Sdfmeta.compsplit(string,comps=('=',))
        if len(things) == 3 and things[1] == '=': 
            self.makenewmeta(*(things[:1]+things[2:]))
        else:
            raise functions.InputException('Invalid new meta definition {}'.format(str(things)))
            
    
    def makenewmeta(self, name, metastatement): #, logicchar = None, value = None):
        '''
        make a new metafield by metastatement and pick only those fullfilling given logical statement
        '''
        newmetas = self.getmollogic(metastatement)
        count = 0
        for molname in newmetas:
            for confn in newmetas[molname]:
                if len(newmetas[molname][confn])>0:
                    count += 1
                    self._dictomoles[molname][confn].addmeta(name, newmetas[molname][confn], overwrite=True)
        if count<len(self):
            warnings.warn('Not all new metas were generated',UserWarning)
        
    def nametometa(self, meta):
        '''Adds a metafield holding the molecule name'''
        for mol in self:
            mol.addmeta(meta, mol.getname(), literal=True)
            
    def changemetaname(self, oldname, newname):
        '''Change the name of a metafield'''
        for mol in self:
            mol.changemetaname(oldname, newname)
            
                
    def dictmaint(self): 
        '''Maintenance of Sdffile dictionaries'''
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
        '''Test wether molecules are "the same" as described in bolist.'''
        return (bolist[0] <= samelist[0]) and (bolist[1] <= samelist[1])
        
    def csvStrSplit(self,stringofmetas, sepa = ',;'):
        '''
        Separate string of metanames by characters in listsepa to a list. If list includes '?', all other metas except those in list are returned.
        '''
        #listofmetas = [met.strip() for met in re.split('\s*{}\s*'.format('|'.join(lmap(lambda x : '('+x+')', listsepa ))) , stringofmetas)]
        listofmetas = [met.strip() for met in re.split('\s*[{}]\s*'.format(sepa) , stringofmetas)]
        return self.getMetas(listofmetas)
        
    def getMetas(self, listofmetas):
        '''
        Returns list of existing metas by argument list. If plain list, 
        removes those that doesn't exist. If includes '?', returns all 
        except for those in that list.
        '''
        metalist = self.listmetas()
        if '?' in listofmetas:
            for meta in listofmetas:
                if meta == '?':
                    continue
                metalist.remove(meta)
            listofmetas = metalist
        else:
            for meta in list(listofmetas):
                if meta not in metalist:
                    listofmetas.remove(meta)
        return listofmetas
    
    
    def makeCsvStr(self, stringofmetas, **kwargs):
        '''
        Generate a csv-file from a string containing list of metanames.
        Setting atomic True generates atomic csv.
        '''
        
        liargs, diargs = kwarger(stringofmetas, [',',';'])
        kwargs = dict(kwargs)
        for key in diargs:
            kwargs[key] = diargs[key]
        
        atomic    = kwargs.get('atomic', False)
        separator = kwargs.get('separator','\t')
        
        #metalist = self.csvStrSplit(stringofmetas)
        metalist = self.getMetas(liargs)
        
        
        if not atomic:
            return self.makeCsv(metalist, separator)
        else:
            return self.makeAtomicCsv(metalist, separator)
        
    '''
    def makeCsvStr(self, stringofmetas, **kwargs):
        \'''
        Generate a csv-file from a string containing list of metanames.
        Setting atomic True generates atomic csv.
        \'''
        
        separator = kwargs.get('separator','\t')
        atomic    = kwargs.get('atomic', False)
        
        metalist = self.csvStrSplit(stringofmetas)
        
        if not atomic:
            return self.makeCsv(metalist, separator)
        else:
            return self.makeAtomicCsv(metalist, separator)
    '''
        
    def makeCsv(self,listofmeta,separator='\t'):
        '''
        Make a csv-list containing all molecules and given metafields as columns. '?' gives all fields
        '''
        ###TESTI!!!
        if len(listofmeta) == 0:
            listofmeta = self.listmetas()
        csv = [separator.join(listofmeta)]
        for info in self._orderlist:
            mol = self._dictomoles[info[0]][info[1]]
            line = []#[mol._name]
            for meta in listofmeta:
                #if meta in mol._meta:
                if mol.hasmeta(meta):
                    memeta = mol.getmeta(meta) #.getmetastr()
                    #if memeta._datastruct != 'single' or  memeta.dtype() == 'str':
                    #print(len(memeta), memeta.isStr())
                    #print('{} : struct={}, type={}, length={}, str?={}'.format(str(memeta), memeta._datastruct,  memeta.dtype(),len(memeta),memeta.isStr()))
                    #if memeta._datastruct != 'single' or  memeta.isStr():
                    if len(memeta)>1 or memeta.isStr() or memeta._datastruct == OrDi:
                        line.append('"'+memeta.getmetastr()+'"')
                    else:
                        line.append(mol.getmeta(meta).getmetastr())
                else:
                    line.append('""')
            csv.append(separator.join(line))
        return csv
    
    def makeAtomicCsv(self,listofmeta,separator='\t'):
        '''
        Make a csv-list containing all molecules and given metafields as 
        columns. '?' gives all fields
        '''
        #listofmeta = [met.strip() for met in re.split('\s*[,;]\s*',stringofmetas)]
        
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
                        #newline.append("NA")
                        newline.append('""')
                    elif mol.getmeta(meta)._datastruct == OrDi:
                        #newline.append(str(mol.getmeta(meta)._data.get(key,"NA")))
                        newline.append(str(mol.getmeta(meta)._data.get(key,'""')))
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
            counts.append((mol,len(self._dictomoles[mol])))
        return counts
    
    def closestStr(self, string):
        #String frontend for closest
        argus = functions.splitter(string)
        kwargs = {}
        if len(argus)>1:
            kwargs['name']=argus[1]
        for arg in argus[2:]:
            #key, op, value = Sdfmeta.compsplit(arg, ('=',))
            key, value = Sdfmeta.compsplit(arg, ('=',))[::2]
            kwargs[key] = value
        self.closest(argus[0],**kwargs)
    
    def closest(self, point, **varargdict): #name, intrests, num
        '''
        Calculates distaces from atoms to point of interest. Distaces by 
        atomnumbers are added to a metafield. point is a coordinate in 3d 
        or a name of metafield. If later, first atom number in meta 
        is picked and it's location is used. 
        '''
        if not 'name' in varargdict:
            name = 'Closest_atoms'
        else:
            name = varargdict['name']
        
        if 'interests' in varargdict:
            switch = {
                      list: lambda : (varargdict['interests'],True),
                      int:  lambda : ([varargdict['interests']],True),
                      str:  lambda : (self.getmollogic(varargdict['interests']), False) }
            intrestings, listbool = switch.get( type(varargdict['interests']), lambda : ([], True) )()
        
        myignores = varargdict.get('ignores', self._ignores)
        
        precalc = dict()
        for molname, conf in self._orderlist:
            if not molname in precalc:
                precalc[molname]=dict()
            mol = self._dictomoles[molname][conf]
            #is coord1 meta or not? If it was a metafield with atomnumber, it would give that one as closest, so no.
            alist = mol.dists(point, precalc=precalc[molname], dictomo=self._dictomoles[molname], ignores = myignores)
            alist=sorted(alist, key=lambda item: item[0])
            moles=[]
            if 'interests' in varargdict:
                mylist = intrestings if listbool else intrestings[molname][conf]._data
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
            mol.addmeta(name,od)
            
    
    def closer(self, point, meta):
        '''
        Calculates how many atoms are closer to the given coordinate
        than the atom in given metafield.
        '''
        metas = self.getmollogic(meta)
        precalc = dict()
        for molname, conf in self._orderlist:
            if not molname in precalc:
                precalc[molname] = dict()
            mol = self._dictomoles[molname][conf]
            alist = sorted(mol.dists(point, dictomo=self._dictomoles[molname], precalc=precalc[molname] ), key  = lambda x: x[0] , ignores=self._ignores )
            for i, atom in enumerate(alist):
                if atom[1] in metas[molname][conf]._data:
                    mol.addmeta('Closer_atoms_than_'+meta.strip(), i)
                    mol.addmeta('Closest_atom_from_'+meta.strip(), atom[1])
                    break
    
    def escapeStr(self, string, inside=False, **kwargs):
        '''
        Generate escapenum field for all molecules. (Number of atoms in self not in range of other molecule.)
        '''
        try:
            from findable import Findable
        except ImportError:
            from sdfconf.findable import Findable
        
        try:
            spli = re.split('\s*,\s*', string)
            otherpath, molnum, maxRange, name = spli[:4]
            manx = spli[4] if len(spli)>4 else 0
        except ValueError:
            return
        if not inside:
            name = name+'_escapenum' if len(name)>0 else 'escapenum'
            manx = -manx
        else:
            name = name+'_insidenum' if len(name)>0 else 'insidenum'
        
        omol = {'sdf':Sdffile, 'mol2':mol2.Mol2File}.get(otherpath.rpartition('.')[2],Sdffile)(otherpath)[int(molnum)]
        
        myignores = kwargs.get('ignores',self._ignores)
        
        matcher = Findable(omol, ignores=myignores)
        
        
        
        for mmol in self:
            numberIn, numberOut = mmol.calcEscapeNumberOrder(matcher,float(maxRange) , maxn=manx, anums = True, ignores=myignores) #also ignores
            if not inside and numberIn is not None:
                if isinstance(numberOut, (list,tuple)) and len(numberOut)==0:
                    continue
                mmol.addmeta(name, numberOut)
            elif inside and numberIn is not None:
                if isinstance(numberIn, (list,tuple)) and len(numberIn)==0:
                    continue
                mmol.addmeta(name, numberIn)
            else:
                warnings.warn('Metafield {} not created for all molecules.'.format(name))
    
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
        '''Test if string is a 3d coordinate. If it is, return it, otherwise return None'''
        cords = lmap(functions.numify,re.split('\s*',canditate))
        if len(cords)==3:
            if  lmap(type,cords)==[float,float,float]:
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
            for mole in Sdffile.xsdfseparator(f):
                self.add(mole)
    
    def selftostring(self, output, **kwargs):
        if output=='getcsv':
            return '\n'.join(self.makeCsvStr(kwargs['getcsv'], atomic = False ))+'\n'
        elif output=='getatomcsv':
            return '\n'.join(self.makeCsvStr(kwargs['getatomcsv'], atomic = True ))+'\n'
        elif output=='metalist':
            return '\n'.join(self.listmetas())+'\n'
        elif output=='counts':
            towrite=[]
            if kwargs['counts']==0 or kwargs['counts']==2:
                towrite = ['Total {}\t{}'.format(len(self._dictomoles),len(self))]
            if kwargs['counts']==1 or kwargs['counts']==2:
                towrite.extend(('\t'.join(lmap(str, count)) for count in self.counts()))
            return '\n'.join(towrite)+'\n'
        elif output=='sdf':
            return str(self)
        elif output == 'donotprint':
            return None
        else:
            return None
    
    def numerizeAll(self):
        '''
        Numerize all possible dumb data in the file.
        '''
        for mol in self:
            mol.numerize()
            for key in mol._meta:
                mol._meta[key].numerize()
    
    def histogrammer(self, Xname, Yname=None, **kwargs):
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        from matplotlib  import cm
        self.plt=plt
        if 'ex' in kwargs:
            sdf = copy.copy(self)
            sdf.mollogicparse(kwargs['ex'])
            del(kwargs['ex'])
        else:
            sdf = self
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
            if 0 in lmap(len, metas): #skips 0 len metas
                continue
            try:
                minlen = min( {len(meta) for meta in metas}-{1} ) # shortest list length, ignore singles
            except ValueError:
                minlen = 1
            try:
                structs = set([meta._datastruct for meta in metas])-{'single'}
            except StopIteration:
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
                            keyorder = list(meta._data.keys())
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
                        print(mol._name)
                        print(meta._data)
                        print(keyorder)
                        raise KeyError(str(key))
            
        self.plt.figure()
        if not Yname:
            larg=[datas[0]]
            if 'bins' in kwargs:
                larg.append(kwargs['bins'])
                del(kwargs['bins'])
            self.plt.hist(*larg,**kwargs) #kwargs?
        else:
            #cmap = mpl.cm.jet
            cmap = cm.jet
            #cmap = cm.get('jet')
            X=self.plt.hist2d(datas[0],datas[1],**kwargs)[0]
            
            ticks=list(numpy.arange(numpy.max(X)+1))
            skip = int(round(float(len(ticks))/20))
            if skip>1:
                ticks = numpy.array(ticks[:1] + ticks[skip:-skip:skip] + ticks[-1:])
            norm = mpl.colors.BoundaryNorm(numpy.arange(-0.5,numpy.max(X)+1.5), cmap.N)
            cbar = self.plt.colorbar()
            cbar = mpl.colorbar.ColorbarBase(cbar.ax , cmap=cmap,norm=norm,ticks=ticks,spacing='proportional')
        
        if 'Xtitle' in newargs:
            self.plt.xlabel(newargs['Xtitle'])
        else:
            self.plt.xlabel(Xname)
        if 'Ytitle' in newargs:
            self.plt.ylabel(newargs['Ytitle'])
        elif Yname:
            self.plt.ylabel(Yname)
        if 'title' in newargs:
            self.plt.title(newargs['title'])
    
    def show(self):
        self.plt.show()
    
    def histogramFromList(self,plots):
        showflag=True
        for plot in plots:
            path=None
            params = functions.splitter(plot)
            larg = []
            darg = dict()
            for para in params:
                things = Sdfmeta.compsplit(para, ('=',))
                if len(things)==3:
                    darg[things[0]]=things[2]
                else:
                    larg.append(para)
                
            for key in darg:
                darg[key]=functions.lister(darg[key])
                darg[key]=functions.numify(darg[key])
            if 'save' in darg:
                path=darg['save']
                del(darg['save'])
            if 'noplot' in larg:
                showflag=False
                larg.remove('noplot')
            self.histogrammer(*larg,**darg)
            if path:
                self.plt.savefig(path, bbox_inches='tight')
        if showflag:
            self.show()
    
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
        porconf = float(len(trues))/( len(trues)+len(falses) )
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
        for info in drops:
            del(self._dictomoles[info[0]][info[1]])
        neworder=[]
        for info in self._orderlist:
            if info in picks:
                neworder.append(info)
        self._orderlist = neworder
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
            trues = []
            falses = []
            
            compares = ( self.getmollogic(opesplit[0]), self.getmollogic(opesplit[2]) )
            
            for info in self._orderlist:
                try: # added so comparisons to nonexistent metas won't crash
                    if Sdfmeta.comps[opera]( compares[0][info[0]][info[1]], compares[1][info[0]][info[1]] ):
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
            tear = Sdfmeta.leveler(string)
            if len(tear) != 2:
                raise ValueError('Weird logic: {}'.format(str(tear)))
            matheus = Sdfmeta.levjoin(tear[1][1])
            rcomindex = matheus.rfind(',')
            if rcomindex == -1:
                raise TypeError('Weird logic. No comma.')
            metatab = matheus[:rcomindex]
            numstring = matheus[rcomindex+1:]
            perindex = numstring.rfind('%')
            if perindex > 0:
                num = functions.numify(numstring[:perindex])
                per = True
            else:
                num = functions.numify(numstring)
                per = False
            trues = []
            falses = []
            
            if tear[0]=='MAX':
                reverse = True
            elif tear[0]=='MIN':
                reverse = False
            values = self.getmollogic(metatab)
            for mol in self._dictomoles:
                try:
                    molec = values[mol]
                except KeyError:
                    falses.extend([[mol, moln] for moln in self._dictomoles[mol]])
                    warnings.warn(mol + ' marked false for not having necesary metafield.')
                    continue
                if len(set(lmap(len,list(molec.values()))) - {0,1})>0:
                    warnings.warn('Comparing lists, might lead to weird results.',UserWarning)
                moles = OrDi(sorted([(key, molec[key]) for key in molec], key= lambda xx: xx[1], reverse = reverse))
                if per:
                    grab = int(math.ceil(len(moles)*num/100.0))
                else:
                    grab = int(num)
                trues.extend( [[mol, item] for item in list(moles.keys())[:grab]] )
                falses.extend( [[mol, item] for item in list(moles.keys())[grab:]] )
            return (trues, falses)
        
        def uniqu( string ):
            '''
            will remove duplicates of conformations by given metastatement
            '''
            metas = self.getmollogic(string[1:])
            trues=[]
            falses=[]
            values = dict()
            for molname, confnum in self._orderlist:
                if not molname in values:
                    values[molname] = []
                try:
                    if metas[molname][confnum]._data not in values[molname]:
                        values[molname].append(metas[molname][confnum]._data)
                        trues.append([molname,confnum])
                        
                    else:
                        falses.append([molname,confnum])
                except KeyError:
                    warnings.warn('No meta {} in molecule {}{{[{}]}}'.format(string[1:], molname, confnum))
                    falses.append([molname,confnum])
            return (trues, falses)
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
        
        opesplit = Sdfmeta.compsplit(string, comps=('==','<=','>=','!=','<','>'))
        
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
        
        metadic = dict()
        mylevel = Sdfmeta.levels(string)
        
        def confloop(mol):
            #Loop through conformations of single type  molecule
            precalc = dict()
            for confnum in self._dictomoles[mol]:
                try:
                    newmet = self._dictomoles[mol][confnum].molelogic(mylevel, precalc, self._dictomoles[mol])
                except (ValueError, AttributeError):
                    newmet = None
                if newmet and len(newmet)>0:
                    metadic[mol][confnum] = newmet
            if len(metadic[mol])==0:
                del(metadic[mol])
        
        def molloop():
            #Loop through types of molecules
            for mol in self._dictomoles:
                metadic[mol] = dict()
                confloop(mol)
        
        molloop()
        return metadic
    
    def orderbymeta(self,sortstring):
        '''
        Compares values of meta. If more than 1 element in requested meta, anything might happen.
        '''
        #sortstring=sortstring.strip()
        sortstring=sortstring.strip('\"|\' ')
        if sortstring[0]=='>':
            rever=True
        elif sortstring[0]=='<':
            rever=False
        else:
            return self._orderlist
        #print rever, sortstring[1:]
        #metacomp = Sdffile.getmollogic( sortstring[1:] )
        metacomp = self.getmollogic( sortstring[1:] )
        #print metacomp
        test = lambda avain: functions.numify(metacomp[avain[0]][avain[1]].getvalues()) if (avain[0] in metacomp and avain[1] in metacomp[avain[0]]) else None
        return sorted(self._orderlist, key=test, reverse=rever)
        #return sorted(self._orderlist, key=lambda avain: numify(self._dictomoles[avain[0]][avain[1]].getmeta(sortstring[1:])[0]) if self._dictomoles[avain[0]][avain[1]].hasmeta(sortstring[1:]) else None, reverse=rever)
        
    def sortme(self, sortstring):
        self._orderlist = self.orderbymeta(sortstring)
    
    @staticmethod
    def readcsv(path, separators = [';',',','\t']):
        #chop = re.compile('\s*{}\s*'.format(''.join(separators)))
        with open(path) as f:
            csvdata = [[cell.strip('"\'') for cell in Sdfmeta.compsplit(line, separators)[::2]] for line in f] #lukee tiedoston, splittaa pilkuista ja poistaa alkioista rivinvaihdot
        return csvdata
    
    @staticmethod
    def xreadcsv(path, separators = [';',',','\t']):
        #chop = re.compile('\s*{}\s*'.format(''.join(separators)))
        with open(path) as f: 
            for line in f:
                yield [cell.strip('"\'') for cell in Sdfmeta.compsplit(line, separators)[::2]] 
    
    @staticmethod
    def headToCol(header, path, head, messageHead=''):
        headn = header.index(head) if isinstance(head, string_types) and head in header else functions.numify( head )
        if not (isinstance(headn, int) and headn >= 0):
            if messageHead != '':
                messageHead = messageHead + ' '
            #print headn
            #print head
            #print header
            raise IndexError('No {}{} in {}!'.format(messageHead,head, path))
        return headn
    
    def addCsvMeta(self, strArg, **kwargs):
        '''
        Add metadata from given csv-file path
        '''
        
        liargs, diargs = kwarger(strArg)
        path = liargs[0]
        kwargs = dict(kwargs)
        for key in diargs:
            kwargs[key] = diargs[key]
        
        csvdata = Sdffile.xreadcsv(path) #separators?
        
        #header = csvdata[0]
        header = next(csvdata)
        
        ckey    = kwargs.get('confkey', Sdfmole.ckey)
        molcol = Sdffile.headToCol(header, path, kwargs.get('molcol', 0), 'molname field')
        
        molcol = header.index(molcol) if isinstance(molcol, string_types) and molcol in header else molcol
        if isinstance(molcol, string_types):
            raise IndexError('No {} in {}!'.format(molcol, path))
        ckeyindex = header.index(ckey ) if ckey in header else None
        
        for csvmol in csvdata:
            m = Sdfmole.confchop.search(csvmol[molcol])
            if m:
                n = m.group(0)[2:-2]
                seeker = csvmol[molcol][:-(len(n)+4)]
            elif ckeyindex:
                seeker = csvmol[molcol]
                n = csvmol[ckeyindex]
            else:
                seeker = csvmol[molcol]
                n = None
            
            if not seeker in self._dictomoles:
                continue
            
            if n is not None:
                keys = [n]
            else:
                keys = list(self._dictomoles[seeker].keys())
            
            for key in keys:
                if key in self._dictomoles[seeker]:
                    mol = self._dictomoles[seeker][key]
                else:
                    continue
                for i, value in enumerate(csvmol):
                    if i != molcol and len(value)!=0:
                        newmeta = Sdfmeta.construct(value)
                        mol.addmeta(header[i], newmeta, overwrite=True)
                        #mol.getmeta(header[i]).cleandelim(True)
        #del(csvdata)
        
    def addAtomicCsvMeta(self, strArg, **kwargs):
        '''
        Add metadata from given atomic csv-file path
        '''
        
        liargs, diargs = kwarger(strArg)
        path = liargs[0]
        kwargs = dict(kwargs)
        for key in diargs:
            kwargs[key] = diargs[key]
        
        csvdata = Sdffile.xreadcsv(path)
        
        header = next(csvdata)
        
        ckey    = kwargs.get('confkey', Sdfmole.ckey)
        atomn   = Sdffile.headToCol(header, path, kwargs.get('atomnumber', 'atom_number'), 'atom number field')
        molcol  = Sdffile.headToCol(header, path, kwargs.get('molcol', 0), 'molname field')
        
        '''
        molcol = header.index(molcol) if isinstance(molcol, string_types) and molcol in header else functions.numify( molcol )
        if not (isinstance(molcol, int) and molcol > 0):
            raise IndexError('No molname field {} in {}!'.format(molcol, path))
        
        atomn = header.index(atomn) if isinstance(atomn, string_types) and atomn in header else functions.numify( atomn )
        if not (isinstance(atomn, int) and atomn > 0):
            raise IndexError('No atom number field {} in {}!'.format(atomn, path))
        '''
        
        ckeyindex = header.index(ckey ) if ckey in header else None
        
        entries = dict()
        
        #if not atomn in header:
        #    raise IndexError('{} not found in {}!'.format(atomn, path))
        #    return
        #else:
        #    atomnindex = header.index(atomn)
        
        for csvatom in csvdata:
            if len(csvatom) != len(header):
                raise IndexError('Number of columns not constant in {}!'.format(path))
            an = functions.numify(csvatom[atomn])
            
            m = Sdfmole.confchop.search(csvatom[molcol])
            if m:
                n = m.group(0)[2:-2]
                seeker = csvatom[molcol][:-(len(n)+4)]
            elif ckeyindex:
                seeker = csvatom[molcol]
                n = csvatom[ckeyindex]
            else:
                seeker = csvatom[molcol]
                n = None
            
            if not seeker in self._dictomoles:
                continue
            
            if n is not None:
                keys = [n]
            else:
                keys = list(self._dictomoles[seeker].keys())
            
            for key in keys:
                if key in self._dictomoles[seeker]:
                    if not (seeker, key) in entries.keys():
                        entries[(seeker, key)] = dict([(head,OrDi()) for head in header[1:] if head not in (ckey, atomn)])
                    for k, head in enumerate(header):
                        if k in (molcol, atomn) or head == ckey:
                            continue
                        if csvatom[k].strip() != '':
                            entries[(seeker, key)][head][an] = csvatom[k]
                else:
                    continue
        
        for (mol, conf) in entries.keys():
            for meta in entries[(mol, conf)]:
                if len(entries[(mol, conf)][meta])!=0:
                    self._dictomoles[mol][conf].addmeta(meta, entries[(mol, conf)][meta], overwrite = True)
    
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
        for mol in self:
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
        if writetype in ('none','donotprint'):
            return
        if not 'split' in kwargs:
            if 'path' in kwargs: #output or overwrite
                if 'makefolder' in kwargs:
                        del(newargs['makefolder'])
                        onepath = kwargs['makefolder']+'/'+kwargs['path']
                        functions.ensure_dir(onepath)
                else:
                    onepath = kwargs['path']
                del(newargs['path'])
                f = open(onepath,'w')
                f.write(self.selftostring(writetype, **newargs))
                f.close()
                
            else: #stdout
                f = sys.stdout
                f.write(self.selftostring(writetype, **newargs))
                f.flush()
        else: #means split
            del(newargs['split'])
            if not 'path' in kwargs:
                self.writer(writetype, **newargs)
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
                        functions.ensure_dir(onepath)
                    newargs['path'] = onepath
                    onesdf.writer(writetype, **newargs)
    
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
        mymol2 = mol2.Mol2File(path)
        for mol in mymol2:
            newdic = mol.pickatomdata(column)
            name = Sdfmole.confchop.sub('',  mol['MOLECULE'][0])
            if name != mol['MOLECULE'][0] :
                conf = Sdfmole.confchop.search(mol['MOLECULE'][0]).group().strip('{[]}')
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
        del(mymol2)
        
    def injectMol2Data(self, path,  outpath, column, metastatement, defaultValue=0.0, precision=4):
        '''
        Inject atomwise information from current .sdf-file metadata to given wanted column in .mol2-file and make wanted output file.
        '''
        molli2 = mol2.Mol2File(path)
        metas = self.getmollogic(metastatement)
        
        for i, mol in enumerate(molli2):
            molname, conf = self._orderlist[i]
            mol.injectatomdata(metas[molname][conf], column, defaultValue, precision)
            del(metas[molname][conf])
        if not outpath:
            outpath=path
        molli2.writefile(outpath)
        
    def injectMol2DataStr(self, string):
        try:
            injargs = re.split('\s*,\s*', string)
        except ValueError:
            raise ValueError('Give: input, output, column, metaname, default, precision')
        
        for i, things in enumerate(zip((str, str, int, str, functions.numify, int),injargs)):
            injargs[i] = things[0](things[1]) 
        self.injectMol2Data(*injargs)
    
#end of Sdffile

class Sdfmole(object):
    '''
    Class containing all information concerning a single conformation
    '''
    mes='Conformation number mismatch! {} vs {}'
    confchop= re.compile('\{\[.+\]\}') #re.compile('\{\[\d+\]\}') #gets conformation number #changed from number to everything
    atomcut=(0, 10, 20, 30, 31, 34, 36, 39, 42, 45, 48, 51, 54, 57, 60, 63, 66, 69)
    ckey = "confnum"
    
    def __init__(self, stringsofone=None, **kwargs):
        '''
        Initialize an empty molecule with no data at all
        Adds data if list of strings containing lines of sdf-file describing single molecule is given
        '''
        self._name = ''
        self._meta = dict()
        self._metakeys = list()
        self._other = list()
        self._numeric = False
        #self._ignoretype = []
        #self._ignoretype = kwargs.get('ignores', ['H'])
        self.setIgnores(kwargs.get('ignores', ['H']))
        
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
        new._name           = self._name
        new._meta           = dict(self._meta)
        new._metakeys       = list(self._metakeys)
        new._numeric        = self._numeric
        new._ignoretype     = list(self._ignoretype)
        new._comment        = list(self._comment)
        if new.numeric:
            new._counts     = list(self._counts)
            new._atoms      = list(self._atoms)
            new._bonds      = list(self._bonds)
            new._properties = list(self._properties)
        else:
            new._other      = list(self._other)
        return new
        
    def __deepcopy__(self,memo):
        '''
        Deep copy method
        '''
        new = Sdfmole()
        new._name           = self._name
        new._meta           = copy.deepcopy(self._meta,memo)
        new._metakeys       = list(self._metakeys)
        new._numeric        = self._numeric
        new._ignoretype     = list(self._ignoretype)
        new._comment        = list(self._comment)
        if new._numeric:
            new._counts     = list(self._counts)
            new._atoms      = copy.deepcopy(self._atoms,memo)
            new._bonds      = copy.deepcopy(self._bonds,memo)
            new._properties = list(self._properties)
        else:
            new._other      = list(self._other)
        return new
    
    def initialize(self, strings):
        '''
        Add data to molecule from list of strings as in .sdf-file
        '''
        self._name = strings[0].strip()
        self._comment = [line.strip() for line in strings[1:3]]
        first = -1
        
        for i, line in enumerate(strings[1:]):
            if Sdfmeta.metaname.match(line):
                first = i+1
                break
            else:
                #self._other.append(line.rstrip('\n'))
                self._other.append(line.rstrip())
        
        add = first
        for i in range(first,len(strings)):
            #if strings[i] == '\n':
            if strings[i].rstrip() == '':
                if i-add < 2:
                    add = i+1
                    continue
                newmeta = Sdfmeta(strings[add:i+1])
                add = i+1
                if newmeta._data is not None:
                    self.addmeta(newmeta.getname(), newmeta)
                
    def dists(self, point1, **kwargs):
        '''
        return list of distances and atom numbers. single line describes distance of an atom in molecule to given point.
        Atom types in ignores are omitted.
        '''
        myignores = [item.upper() for item in kwargs.get('ignores',self._ignoretype)]
        dictomo = kwargs.get('dictomo',dict())
        precalc = kwargs.get('precalc',dict())
        
        self.numerize()
        coord1 = Sdfmole.coorder(point1)
        #if coord1 == None:
        if coord1 is None:
            self.molelogic(Sdfmeta.levels(point1+'[]'), precalc, dictomo)
            coord1 =  self.getatomloc(1) #Why 1? #FIXME
        alist=[]
        for anum, coord2 in self.atomsGenerator(ignores=myignores):
            dist = sum((coord1-coord2)**2)**0.5
            alist.append([float(dist),anum])
        return alist
    
    def calcEscapeNumberOrder(self,matcher,maxRange=2.0, **kwargs):
        '''
        calculate number of atoms in [self] that are farther than [range] from at least one atom in [other]. Returns number of atoms inside and outside.
        '''
        try:
            from findable import Findable
        except ImportError:
            from sdfconf.findable import Findable
        
        #ignores=['H'] if not 'ignores' in kwargs else kwargs['ignores']
        myignores = kwargs.get('ignores',['H'])
        #maxn=0 if not 'maxn' in kwargs else kwargs['maxn']
        maxn= kwargs.get('maxn',0)
        
        finder = {Findable:lambda : matcher, Sdfmole: lambda : Findable(matcher, ignores=myignores), mol2.Mol2Mol: lambda : Findable(matcher, ignores=myignores)}.get(type(matcher), lambda : None)()
        if not finder:
            raise TypeError('Type must be of type Findable, Sdfmole or Mol2mol.')
        
        if kwargs.get('anums',False):
            adder = lambda add: [add]
            outCount = []
            inCount = []
            tester = lambda count : len(count)
        else:
            adder = lambda add: 1
            outCount = 0
            inCount = 0
            tester = lambda count : count
        
        if maxn>0:
            test = lambda : tester(inCount)>=maxn
        elif maxn<0:
            test = lambda : tester(outCount)>=-maxn
        else:
            test = lambda : False
        
        for sAn, sCoord in self.atomsGenerator(ignores=myignores):
            atoms = finder.findinrange(sCoord, maxRange)
            if len(atoms)>0:
                inCount += adder(sAn) #1
            else:
                outCount += adder(sAn) #1
            if test():
                break
        return inCount, outCount
    
    
    def atomsGenerator(self,**kwargs):
        '''
        generator that yields tuples including (atom number, coordinates of atom) for atoms in self, except for those of type represented in ignores.
        '''
        self.numerize()
        gettab = [lambda i: i+1, lambda i: self.getatomloc(i+1)]

        ignores = kwargs.get('ignores', ['H'])
        
        if kwargs.get('types',False) :
            gettab.append(lambda i: self.gettype(i+1).strip())
        
        if kwargs.get('tabs',False):
            tablist = kwargs.get('tabs')
            if not isinstance(tablist, (list, tuple)):
                tablist = (tablist, )
            self.numerize()
            gettab[1:]=[lambda i: functions.numify(self._atoms[i][j-1]) for j in tablist]
            self.numerize()
        
        for i in range(len(self._atoms)):
            if self.gettype(i+1).strip() in ignores:
                continue
            else: 
                yield tuple(gette(i) for gette in gettab)

    def numerize(self):
        '''
        Turn initial textual format of the atomblock to listlike and numerical format.
        '''
        if not self._numeric:
            self._counts = [self._other[2][i:i+3].strip() for i in range(33)[::3]]+[self._other[2][33:39].strip()]
            blockdiv = [3,3+int(self._counts[0]),3+int(self._counts[0])+int(self._counts[1])]
            self._atoms = [[atom[Sdfmole.atomcut[i]:Sdfmole.atomcut[i+1]] for i in range(len(Sdfmole.atomcut)-1)] for atom in self._other[3:blockdiv[1]]]
            self._bonds = [[bond[i:i+3].strip() for i in range(7*3)[::3]] for bond in self._other[blockdiv[1]:blockdiv[2]]]
            self._properties = self._other[blockdiv[2]:]
            del(self._other)
            self._numeric = True
        
    def getconf(self):
        '''
        return conformation number in the name and meta for the conformation
        '''
        ans = [None, None]
        m = Sdfmole.confchop.search(self._name)
        if m:
            ans[0] = str( m.group(0)[2:-2] )
        if Sdfmole.ckey in self._meta:
            ans[1] = str( self.getmeta(Sdfmole.ckey)[0] )
        if ans[0] and ans[1]:
            if ans[0]!=ans[1]:
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
                    print(Sdfmole.mes.format(already[0],already[1]))
            else:
                self.addmeta(Sdfmole.ckey, conf, literal = True)
        if bolist[0]:
            if already[0]:
                if already[0] != conf:
                    print(Sdfmole.mes.format(already[0], conf))
            else:
                self._name = self._name + '{[' + conf + ']}'
                
    def remconf(self, bolist = [True, True]):
        '''
        remove conformation number
        bolist =(fromName, fromMeta confnum)
        '''
        already = self.getconf()
        if bolist[0] and already[0]:
            self._name = Sdfmole.confchop.sub('',self._name)
        if bolist[1] and already[1]:
            self._metakeys.remove(Sdfmole.ckey)
            self._meta.pop(Sdfmole.ckey)
            
    def setIgnores(self, ignores=['H']):
        #print ignores
        self._ignoretype = list(ignores)
            
    def getname(self):
        '''
        return name of molecule without possible conformation number
        '''
        return Sdfmole.confchop.sub('', self._name)
    
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
        same = [Sdfmole.confchop.sub('', self._name) == Sdfmole.confchop.sub('', othersdfmol._name), self.getconfn() == othersdfmol.getconfn()]
        return same
    
    def make_other(self):
        '''
        return reconstructed atomblock from numerized molecule
        '''
        other = []
        for line in self._comment:
            other.append(line + '\n')
        other.append(functions.listtostring(self._counts[:-1],3)+'{:>6}'.format(self._counts[-1])+'\n')
        for atom in self._atoms:
            other.append(Sdfmole.atomtostring(atom)+'\n')
        for bond in self._bonds:
            other.append(functions.listtostring(bond, 3)+'\n')
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
            atom=functions.numify(self.getmeta(metaatom))
            manyatom=self.getmeta(metalist)
        except KeyError:
            return None
        atoms = [functions.numify(one) for one in re.split('\s+',manyatom.strip())]
        dists = []
        for oneatom in atoms:
            dists.append(self.atomdist(atom,oneatom))
        return (atoms,dists)
    
    @staticmethod
    def atomtostring(tab):
        d=[]
        for i in range(len(Sdfmole.atomcut)-1):
            n=Sdfmole.atomcut[i+1]-Sdfmole.atomcut[i]
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
        return numpy.array(lmap(functions.numify, self._atoms[n-1][0:3]))
    
    def atomdist(self, atom1, atom2):
        '''
        return distance between 2 atoms
        '''
        self.numerize()
        return sum((self.getatomloc(atom1)-self.getatomloc(atom2))**2)**0.5
    
    def bonddists(self, atom1, intrests=None):
        '''
        Doesn't work right now
        #dijikstra algorithm for finding path length to all atoms, or to list of intresting atoms
        '''
        self.numerize()
        n=len(self._atoms)
        visited = list(numpy.zeros(n,dtype=bool))
        distances = list(numpy.ones(n,dtype=int)*numpy.inf)
        atom=functions.numify(atom1)
        distances[atom-1]=0
        if intrests:
            intrests = lmap(functions.numify, intrests)
        while True:
            links = []
            for bond in self._bonds:
                if functions.numify(bond[0])==atom:
                    links.append(functions.numify(bond[1]))
                if functions.numify(bond[1])==atom:
                    links.append(functions.numify(bond[0]))
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
    
    def getColumn(self, meta):
        if isinstance(meta, (list, tuple)):
            columns = meta
        elif isinstance(meta, int):
            columns = (meta,)
        elif isinstance(meta, Sdfmeta):
            columns = meta.getvalues()
        columns = [item+1 for item in columns]
        
        reva = OrDi()
        for item in self.atomsGenerator(tabs=columns, ignores = self._ignoretype):
            if len(item)==2:
                reva[item[0]]=item[1]
            else:
                reva[item[0]]=str(item[1:]).strip('() ')
        
        for key in reva:
            reva[key] = functions.numify(reva[key].strip()) if isinstance(reva[key], string_types) else functions.numify(reva[key]) 
        return reva
    
    def molelogic(self, mylevel, precalc=dict(), dictomo=dict()):
        
        if len(dictomo)==0:
            dictomo={'me':self}
        
        def mmaxmin(meta, getmax=True):
            fu = max if getmax else min
            if meta._datastruct == OrDi:
                items = meta._data.iteritems()
                la = lambda x : x[1]
                rety = OrDi
            else:
                items = meta._data
                la = lambda x : x
                rety = list
            
            return  rety((fu(items, key = la),))
        
        sortfunx = { 
                    'asc':True,
                    'des':False 
                    }
        
        def mlen(meta):
            try:
                return len(meta)
            except TypeError:
                return 0
            
        metafunx = {'len':mlen, 
                    'avg':lambda meta: functions.avg([thing for thing in meta]), 
                    'max':lambda meta: mmaxmin(meta, True), 
                    'min':lambda meta: mmaxmin(meta, False), 
                    'confcol': self.getColumn, 
                    'str': lambda meta: Sdfmeta.construct( meta.getmetastr()),
                    'int': lambda meta: Sdfmeta.metamap(int, meta),
                    'sum': lambda meta: sum((thing for thing in meta)), 
                    'prod': lambda meta : numpy.asscalar(numpy.prod([thing for thing in meta])), 
                    'rdup' : lambda meta : meta.withoutDuplicates() ,
                    'getmeta' : lambda meta : self.getmeta( ''.join(meta.getvalues()) ) ,
                    } 
        
        molfunx = { 
                   'mmax': lambda mets: max([mmaxmin(met, True) for met in mets]),
                   'mmin': lambda mets: min([mmaxmin(met, False) for met in mets]),
                   'mavg': lambda mets: functions.avg([metafunx['avg'](met) for met in mets]),
                   }
        
        maths = {'+' :sum,
                 '-' :functions.sub,
                 '*' :functions.myprod,
                 #'*' :numpy.prod,
                 #'*' :lambda x: functions.numify(numpy.prod(x)),
                 '/' :functions.div,
                 '%' :functions.remainder ,
                 '**':functions.mypow,
                 '++':Sdfmeta.metajoiner,
                 '--':Sdfmeta.metacut,
                 }
        
        pars = ('(','[','{','"',"'")
        
        def tabiter(conf, tab, par = None):
            '''
            collapse given levels tab generated from some meta expression into a single Sdfmeta or None
            '''
            
            def listope(conf,tab,par=None):
                #Work with lists in given structures
                
                if len(tab) > 1 and isinstance(tab[0], string_types):
                    if tab[0] in Sdfmeta.comps.keys() and (not isinstance(tab[1],Sdfmeta)) and tab[1] is not None:
                        slitab = tabiter(conf,tab[1:])
                        return tabiter(conf, [tab[0], slitab])
                
                if len(tab)==1: #evaluate
                    return tabiter(conf, tab[0], par)
                elif len(tab)==2: #it a slice / asc/max, or something like that
                    if isinstance(tab[0], string_types):
                        if tab[0] in pars:
                            #do slicing, etc.
                            if not par or tab[0] in ('"',"'"):
                                return tabiter(conf, tab[1], tab[0])
                        elif tab[0] in molfunx and tab[1][0] == '(':
                            #do precalc, etc
                            metastri = Sdfmeta.levjoin(tab)
                            if not metastri in precalc:
                                metas = []
                                for confi in dictomo:
                                    metas.append(tabiter(dictomo[confi],tab[1]))
                                precalc[metastri] = Sdfmeta.construct( molfunx[tab[0]](metas) ) 
                                del(metas)
                            if metastri in precalc:
                                return precalc[metastri] #start and end missing
                        elif tab[0] in sortfunx and tab[1][0] == '(':
                            #sort next tuple, etc.
                            meta = tabiter(conf, tab[1] )
                            if meta:
                                meta.sortme(sortfunx[tab[0]])
                                return meta
                            else:
                                return None
                            
                        elif tab[0] in Sdfmeta.comps : 
                            return tab
                            
                        elif tab[0] in metafunx and tab[1][0] == '(':
                            return Sdfmeta.construct( metafunx[tab[0]](tabiter(conf, tab[1] )) )
                    #slice
                    meta = tabiter(conf, tab[0])
                    if meta:
                        sli = tabiter(conf, tab[1]) #assumes tuple
                        if isinstance(sli, (slice,Sdfmeta)):
                            return meta.slicer(sli, tab[1][0])
                        elif sli: 
                            if tab[1][0] == '(':
                                meta.pickvalues(tabiter(conf, sli[1]), Sdfmeta.comps.get(sli[0], None), bykeys=False)
                                return meta
                            elif tab[1][0] == '{':
                                meta.pickvalues(tabiter(conf, sli[1]), Sdfmeta.comps.get(sli[0], None), bykeys=True)
                                return meta
                    return None
                elif len(tab)>2:
                    if len(tab)==3 and isinstance(tab[0], string_types) and tab[0] in maths:
                        #do math
                        f1 = tabiter(conf, tab[1])
                        if not f1:
                            if tab[0] not in ('++','--'):
                                f1 = Sdfmeta.construct(0)
                            #elif tab[0] == '--':
                            #    f1 = Sdfmeta()
                            else:
                                print "no f1"
                                f1 = Sdfmeta()
                            #else:
                            #    f1 = Sdfmeta()
                        if tab[0] in ('++','--'):
                            print "join/cut {}".format(tab)
                        return Sdfmeta.metaoper(maths[tab[0]],[f1,tabiter(conf, tab[2])])
                    else:
                        metaus = tabiter(conf, tab[:2], None)
                        for item in tab[2:]:
                            metaus = tabiter(conf, [metaus,item])
                        return metaus
                else:
                    return None
                
            def striope(conf,tab,par=None):
                #Work with strings in given structure
                if par in ('"',"'"):
                    return Sdfmeta.construct(tab,literal=True)
                elif tab in (sortfunx, molfunx, metafunx) :
                    return tab
                elif conf.hasmeta(tab):
                    return copy.copy(conf.getmeta(tab))
                elif par in ('(','[','{'):
                    trytab = [functions.numify(i) for i in re.split('\s*[ ,]{0,1}\s*', tab )]
                    if not str in lmap(type, trytab):
                        return Sdfmeta.construct( trytab )
                    else:
                        try:
                            return slice(*[{True: lambda n: None, False: int}[x == ''](x) for x in (tab.split(':') + ['', '', ''])[:3]])
                        except ValueError:
                            rip = Sdfmeta.compsplit(tab)
                            if len(rip)>1 and par in ('(','{'):
                                return rip
                            raise ValueError('Your logic makes no sense. '+str(tab))
                else:
                    trytab = [functions.numify(i) for i in re.split('\s*[ ,]{0,1}\s*', tab )]
                    if not str in lmap(type, trytab):
                        return Sdfmeta.construct( trytab )
                return None
                    
            def raiser(conf, tab, par=None):
                raise TypeError('Bad levels: '+str(tab))
            
            #testdi = {list: listope, tuple: listope, str: striope, Sdfmeta: lambda conf, me, par=None : me if len(me)>0 else None, type(None) : lambda conf, me, par=None : None }
            testdi = {list: listope, tuple: listope, string_types: striope, Sdfmeta: lambda conf, me, par=None : me if len(me)>0 else None, type(None) : lambda conf, me, par=None : None }
            for key in testdi:
                if isinstance(tab, key):
                    return testdi[key](conf, tab, par)
            return raiser(conf, tab, par)
            #return testdi.get(type(tab), raiser )(conf, tab, par)
        
        return tabiter(self, mylevel)
    
    def levopemap(self, tab, par=None, length=None):
        '''
        partly deprecated
        Accepts lists made by leveler method. Maps this list for mathematical operators and metafield names. Also slices.
        '''
        if isinstance(tab, tuple):
            #if par and par in ('"',"'"):
            if tab[0] in ('"',"'"):
                return Sdfmeta.construct(Sdfmeta.levjoin(tab[1]))
                #return levjoin(tab[1])
            else:
                return (tab[0], self.levopemap(tab[1], tab[0]))
        elif isinstance(tab, list):
            for j, thing in enumerate(tab):
                if not isinstance(thing, string_types):
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
                tab = functions.numify(tab)
                if not isinstance(tab, string_types):
                    return Sdfmeta.construct(tab)
                else:
                    raise ValueError('Your logic makes no sense: '+tab)
    
    @staticmethod
    def coorder(point):
        goodchop = re.compile('\s*,{0,1}\s*') #CSV-separator
        coord = numpy.array( lmap(functions.numify, goodchop.split(point.strip('{[()]}'))) )
        if type(coord) == numpy.ndarray and len(coord) > 1:
            return coord
        else:
            return None
    
    def stripbutmeta(self, metaobj):
        '''
        Remove atoms except those defined in given liststructure
        '''
        self.numerize()
        #atoms = self.logicgetmeta(levemeta)
        #atoms = sorted(metaobj.getvalues())
        atoms = metaobj.getvalues()
        
        newatoms = []
        
        toNew = OrDi()
        toOld = OrDi()
        
        #for n in atoms._data:
        for i, n in enumerate(atoms):
            toNew[int(n)]=i+1
            toOld[i+1]=int(n)
            newatoms.append(self._atoms[n-1])
        self._atoms = newatoms
        self._bonds = []
        self._counts[0] = len(newatoms)
        self._counts[1] = 0
        
        self.addmeta('toNew', toNew)
        self.addmeta('toOld', toOld)
        

#end of Sdfmole

class Sdfmeta(object):
    '''
    Class including information in single metafield
    '''
    
    metaname = re.compile('\>(.*)\<(.+)\>') #Match gets metafield name
    ematch =   re.compile('[ ]{0,2}')
    fisep  =   re.compile('[ ,;]')
    comps = OrDi( ( ('>=',operator.ge), ('<=',operator.le), ('<',operator.lt), ('>',operator.gt), ('==',operator.eq), ('!=',operator.ne) , ('=',operator.eq) ) )
    
    def __init__(self, listofstrings=None):
        '''
        Initializes an empty metafield
        If listofstrings is given (a single metafield in .sdf-file) adds that data
        '''
        self._name = None #Metafield name
        self._datatype = None #int, float, 'str'
        self._datastruct = None #list, dict, single
        self._data = None #The actual data
        self._delims = []
        self._dumb = True
        #self._dumbcontent = []
        if listofstrings:
            self.initialize(listofstrings)
        
    def __getitem__(self, ind):
        '''
        get a slice from list or dict Sdfmeta
        '''
        mylist = self.getvalues()
        if isinstance( ind, (slice, int) ):
            return mylist[ind]
        else:
            raise TypeError("Invalid argument type.")
            
    def __iter__(self):
        return iter(self.getvalues()
    )
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
            othermeta = Sdfmeta.construct( functions.numify(other))
        else:
            othermeta = other
        li1 = list(self.getvalues())
        li2 = list(othermeta.getvalues())
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
            othermeta = Sdfmeta.construct( functions.numify(other))
        else:
            othermeta = other
            othermeta.numerize()
        
        li1 = list(self.getvalues())
        li2 = list(othermeta.getvalues())
        
        #print '{}, {}'.format(self.dtype(), other.dtype())
        #if self.dtype() == 'str' and other.dtype() == 'str':
        if self.isStr() and other.isStr():
            #print 'bang'
            def compar(string1, string2):
                test = lambda x: x[:6].upper() == 'REGEX:'
                if test(string1):
                    return bool(re.search(string1[6:], string2 ))
                elif test(string2):
                    return bool(re.search(string2[6:], string1 ))
                else:
                    return operator.eq(string1, string2)
            if oper == operator.eq:
                oper = compar
            elif oper == operator.ne:
                oper = lambda x,y : not compar(x,y)
            #else crash
        for item2 in li2:
            for item1 in li1:
                if oper(item1, item2):
                    return True
        return False
    
    def __str__(self):
        return "Sdfmeta({})".format(str(self._data))
    
    def __repr__(self):
        return "Sdfmeta,{},{},{}".format(self.dtype(), self._datastruct, str(self._data))
        
    def initialize(self, listofstrings):
        '''
        parse the metadata from list of strings
        '''
        #get name of metafield
        if listofstrings[0][0] != '>' :
            self._name = None
            fi = 0
        else:
            self._name = Sdfmeta.metaname.match(listofstrings[0]).groups()[1]
            fi = 1
        if len(listofstrings[-1].strip())==0:
            self._data = [listofstrings[i].strip('\n ') for i in iterange(fi,len(listofstrings)-1) ]
        
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
                if not Sdfmeta.fisep.match(line[-1]):
                    newlines.append(' ')
            newlines.append(mylines[-1])
            (data,dtype,delims) = Sdfmeta.whattype(''.join(newlines))
    
            #if string, it's special
            if dtype == str and type(data) != OrDi:
                #self._datatype = 'str'
                self.setType(dtype)
                #self._data = [line.strip('\n') for line in mylines]
                self._data = [line.rstrip() for line in mylines]
                self._datastruct = list
                self._delims = ['' for line in mylines[:-1]]
                #return
            #not just string
            else:
                self._data = data
                #self._datatype = dtype if dtype != str else 'str'
                self.setType(dtype)
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
        else:
            name = None
        new.setname(name)
        
        if isinstance(data, string_types):
            
            if 'literal' in dictarg and dictarg['literal']:
                new._datastruct = 'single'
                data = [data]
            else:
                (newdata, newtype, newdelims) = Sdfmeta.whattype(data)
                if newtype != str:
                    return Sdfmeta.construct(newdata, delims = newdelims, name = name)
                else:
                    new._datastruct = 'single'
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
        elif isinstance(data, Sdfmeta):
            out = copy.deepcopy(data)
            #out.numerize()
            out.setname(name)
            if 'delims' in dictarg:
                out._delims = dictarg['delims']
            return out
        else:
            print(type(data))
            raise TypeError('Datastructure type not list, dict, OrderedDict, str, int or float.')
        if not dictarg.get('literal',False): #('literal' in dictarg and dictarg['literal']):
            if type(data) == list:
                data = lmap(functions.numify, data)
                
            elif type(data) == OrDi:
                for key in data:
                    newkey = functions.numify(key)
                    data[newkey] = functions.numify(data[key])
                    if key != newkey:
                        del(data[key])
                
        if type(data) == OrDi:
            types = {type(data[key]) for key in data}
        else:
            types = set(lmap(type, data))
        
        new._data = data
        if not all({typ in (int, float, str) for typ in types}):
            raise TypeError('Data type not str, int or float.')
        if len(types)==1:
            #new._datatype = iter(types).next()
            #new._datatype = next( iter(types) )
            new.setType(next( iter(types) ))
            #if new._datatype == str:
            #    new._datatype = 'str'
        elif len(types)>1 and str in types:
            raise TypeError('Mixed datatypes. Strings and numbers.')
        else:
            #new._datatype = float
            new.setType(float)
            if type(new._data) == list:
                new._data = lmap(float, new._data)
            elif type(new._data) == OrDi:
                for key in new._data:
                    new._data[key] = float(new._data[key])
            else:
                raise TypeError('DaFaq?')
        if type(new._data) in (list, OrDi) and not 'delims' in dictarg:
            new._delims =  [ ', ' ]*(len(new._data)-1) 
        return new
        
    def typeCaster(self, value, casttype=None):
        if not casttype:
            casttype = self.dtype()
        if casttype == 'str' or not casttype:
            return str
        else:
            return self.dtype()
    
    @staticmethod
    def leveler(string):
        parre = re.compile('[\{\[\(\)\]\}\"\']') #Matches all parentheses
        pair = {'"':'"',"'":"'",'(':')','[':']','{':'}'}
        pars = [(item.group(), item.start()) for item in parre.finditer(string)]
        
        def deslash(pali,string):
            for i, par in enumerate(pali):
                if par[1]-1>=0  and string[par[1]-1]=='\\':
                    newpali, string = deslash(lmap(lambda x: (x[0],x[1]-1),pali[i+1:]),string[:par[1]-1]+string[par[1]:])
                    pali = pali[:i] + newpali
                    break
            return pali, string
        
        pars, string = deslash(pars, string)
        
        tab = []
        
        if len(pars)<2:
            return [string]
        
        curpar = None
        backpar = None
        
        if pars[0][1]>0:
            stastring = string[:pars[0][1]]
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
                    if curpar[0] in ("'",'"'):
                        tab.append((curpar[0],string[curpar[1]+1:item[1]]))
                        curpar = None
                        backpar = item
                    else:
                        tab.append((curpar[0],Sdfmeta.leveler(string[curpar[1]+1:item[1]])))
                        curpar = None
                        backpar = item
            elif item[0] == curpar[0]:
                level += 1
        endstring = string[pars[-1][1]+1:]
        if len(endstring)>0:
            tab.append(endstring)
        return tab
    
    @staticmethod
    def dumb_levopemap(tab, par=None, length=None):
        '''
        Accepts lists made by leveler method. Maps this list for mathematical operators and metafield names. Doesn't understand local metafields. a dumb version.
        '''
        #print tab
        if isinstance(tab, tuple):
            if tab[0] in ('"',"'"):
                return (tab[0], Sdfmeta.levjoin(tab[1]))
            else:
                return (tab[0], Sdfmeta.dumb_levopemap(tab[1], tab[0]))
        elif isinstance(tab, list):
            
            searches = lmap(re.compile, ('[^\+](\+)[^\+]|[^-](-)[^-]', '[^\*](\*)[^\*]|(/)|(%)', '(\*{2})|(\+{2})|(-{2})'))
            
            mycomps = list(Sdfmeta.comps.keys())
            mycomps.remove('=')
            
            for j, thing in enumerate(tab):
                if not isinstance(thing, string_types):
                    continue
                #for comp in Sdfmeta.comps:
                for comp in mycomps:
                    ma = re.match(comp,thing)
                    if ma and len(thing) > len(ma.group()):
                        interpret = Sdfmeta.dumb_levopemap([thing[len(ma.group()):]] + tab[j+1:])
                        if not isinstance(interpret, list):
                            interpret = [interpret]
                        return tab[:j] + [ma.group(), interpret ]
                    
                thing = thing.strip()
                thing = ' ' + thing + ' '
                sear=None
                for s in searches:
                    sear = list(s.finditer(thing))
                    if len(sear)>0:
                        mygroup = next(item for item in sear[-1].groups() if item is not None)
                        i = sear[-1].start() + sear[-1].group().find(mygroup)
                        ie = i +len(mygroup)
                        break
                
                if sear:
                    if i>0 and thing[i-1]=='\\' :  #Escape character for math operators
                        tab[j] = thing[1:i-1]+thing[i:-1]
                        continue
                    tab1=[]
                    tab1.extend(tab[:j])
                    if i > 1:
                        tab1.append(thing[1:i])
                    if ie < len(thing)-1:
                        tab2 = [thing[ie:-1]]
                    else:
                        tab2 = []
                    tab2.extend(tab[j+1:])
                    return (thing[i:ie], Sdfmeta.dumb_levopemap(tab1), Sdfmeta.dumb_levopemap(tab2))
                #new ends
            return [Sdfmeta.dumb_levopemap(item, par, len(tab)) for item in tab]
        else:
            return tab.strip()
    
    @staticmethod
    def levels(string):
        return Sdfmeta.dumb_levopemap(Sdfmeta.leveler(string))
    
    @staticmethod
    def levjoin(taber):
        pair = {'(':')','[':']','{':'}',"'":"'",'"':'"'}
        string = []
        for item in taber:
            if isinstance(item, string_types):
                string.append(item)
            elif isinstance(item, list):
                string.append(Sdfmeta.levjoin(item))
            elif isinstance(item, tuple):
                if item[0] in pair:
                    string.extend([ item[0], Sdfmeta.levjoin(item[1]), pair[item[0]] ])
                else:
                    string.extend([Sdfmeta.levjoin(item[1][0]), item[0], Sdfmeta.levjoin(item[1][1]) ])
        return ''.join(string)
    
    @staticmethod
    def compsplit(string, comps = ('==','<=','>=','!=','=>','=<','<','>','=')):
        tab = Sdfmeta.leveler(string)
        ret = []
        for i, item in enumerate(tab):
            if isinstance(item, string_types):
                j=-1
                for comp in comps:
                    j=item.find(comp)
                    if j>=0:
                        #return Sdfmeta.compsplit(Sdfmeta.levjoin(tab[:i])+item[:j]) + [ item[j:j+len(comp)]] + Sdfmeta.compsplit(item[j+len(comp):]+Sdfmeta.levjoin(tab[i+1:]))
                        return [(Sdfmeta.levjoin(tab[:i])+item[:j]).strip() , (item[j:j+len(comp)]).strip() ] + Sdfmeta.compsplit( item[j+len(comp):]+Sdfmeta.levjoin(tab[i+1:]) , comps)
        if string.strip() != '':
            ret.append(string.strip())
        return ret
    
    @staticmethod
    def metajoiner(metas, **params):
        '''
        Joins multiple metavalues into one
        params : name, delims
        '''
        newmeta = None
        for meta in metas:
            if meta is not None and len(meta)>0:
                if not newmeta:
                    newmeta = copy.copy(meta)
                else:
                    newmeta.extend(meta)
        newmeta.cleandelim(True)
        if 'name' in params and newmeta:
            newmeta.setname(params['name'])
        else:
            newmeta.setname(None)
        return newmeta
    
    @staticmethod
    def metacut(metas):
        newmeta = copy.copy(metas[0])
        for meta in metas[1:]:
            newmeta.cut(meta)
        newmeta.cleandelim(True)
        return newmeta
    
    def __copy__(self):
        '''
        Shallow copy method
        '''
        new = Sdfmeta()
        new._name = self._name
        #new._datatype = self.dtype()
        new.setType(self.dtype())
        new._datastruct = self._datastruct
        new._dumb = self._dumb
        if len(self)>0:
            new._data = type(self._data)(( self._data ))
            new._delims = list( self._delims )
        return new
        
    def __deepcopy__(self, memo):
        '''
        Deep copy method is the same as __copy__
        '''
        
        '''
        new = Sdfmeta()
        new._name = self._name
        #new._datatype = self.dtype()
        new.setType(self.dtype())
        new._datastruct = self._datastruct
        new._dumb = self._dumb
        if len(self)>0:
            new._data = type(self._data)(( self._data ))
            new._delims = list( self._delims )
        return new
        '''
        return self.__copy__()
    
    def getname(self):
        '''
        return name of Sdfmeta
        '''
        return self._name #[-1]
    
    def getvalues(self):
        '''
        return values of meta. For OrDi keys are neglected.
        '''
        if self._datastruct == OrDi:
            return list(self._data.values())
        else:
            return list(self._data)
    
    def isdumb(self):
        '''
        Has Sdfmeta been numerized yet
        '''
        return self._dumb
    
    def setname(self, newname):
        '''
        Change the name of Sdfmeta
        '''
        if isinstance(newname, string_types):
            self._name = newname
        elif not newname:
            self._name = None
        elif isinstance(newname, (tuple,list)):
            self.setname(newname[-1])
            #historical reasons
        else:
            raise TypeError('Name defined incorrectly! '+str(newname))
        
    def dtype(self):
        return self._datatype
    
    def setType(self, datatype):
        '''
        set datatype for meta.
        '''
        try:
            if datatype == 'str' or isinstance(datatype(''), string_types):
                datatype = 'str'
        except ValueError:
            pass
        except TypeError:
            pass
        self._datatype = datatype
        
    def isStr(self):
        '''
        string test
        '''
        return self.dtype() == 'str'
        
    def extend(self, other):
        '''
        Used to merge more data to a Sdfmeta
        '''
        if other is None:# or len(other)==0:
            return
        
        if other._datastruct == 'single':
            other._datastruct = list
        if self._datastruct == 'single':
            self._datastruct = list
        
        #If new Sdfmeta, do this
        if not self._datastruct:
            self._data =       copy.deepcopy( other._data )
            self._datastruct = other._datastruct
            #self._datatype =   other.dtype()
            self.setType( other.dtype() )
            self._delims =     list( other._delims )
            self._dumb =       other._dumb 
            #self._dumbcontent =list( other._dumbcontent )
            return
        
        floatflag = False
        if self.dtype() != other.dtype():
            #if self.dtype() == 'str':
            if self.isStr() :
                #other._datatype = 'str'
                other.setType(str)
                other._data = OrDi(((key,str(other._data[key])) for key in other._data)) if other._datastruct == OrDi else lmap(str, other._data) 
                self.extend(other)
            #elif other.dtype() == 'str':
            elif other.isStr() :
                #self._datatype = 'str'
                self.setType(str)
                self._data = OrDi(((key,str(self._data[key])) for key in self._data)) if self._datastruct == OrDi else lmap(str, self._data)
                self.extend(other)
                
            else:
                floatflag = True
        
        elif self._datastruct == other._datastruct:
            if self._datastruct == list:
                self._data.extend(other._data)
            elif self._datastruct == OrDi:
                self._data.update(other._data)
        
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
            if self._datastruct == OrDi:
                #self._data = {k: float(v) for k, v in self._data.items()}
                for key in self._data:
                    self._data[key] = float(self._data[key])
            else:
                self._data = lmap(float, self._data)
                
    def cut(self, other):
        '''
        Used to cut data in other Sdfmeta from meta of Sdfmeta
        '''
        
        if other._datastruct == 'single':
            other._datastruct = list
        if self._datastruct == 'single':
            self._datastruct = list
        
        #If new Sdfmeta, do this
        if not self._datastruct:
            return Sdfmeta()#.construct([])
        
        floatflag = False
        if self.dtype() != other.dtype():
            #if self.dtype() == 'str' or other.dtype() == 'str':
            if self.isStr() or other.isStr():
                return None
            else:
                floatflag = True
        
        elif self._datastruct == other._datastruct:
            if self._datastruct == list:
                for rem in other._data:
                    while True:
                        try:
                            index = self._data.index(rem)
                        except ValueError:
                            break
                        self._data = self._data[:index] + self._data[index+1:]
                        self._delims = self._delims[:index] + self._delims[index+1:]
                        
            elif self._datastruct == OrDi:
                for key in other._data:
                    try:
                        index = list(self._data.keys()).index(key)
                        del(self._data[key])
                        self._delims = self._delims[:index] + self._delims[index+1:]
                    except ValueError:
                        continue
        else:
            raise TypeError('Mixed datastructures')
        
        if floatflag:
            '''
            if self._datastruct == list:
                self._data = lmap(float, self._data)
            else:
                self._data = {k: float(v) for k, v in self._data.items()}
            '''
            if self._datastruct == OrDi:
                for key in self._data:
                    self._data[key] = float(self._data[key])
            else:
                self._data = lmap(float, self._data)
    
    def withoutDuplicates(self):
        if self._datastruct != OrDi:
            newlist = []
            for item in self.getvalues():
                if item not in newlist:
                    newlist.append(item)
        else:
            newlist = OrDi()
            for key in list(self._data.keys()):
                if self._data[key] not in newlist:
                    newlist[key] = self._data[key]
        return newlist
        
    def getmetastrings(self, length=float('inf')):
        def floattosting(flo):
            try1 = str(round(flo,4))
            if try1 == '0.0' and flo != 0:
                return '{:.4e}'.format(flo)
            else:
                return try1
        
        stringers = {int:str, float:floattosting, 'str' : str}
        strifu = stringers[self.dtype()]
        
        dictflag = self._datastruct == OrDi
        
        if dictflag:
            key = list(self._data.keys())[0]
            tmp=[str(key)+':'+strifu(self._data[key])]
            del(key)
            itera = enumerate(list(self._data.keys())[1:])
        else:
            tmp=[strifu(self._data[0])]
            itera = enumerate(self._data[1:])
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
                print(self._name)
                print(self._data)
                print(self._delims)
            if l +lde > length:
                strings.append(''.join(tmp))
                tmp=[]
                l=0
            l += lde
            tmp.append(self._delims[i])
            #add data
            lda = len(stuff)
            if l +lda > length and lda < length:
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
    
    @staticmethod
    def metamap(fu, meta):
        if fu not in (str, int, float):
            raise TypeError('Given type must be int, str or float! It was {}!'.format(fu))
        meta = Sdfmeta.construct(meta)
        if meta._datastruct == OrDi:
            for key in meta._data.keys():
                try:
                    meta._data[key] = fu(meta._data[key])
                except ValueError:
                    del(meta._data[key])
        else:
            new = []
            for item in meta.getvalues():
                try:
                    new.append(fu(item))
                except ValueError:
                    continue
            meta._data = new
        meta.setType(fu)
        meta.cleandelim()
        return meta
        
    
    def selftolistofstrings(self):
        '''return in .sdf-file format. No linechanges'''
        #Make list
        strings = []
        #Append name
        strings.append('>  <'+self._name+'>')
        
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
            self._delims = [newdelim] * (len(self._data)-1)
                
        
    @staticmethod
    def metaoper(oper, metas):
        '''
        method to make operations for Sdfmeta objects
        return resulting metavalue of operation.
        resulting Sdfmeta is nameless
        '''
        structs = [meta._datastruct for meta in metas]
        types   = [meta.dtype()   for meta in metas]
        #mathopers = (sum, functions.sub, numpy.prod, max, min, functions.avg, functions.div, functions.mypow, functions.remainder)
        mathopers = (sum, functions.sub, functions.myprod, max, min, functions.avg, functions.div, functions.mypow, functions.remainder)
        workmetas = [copy.copy(meta) for meta in metas] 
        if oper in mathopers:
            
            
            if len(set(structs))>2:
                raise TypeError('Mixed datastructs')
            
            kwargs = {}
            if 'str' in types:
                kwargs['literal'] = True
                if oper == sum:
                    for mtype in types:
                        if mtype != 'str':
                            return None
                    oper = lambda x : ''.join(x)
                elif oper == functions.sub:
                    oper = lambda x : re.sub(x[1], '', x[0])
                else:
                    return None
                
            try:
                minlen = min( {len(meta._data) for meta in workmetas}-{1} ) # shortest list length, ignore singles
            except ValueError:
                minlen = 1
                
            try:
                ostru = next(iter(set(structs)-{'single'})) #the other structuretype
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
                return Sdfmeta.construct( Sdfmeta.listoper(oper, workmetas, singles), **kwargs ) #Nameless meta
            
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
                            keyorder = list(meta._data.keys())
                
                for i,meta in enumerate(workmetas):
                    if meta._datastruct == 'single':
                        newdict = OrDi()
                        for key in keyorder:
                            if key in keys:
                                newdict[key]=meta._data[0]
                        workmetas[i]=Sdfmeta.construct(newdict, name = meta.getname(), **kwargs)
                #'''
                return Sdfmeta.construct(Sdfmeta.OrDioper(oper,workmetas), **kwargs)
        
        elif oper in (Sdfmeta.metajoiner,Sdfmeta.metacut): #Therefore, it's joiner
            newmeta = oper(metas)
            newmeta.cleandelim(True)
            return newmeta
        
        else:
            raise ValueError('Non-existent operator ')
        
    @staticmethod
    def listoper(oper, metas, singles=True):
        tab = []
        i=0
        while i<max(lmap(len, metas)):
            caltab = []
            for meta in metas:
                try:
                    caltab.append(meta._data[i])
                except IndexError:
                    if meta._datastruct == 'single' and singles: #if singles == False, chop the meta even if length == 1
                        caltab.append(meta._data[0])
                    else:
                        return tab
            tab.append(oper(caltab))
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
            dic[key] = oper(caltab) 
        return dic #OrDi(dic)
        
    @staticmethod
    def singleoper(oper, metas):
        if len({meta._datastruct for meta in metas} - {'single'})>0:
            return None
        else:
            return [oper([meta._data[0] for meta in metas])]
        
    def pickvalues(self, value, oper, **kwargs):
        if type(value)!=Sdfmeta:
            value=Sdfmeta.construct(value)
        
        
        if self.isdumb():
            self.numerize()
        
        if self._datastruct == OrDi :
            skiplist = []
            compbykey = kwargs.get('bykeys', False)
            if value._datastruct == OrDi:
                for key in self._data.keys():
                    if not key in value._data.keys():
                        skiplist.append(key)
                        continue
        
        elif self._datastruct in (list, 'single') :
            if ( not  ( (len(value) == 1) or (len(value) == len(self))) ) or value._datastruct == OrDi :
                raise functions.InputException('Valuepick: comparing list to non-single, different lenght list or OrDi-meta.')
            
        #if oper in (operator.eq, operator.ne) and (self.dtype() == 'str' and value.dtype() == 'str'):
        if oper in (operator.eq, operator.ne) and (self.isStr() and value.isStr()):
            def comparator(string1, string2):
                test = lambda x: x[:6].upper() == 'REGEX:'
                if test(string1):
                    return bool(re.search(string1[6:], string2 ))
                elif test(string2):
                    return bool(re.search(string2[6:], string1 ))
                else:
                    return (string1 == string2)
            
            if oper == operator.eq:
                compar = comparator
            elif oper == operator.ne:
                compar = lambda x,y : not comparator(x,y)
            #else crash
            #compar = comparator if oper == operator.eq else lambda x,y : not comparator(x, y)
            
        else:
            compar = oper
        
        if self._datastruct == OrDi:
            newdata = OrDi()
            if value._datastruct == OrDi:
                oval = lambda x : (value._data[x], )
            else:
                oval = lambda x : (val for val in value._data)
            
            myval = (lambda x : x) if compbykey else (lambda x: self._data[x])
            
            for key in self._data.keys():
                if key in skiplist:
                    continue
                flag = True
                myva = myval(key)
                for ova in oval(key):
                    if not compar(myva, ova):
                        flag = False
                        break
                if flag:
                    newdata[key] = self._data[key]
        else:
            newdata = []
            if len(value)>1:
                oval = lambda x : value._data[x]
            else:
                oval = lambda x : value._data[0]
            for i, thing in enumerate(self._data):
                if compar(thing, oval(i)):
                    newdata.append(thing)
        
        self._data = newdata
        self._delims = self._delims[:len(self._data)-1]
        
    def slicer(self, sliceorindex, paren):
        
        def itsslice( toget, slic ):
            if paren == '(':
                return Sdfmeta.construct(  OrDi( [ (i, toget[i]) for i in self._data.keys()[slic] ] ) )
            else:
                return Sdfmeta.construct(  toget[slic] )
        
        def itsmeta( toget, meta ):
            indexes = meta.getvalues()
            if paren == '(':
                if type(self._data) != OrDi:
                    raise TypeError('"(" not applicaple for lists')
                return Sdfmeta.construct(  OrDi([(i, toget[i]) for i in indexes if i in toget]) )
            le = len(toget)
            return Sdfmeta.construct(  [toget[i] for i in indexes if i < le ] )
        
        if paren == '{':
            toget = list(self._data.keys())
        elif paren == '[' and self._datastruct == OrDi:
            toget = list(self._data.values())
        else:
            toget = self._data
        
        if not type(toget) in (list, OrDi):
            raise TypeError('There is something wrong with our toget going to slicer') 
        slices = {Sdfmeta:itsmeta, slice:itsslice, type(None): lambda toget, sliceorindex: None }
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
        if isinstance(onestring, string_types): #type(onestring) == str:
            if functions.numitest(onestring):
                numoutof = functions.numify(onestring)
            else:
                numoutof = onestring.strip()
        else:
            numoutof = functions.numify(onestring)
        
        ty = type(numoutof)
        if ty in (int, float):
            #Number, return
            return ([numoutof],ty,[])
        #not just number, go on.
        del(numoutof,ty)
        
        
        def listtypetest(onestring, delimiter):
            #Is it a list delimited by given delimiter ?
            splitted = functions.listsep(onestring,delimiter)
            #is dict?
            dicti = OrDi() #OrderedDict
            notdict = False
            for cell in splitted[0]:
                dicted = cell.split(':')
                if len(dicted)==2:
                    dicti[functions.numify(dicted[0].strip())] = functions.numify( dicted[1].strip() )
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
                lister = [ functions.numify( cell.strip() ) for cell in splitted[0] ]
                mytype = functions.allsame(lister)
                if mytype: #all same
                    if mytype == str:
                        return None
                    elif mytype == float:
                        lister = lmap(float,lister)
                    return (lister,mytype,splitted[1])
                else:
                    pass
            else: #is dict
                mytype = functions.allsame(dicti)
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

def kwarger(argumentstr, separators = [',',';']):
    '''
    Reads a string of arguments and parses list like *args and dict like **kwargs.
    '''
    arguments = Sdfmeta.compsplit(argumentstr, separators )[::2]
    boostrings = {'True':True, 'T':True, 'False':False, 'F':False}
    args = []
    kwargs = dict()
    for arg in arguments:
        split = Sdfmeta.compsplit(arg, '=')[::2]
        if len(split) == 1:
            '''
            try:
                data = eval(split[0])
            except NameError:
                data = split[0].strip()
            '''
            data = functions.numify(split[0])
            data = boostrings.get(data, data)
            args.append(data)
            #args.append(split[0].strip())
        elif  len(split) == 2:
            '''
            try:
                data = eval(split[1])
            except NameError:
                data = split[1].strip()
            '''
            data = functions.numify(split[1])
            data = boostrings.get(data, data)
            kwargs[split[0].strip()] = data  #functions.numify(split[1])
    return args, kwargs
    
