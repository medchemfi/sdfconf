#!/usr/bin/env python
# -*- coding: latin-1 -*-

#a test brach similiar to metatolist. Just git-training for myself...
#test

import os
import sys
import re
import math
import argparse
import numpy
import time
import pylab
import copy
import operator
from collections import OrderedDict as OrDi

#Common regular expressions used.  
confchop=re.compile('\{\[\d+\]\}') #gets conformation number
molchop = re.compile('^\${4}') #Separates molecules in file
#metachop = re.compile('>\s+<') #Detects beginning of metafield name
stapa=re.compile('\{|\[|\(') #Starting parentesis
endpa=re.compile('\}|\]|\)') #Ending parenthesis
goodchop = re.compile('\s*,{0,1}\s*') #CSV-separator

metaname = re.compile('\>(.*)\<(.+)\>') #Match gets metafield name
#metadict = re.compile('[;,^]{0,1}\s*(\w+):(\w+)\s*[;,$]{0,1}' #re.compile('\s*[,;]{0,1}\s*(.*):(.*)\s*[,;]{0,1}\s*')

parre = re.compile('[\{\[\(\)\]\}]') #Matches all parentheses

ckey = "confnum"
atomcut=[0, 10, 20, 30, 31, 34, 36, 39, 42, 45, 48, 51, 54, 57, 60, 63, 66, 69]


#Parsing a slice from string
#slice(*[{True: lambda n: None, False: int}[x == ''](x) for x in (':3'.split(':') + ['', '', ''])[:3]])


class Sdffile(object):
    def __init__(self, path=''):
        self._dictomoles = dict() # {'aspirin':{1:aspi1, 2:aspi2, ...}, 'bentzene':{1:benz1, 2:benz2, ...}}
        self._orderlist = list()
        if path != '':
            self.readself(path)
            
    def __copy__(self):
        new = Sdffile()
        new._dictomoles = copy.copy(self._dictomoles)
        new._orderlist = copy.copy(self._orderlist)
        return new
        
    def __deepcopy__(self,memo):
        new = Sdffile()
        new._dictomoles = copy.deepcopy(self._dictomoles,memo)
        new._orderlist = copy.deepcopy(self._orderlist,memo)
        return new

    def __iter__(self): #iterator returns actual objects! not copies
        self.iteone = iter(self._dictomoles)
        try:
            self.itehelpone = self.iteone.next()    #CRASHES IF dictomoles is empty. Try-except takes care of it... this is awful
        except StopIteration:
            return iter([])
        self.itetwo = iter(self._dictomoles[self.itehelpone])
        return self
        

    def next(self):
        try:
            self.itehelptwo = self.itetwo.next()
            return self._dictomoles[self.itehelpone][self.itehelptwo]
            #return (self.itehelpone, self.itehelptwo, self._dictomoles[self.itehelp][self.itehelptwo]) #tuple with two keys
        except StopIteration:
            self.itehelpone=self.iteone.next()
            self.itetwo = iter(self._dictomoles[self.itehelpone])
            self.itehelptwo = self.itetwo.next()
            return self._dictomoles[self.itehelpone][self.itehelptwo]
            #return (self.itehelpone, self.itehelptwo, self._dictomoles[self.itehelpone][self.itehelptwo]) #tuple with two keys

    def __len__(self):
        return len(self._orderlist)


    def __str__(self):
        tab = []
        for molinfo in self._orderlist:
            for line in self._dictomoles[molinfo[0]][molinfo[1]].selftolist():
                tab.append(line)
        return ''.join(tab)
        
    def __getitem__(self, ind):
        if isinstance( ind, slice ) :
            mols = self._orderlist[ind]
            return [self._dictomoles[inf[0]][inf[1]] for inf in mols]
        elif isinstance( ind, int ) :
            mols = self._orderlist[ind]
            return self._dictomoles[mols[0]][mols[1]]
        else:
            raise TypeError, "Invalid argument type."
        
        

    def sdfseparator(self, strings):
        #Hajoittaa .sdf tiedoston listaksi listoja, joista yksi sisaltaa yhden tiedoston molekyyleista
        _moles = []
        _lines = []
        for line in strings:
            if molchop.match(line):
                _moles.append(_lines)
                _lines = []
            else:
                _lines.append(line)
        del(_lines)
        return _moles
    
    def add(self, stringsofone):
        #Add a molecule
        new = Sdfmole(stringsofone)
        name = new.getname()

        if not name in self._dictomoles:
            self._dictomoles[name]=dict()
        n = self.tempconfn(new, self._dictomoles[name])
            
        self._dictomoles[name][n] = new
        self._orderlist.append([name,n])
        
    def remove(self, name, confn):
        #Remove a molecule
        del(self._dictomoles[name][str(confn)])
        self._orderlist.remove([name,confn])

                        
    def sdfmetacombi(self, other, bolist=[True, True], overwrite=False):
        #Combine metadata from another sdf-file.
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
        #Remove that are present in the other file. Cut operator.
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
        '''
        levels = leveler(metalogic)
        for mol in self:
            mol.stripbutmeta(levels)
        
    
    #conformation numbers
    
    def tempconfn(self, mol, dictomol=None):
        #Returns a temporary conformation number
        n = mol.getconfn()
        if not n:
            i=-1
            while dictomol:
                j=str(i)
                if not j in dictomol:
                    return j
                else:
                    i -= 1
            n = str(i)
        return n
        
    def uniconfn(self, mol, n='-1'):
        #Returns a unique conformation number
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
        for i, molinfo in enumerate(self._orderlist):
            mol = self._dictomoles[molinfo[0]][molinfo[1]]
            c = mol.getconfn()
            if not c:
                c = self.uniconfn(mol, 1)
                #c = self.uniconfn(mol, molinfo[1])
            #self._dictomoles[molinfo[0]][molinfo[1]].
            mol.addconf(c, bolist)
            del(self._dictomoles[molinfo[0]][molinfo[1]])
            self._dictomoles[molinfo[0]][c]=mol
            self._orderlist[i]=[molinfo[0],c]
            
    def remconfs(self, bolist = [True, True]):
        for mol in self:#self._dictomoles:
            mol.remconf(bolist)
            
    #misc

    def metatoname(self, meta, joiner='_'):
        '''
        Changes the name of molecules to whatever found in given metafield.
        Applies only for metafields of type str
        '''
        #for i, molord in enumerate(self._orderlist):
        for molord in self._orderlist:
            olname = molord[0]
            n = molord[1]
            mol = self._dictomoles[olname][n]
            
            name = joiner.join(mol._meta[meta].getmetastrings())
            #index = self._orderlist.index([olname,n])
            #self._orderlist[]=
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
        Method to merge multiple metafields into one. Comes with a few operations. Frontend to mergenewmeta
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
                    #newmeta.extend(mol._meta[meta])
                    tab.append(mol._meta[meta]) #add meta instead of number
                else:
                    try:
                        num =  mol.logicgetmeta(leveler(meta)) #numify(meta)
                    except ValueError:
                        #meta isn't a statement
                        continue
                    if type(num) == Sdfmeta:
                        tab.append(num)
                    #if type(num) in (float, int):
                        #tab.append(Sdfmeta.construct(num))
                        #newmeta.extend(Sdfmeta.construct('',num)) #You can use numbers instead of metafield too.
            #mol.addmeta(newmetaname, Sdfmeta.metaoper(oper, tab))
            newmet = Sdfmeta.metaoper(oper, tab)
            #if not newmet:
                #for met in tab:
                #    print met._datatype
                #print tab
                #pass
            newmet.cleandelim(True)
            #newmet.setname()
            mol.addmeta(newmetaname, newmet)
            #mol._meta[]
            
            '''
            try:
                mol.addmeta(newmeta, [str(oper(tab))])
            except TypeError:
                mol.addmeta(newmeta, [str(oper(map(str,tab)))])
            '''
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
                #break
                return
        self.makenewmeta(name, ostri.strip(), None, None )
        
    
    def makenewmeta(self, name, metastatement, logicchar = None, value = None):
        comps = {'>=':operator.ge, '<=':operator.le, '<':operator.lt, '>':operator.gt, '==':operator.eq, '=':operator.eq, '!=':operator.ne }
        level = leveler(metastatement)
        if logicchar:
            for mol in self:
                newmeta = copy.deepcopy( mol.logicgetmeta(level) )  #\FIXME
                newmeta.pickvalues( numify(value), comps[logicchar] )
                if len(newmeta)>0:
                    mol.addmeta(name, newmeta)
        else:
            for mol in self:
                newmeta = copy.deepcopy( mol.logicgetmeta(level) )  #\FIXME
                #newmeta.pickvalues( numify(value), comps[logicchar] )
                #if len(newmeta)>0:
                mol.addmeta(name, newmeta)
    
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
    
    def sametest(self, bolist, samelist):
        #Test wether molecules are "the same" as described in bolist.
        return (bolist[0] <= samelist[0]) and (bolist[1] <= samelist[1])
        
    def makecsv(self,stringofmetas,separator=';\t'):
        #Make a csv-list containing all molecules and given metafields as columns. '?' gives all fields
        listofmeta = [met.strip() for met in re.split('\s*,|;\s*',stringofmetas)]
        if listofmeta[0]=='?':
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
                    line.append('"'+mol._meta[meta].getmetastr()+'"')
                else:
                    line.append('""')
            csv.append(separator.join(line))
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
            counts.append('{}\t{}'.format(mol,len(self._dictomoles[mol])))
        return counts
    
    #should work
    def closest(self, point, **varargdict): #name, intrests, num
        #coord1 = numpy.array(coorder(point)) #map(numify, goodchop.split(point.strip('{[()]}')))
        #coord1 = numpy.array(coord1)
        if not 'name' in varargdict:
            name = ''
        else:
            name = varargdict['name']+'_'
        for mol in self:
            #is coord1 meta or not? If it was a metafield with atomnumber, it would give that one as closest, so no.
            #alist = mol.dists(coord1)
            alist = mol.dists(point)
            alist=sorted(alist, key=lambda item: item[0])
            '''
            mol.addmeta(name+'Closest_atom',str(alist[0][1]))
            mol.addmeta(name+'Closest_distance',str(alist[0][0]))
            '''
            newlist = []
            moles=[]
            if 'intrests' in varargdict:
                switch = {list:varargdict['intrests'],int:[varargdict['intrests']],str:mol.logicgetmeta(leveler(varargdict['intrests']))._data}
                '''
                if type(varargdict['intrests']) == list:
                    mylist = varargdict['intrests']
                elif type(varargdict['intrests']) == int:
                    mylist = [varargdict['intrests']]
                elif type(varargdict['intrests']) == str:
                    mylist = [varargdict['intrests']]
                '''
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
            

    '''old version
    def closer(self, point, meta, name=None):
        coord1 = coorder(point) #map(numify, goodchop.split(point.strip('{[()]}')))
        if type(coord1)==list:
            coord1 = numpy.array(coord1)
        else:
            #raise error?
            #metar = coord1
            pass
        for mol in self:
            #is coord1 meta or not?
            alist = mol.dists(coord1)
            alist=sorted(alist, key=lambda i: i[0])
            for i in range(len(alist)):
                alist[i].append(i)
            #alist=sorted(alist, key=lambda i: i[2])
            try:
                intrests = map(numify, re.split('\s+',mol._meta[meta.strip()]))
            except TypeError:
                continue
            closers=[]
            for atomn in intrests:
                closers.append(alist[nestfind(alist, 1, atomn)])
            if len(closers)>0:
                closers = sorted(closers, key=lambda i: i[0])
                mol.addmeta('Closest_atom_from_{}'.format(meta),str(closers[0][1]))
                mol.addmeta('Closer_atoms_than_{}'.format(meta),str(closers[0][2]))
    '''
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
        dot = coorder(point)#self.coordormeta(point)
        if type(dot) != list:
            #print dot
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
                #print 'point: '+point
                #print 'metaf: '+metafield
                #print 'dot:   '+str(dot)
                continue
            mol.addmeta(metafield+'_closest_distances',' '.join([str(cell) for cell in mind]))
            ind = mind.index(min(mind))
            mol.addmeta(metafield+'_closest_atom',str(mins[ind]))
            mol.addmeta(metafield+'_closest_distance',str(mind[ind]))

    def closestbybonds(self,fromwfield,towfield,newfield):
        for mol in self:
            if towfield in mol._meta and fromwfield in mol._meta:
                home = mol._meta[fromwfield]
                target = mol._meta[towfield]
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
        #meta = True
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
        moles = self.sdfseparator(data)
        del (data)
        for mole in moles:
            self.add(mole)
        
        
    def selftostring(self, output, **kwargs):
        if output=='getcsv':
            return '\n'.join(self.makecsv(kwargs['csv']))+'\n'
        elif output=='metalist':
            return '\n'.join(self.listmetas())+'\n'
        elif output=='counts':
            towrite=[]
            if kwargs['counts']==0 or kwargs['counts']==2:
                towrite = ['Total {}\t{}'.format(len(self._dictomoles),len(self))]
            if kwargs['counts']==1 or kwargs['counts']==2:
                towrite.extend(self.counts())
            return '\n'.join(towrite)+'\n'
        elif output=='sdf':
            return str(self)
        else:
            return None

    def mehist(self, Xname, Yname=None, **kwargs):
        #Histograms work again!
        #kwargs: title, Xtitle, Ytitle, bins
        newargs = dict()
        titargs = ['Xtitle','Ytitle','title']
        for key in kwargs:
            if key in titargs:
                newargs[key]=kwargs[key]
                #del(kwargs[key])
        for key in newargs:
            del(kwargs[key])
        
        #New implementation
        datas = [[]]
        levels = [leveler(Xname)]
        if Yname:
            datas.append([])
            levels.append(leveler(Yname))
        #YData = None
        
        for mol in self:
            #Xcan = mol.logicgetmeta(Xlevel)
            metas = tuple( mol.logicgetmeta(level) for level in levels )

            try:
                minlen = min( {len(meta) for meta in metas}-{1} ) # shortest list length, ignore singles
                #structs = set([meta._datastruct for meta in metas])-{'single'}
                #ostru = iter(structs).next()
                #ostru = iter(set(tuple(meta._datastruct for meta in metas))-{'single'}).next() #the other structuretype
            except ValueError:
                #Only singles
                minlen = 1
                #structs = {list}
                #ostru = list
                #singles = False
            try:
                structs = set([meta._datastruct for meta in metas])-{'single'}
                ostru = iter(structs).next()
            except StopIteration:
                ostru = list
                
            
            
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
            #    print keyorder
                
            #if ostru == list:
            #make lists have the same length
            #print meta._data
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
            '''
            if len(structs)>1:
                for i,meta in enumerate(metas):
                    if meta._datastruct == OrDi:
                        meta = copy.deepcopy(meta)
                        meta._data = tstmeta._data.values()
                        metas[i] = meta
            elif iter(structs).next() == OrDi:
                keys = metas[0]._data.keys()[:minlen]
                for i, meta in enumerate(metas):
                    for key in keys:
                        try:
                            datas[i].append(meta._data[key])
                        except AttributeError:
                            datas[i].append(meta._data[0])
            else:
                for i, meta in enumerate(metas):
                    if len(meta)<minlen:
                        datas[i].extend(meta._data[:1]*minlen)
                    elif meta._datastruct==list:
                        datas[i].extend(meta._data[:minlen])
            '''
            
        pylab.figure()
        if not Yname:
            larg=[datas[0]] #[getcolumn(data,Xname)]
            if 'bins' in kwargs:
                larg.append(kwargs['bins'])
                #print kwargs['bins']
                del(kwargs['bins'])
            (n, bins, patches)=pylab.hist(*larg,**kwargs) #kwargs?
        else:
            X,xe,ye=pylab.histogram2d(datas[0],datas[1],**kwargs)
            ex=[xe[0],xe[-1],ye[-1],ye[1]]
            pylab.imshow(numpy.transpose(X),extent=ex,interpolation='nearest',aspect='auto')
            cb=pylab.colorbar()
            #cb.spacing
            
        ''' Old implementation
        data = csvtomatrix(self.selftostring('getcsv',csv='?').split('\n'))[:-2]
        pylab.figure()
        if not Yname:
            larg=[getcolumn(data,Xname)]
            if 'bins' in kwargs:
                larg.append(kwargs['bins'])
                #print kwargs['bins']
                del(kwargs['bins'])
            (n, bins, patches)=pylab.hist(*tuple(larg),**kwargs) #kwargs?
        else:
            X,xe,ye=pylab.histogram2d(getcolumn(data,Xname),getcolumn(data,Yname),**kwargs)
            ex=[xe[0],xe[-1],ye[-1],ye[1]]
            pylab.imshow(numpy.transpose(X),extent=ex,interpolation='nearest',aspect='auto')
            cb=pylab.colorbar()
            #cb.spacing
        '''
        
        if 'Xtitle' in newargs:
            pylab.xlabel(newargs['Xtitle'])
        else:
            pylab.xlabel(Xname)
        if 'Ytitle' in newargs:
            pylab.ylabel(newargs['Ytitle'])
        elif Yname:
            pylab.ylabel(Yname)
        if 'title' in newargs:
            pylab.title(newargs['title'])
        
    
    def show(self):
        pylab.show()
        
   
       
    def logicchop(self, string): #This is horrible. Don't even try if you value your sanity...
        exp = '<|>|%|=|~|!' 
        ope = re.findall(exp, string)
        #ope = mope.group(0)
        things = re.split(exp, string)
        foreach = None
       
        #print ope
        #print things

        if len(ope)==1:
            if ope[0]!='!':
                value = numify(things[1].strip())
        elif len(ope)==2:
            foreach = dict()
            for mol in self:
                name = mol.getname()
                if not name in foreach:
                    foreach[name] = list()
                foreach[name].append(numify(mol._meta[things[0]]))
            if ope[1] =='%':
                if things[2].upper() == 'T':
                    fu = max
                elif things[2].upper() == 'B':
                    fu = min
                elif things[2].upper() == 'A':
                    fu = avg
                elif re.match('M|N',things[2].upper()): #things[2].upper() == 'M':
                    if ope[0] == '<':
                        reve=False
                    elif ope[0] == '>':
                        reve=True
                    for key in foreach:
                        values = foreach[key]#[numify(mol._meta[things[0]]) for mol in foreach[key]]
                        if things[2].upper() == 'M':
                            N = int(ceil(len(values)*float(things[1])/100.0))
                        else:
                            N = int(things[1])
                        try:
                            foreach[key] = sorted(values,reverse=reve)[N]
                        except IndexError:
                            if len(values)<2:
                                if reve:
                                    foreach[key] = values[0]-1
                                else:
                                    foreach[key] = values[0]+1
                            else:
                                foreach[key] = sorted(values,reverse=reve)[-1]
                            #if args.verbose:  #TODO
                            #    print 'Too few conformations for {}'.format(key)
                else:
                    print 'Invalid statement. End with [T]op, [B]ottom or [A]verage'
            try:
                for key in foreach:
                    foreach[key]=fu(foreach[key])*(numify(things[1])/100)
            except UnboundLocalError:
                pass
        deles = []
        for i, inf in reversed(list(enumerate(self._orderlist))):
            try:
                meta = numify(self._dictomoles[inf[0]][inf[1]]._meta[things[0].strip()].strip())
                if foreach:
                    value = foreach[inf[0]]
                if ope[0] == '<':
                    if not meta < value:
                        deles.append(inf)
                elif ope[0] == '~':
                    if not re.search('(^|\s|,|;)'+re.escape(self._dictomoles[inf[0]][inf[1]]._meta[things[1]])+'($|\s|,|;)', str(meta)):
                        deles.append(inf)
                elif ope[0] == '=':
                    if not re.search('(^|\s|,|;)'+re.escape(str(things[1]).strip())+'($|\s|,|;)', str(meta)):
                        deles.append(inf)
                elif ope[0] == '>':
                    if not meta > value:
                        deles.append(inf)
                elif ope[0] == '!':
                    #l=len(self)
                    for j in range(0,i):#leaves the first one intact #for j in range(i+1,l): #Leaves the last one in place!
                        meinf = self._orderlist[j]
                        if inf[0]!=meinf[0]:
                            continue
                        if meta == numify(self._dictomoles[meinf[0]][meinf[1]]._meta[things[0].strip()].strip()):
                            deles.append(inf)
                            break
            except KeyError:
                deles.append(inf)
        for key in deles:
            try:
                self.remove(key[0],key[1])
            except KeyError:
                print key
                continue
        self.dictmaint()
        
    
    def logicparse(self, string):
        
        def dealer(picks, drops):
            #if pick:
            for info in drops:
                del(self._dictomoles[info[0]][info[1]])
            self._orderlist = picks
            self.dictmaint()
            '''
            else:
                for info in trues:
                    del(self._dictomoles[info[0]][info[1]])
                self._orderlist = falses
                self.dictmaint()
            '''
        string = string.strip()
        if string[:2] == '+ ':
            pick = True
            string = string[2:]
        elif string[:2] == '- ':
            pick = False
            string = string[2:]
        else:
            pick = True
        
        #pick = True
        opesplit=[item.strip() for item in re.split('(>=|<=|<|>|==|=|!=)',string)]
        comps = {'>=':operator.ge, '<=':operator.le, '<':operator.lt, '>':operator.gt, '==':operator.eq, '=':operator.eq, '!=':operator.ne }
        
        if len(opesplit)==3:
            opera = opesplit[1]
            tocompare = [ leveler(opesplit[0]) , leveler(opesplit[2]) ]
            trues = []
            falses = []
            for info in self._orderlist:
                mole = self._dictomoles[info[0]][info[1]]
                if comps[opera]( mole.logicgetmeta(tocompare[0]), mole.logicgetmeta(tocompare[1]) ):
                    trues.append(info)
                else:
                    falses.append(info)
                
                
        elif len(opesplit)==1:
            tear = leveler(string)
            if len(tear) != 2:
                raise Error('Weird logic.')
            #funks = {'min':dicmin, 'max':dicmax}
            funk = tear[0]
            matheus = tear[1][1]
            
            rcomindex = matheus[-1].rfind(',')
            if rcomindex == -1:
                raise TypeError('Weird logic. No comma.')
            
            if len(matheus) > 1:
                metatab = matheus[:-1]
            else:
                metatab = matheus[-1][:rcomindex]
            
            #metastring = tear[1][1][-1][:rcomindex]
            #numstring = tear[1][1][-1][rcomindex+1:]
            numstring = matheus[-1][rcomindex+1:]
            
            perindex = numstring.rfind('%')
            if perindex > 0:
                num = numify(numstring[:perindex])
                per = True
            else:
                num = numify(numstring)
                per = False
            trues = []
            falses = []
            for molec in self._dictomoles:
                moles = OrDi()
                for conf in self._dictomoles[molec]:
                    moles[conf] = self._dictomoles[molec][conf].logicgetmeta(metatab)
                    #metas.append( self._dictomoles[molec][conf].logicgetmeta(tear[1][1]) )
                if tear[0]=='max':
                    reverse = True
                else:
                    reverse = False
                moles = OrDi(sorted(moles.iteritems(), key= lambda xx: xx[1], reverse = reverse))
                if per:
                    grab = int(math.ceil(len(moles)*num/100.0))
                    #if grab == 0:
                    #    grab = 1
                else:
                    grab = int(num)
                trues.extend( [[molec, item] for item in moles.keys()[:grab]] )
                falses.extend( [[molec, item] for item in moles.keys()[grab:]] )
            
        if pick:
            dealer(trues, falses)
        else:
            dealer(falses, trues)
            #find min or max
        
        

    def sorter(self,sortstring):
        '''
        if datastructure has more than 1 entry, picks the 1st.
        '''
        sortstring=sortstring.strip()
        sortstring=sortstring.strip('\"|\'')
        if sortstring[0]=='<':
            rever=True
        elif sortstring[0]=='>':
            rever=False
        else:
            #if args.verbose:
            #    print 'Bad sortstring'
            return self._orderlist
        return sorted(self._orderlist, key=lambda avain: numify(self._dictomoles[avain[0]][avain[1]]._meta[sortstring[1:]][0]),reverse=rever)
        
    
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
                mol._meta[tosort].sortme(ascending,byValue)
        #return sorted(self._orderlist, key=lambda avain: numify(self._dictomoles[avain[0]][avain[1]]._meta[sorty[1:]][0]),reverse=rever)
    
    def addcsvmeta(self, path, verbose=False):
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
                pass
                #if verbose:
                #    print seeker+' not found!'
                #continue
            if m:
                keys = [n]
            else:
                keys = self._dictomoles[seeker].keys()
            for key in keys:
                if key in self._dictomoles[seeker]:
                    mol = self._dictomoles[seeker][key]
                else:
                    continue
                for i, newmeta in enumerate(csvmol[1:]):
                    if len(newmeta)!=0:
                        mol.addmeta(header[i+1], newmeta, overwrite=True)
                        mol._meta[header[i+1]].cleandelim(True)
        del(csvdata)
                    
    def removemeta(self, metaliststring, pick = False):
        remoli = re.split('\s*,|;\s*',metaliststring)
        if pick:    #keep selected metafields
            '''if len(metaliststring)==0:
                #remoli =['Nonspecifiedfield']
                for mol in self:
                    for thing in mol._meta:
                        del(mol._meta[thing])
                        mol._metakeys.remove(thing)
            else:'''
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
            #noff=ceil(nofm/n)
        elif n>0:
            noff=n
            n=int(math.ceil(float(nofm)/noff))
            #if nofm%noff!=0:
            #    n+=1
        files = [Sdffile()]
        #filen=0
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
                #del(newargs['path'])
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
                        #onepath = kwargs['path'][:dotplace]+'_'+str(i)+'/'+kwargs['path']
                        onepath = kwargs['makefolder']+'_'+str(i)+'/'+kwargs['path']
                        ensure_dir(onepath)
                    newargs['path'] = onepath
                    onesdf.writer(writetype, **newargs) #TODO
                    
    def getMol2Data(self, metaname, column, path):
        mol2 = Mol2File(path)
        for i, mol in enumerate(mol2):
            newdic = mol.pickatomdata(column)
            self[i].addmeta(metaname, newdic)
        del(mol2)
        
    def injectMol2Data(self, metaname, column, path, defaultValue=0.0, precision=4, outpath = None):
        mol2 = Mol2File(path)
        metaname = metaname.strip()
        for i, mol in enumerate(mol2):
            mol.injectatomdata(self[i]._meta[metaname], column, defaultValue, precision)
        if not outpath:
            outpath=path
        mol2.writefile(outpath)
        
#end of Sdffile

class Sdfmole(object):
    mes='Conformation number mismatch! {} vs {}'
    
    def __init__(self, stringsofone=None):
        self._name = '' #nimi, ei ID
        self._meta = dict() #sisältää metadatan
        self._metakeys = list() #sisältää metadatan avaimet ja pitää niille järjestyksen
        self._other = list() #sisältää kaiken nimen ja ekan metadatan välissä #myöhemmin ainoastaan headerin nimen ja atomien välissä
        self._numeric = False
        self._ignoretype = []
        
        if stringsofone:
            self.initialize(stringsofone)
            
    def __str__(self):
        return ''.join(selftolist())
    
    def __copy__(self):
        new = Sdfmole()
        new._name           = copy.copy(self._name)
        new._meta           = copy.copy(dicself._meta)
        new._metakeys       = copy.copy(self._metakeys)
        new._numeric        = copy.copy(self._numeric)
        new._ignoretype     = copy.copy(self._ignoretype)
        new._comment        = copy.copy(self._comment)
        if new.numeric:
            #new._comments   = copy.copy(self._comments)
            new._counts     = copy.copy(self._counts)
            new._atoms      = copy.copy(self._atoms)
            new._bonds      = copy.copy(self._bonds)
            new._properties = copy.copy(self._properties)
        else:
            new._other      = copy.copy(liself._other)
        return new
        
    def __deepcopy__(self,memo):
        new = Sdfmole()
        new._name           = copy.deepcopy(self._name,memo)
        new._meta           = copy.deepcopy(self._meta,memo)
        new._metakeys       = copy.deepcopy(self._metakeys,memo)
        new._numeric        = copy.deepcopy(self._numeric,memo)
        new._ignoretype     = copy.deepcopy(self._ignoretype,memo)
        new._comment        = copy.deepcopy(self._comment)
        if new._numeric:
            #new._comments   = copy.deepcopy(self._comments)
            new._counts     = copy.deepcopy(self._counts,memo)
            new._atoms      = copy.deepcopy(self._atoms,memo)
            new._bonds      = copy.deepcopy(self._bonds,memo)
            new._properties = copy.deepcopy(self._properties,memo)
        else:
            new._other      = copy.deepcopy(self._other,memo)
        return new
    
    ''' deprecated
    def metahand(self, strilist):
            return strilist[0]
            
    '''
    
    def initialize(self, strings):
        self._name = strings[0].strip()
        self._comment = [strings[i].strip() for i in [1,2]]
        first = -1
        
        for i, line in enumerate(strings[1:]): #notice index shift: i=0; line=strings[1]
            if metaname.match(line):#if line[:3] == '> <':
                first = i+1
                break
            else:
                self._other.append(line.rstrip('\n'))
                #self._other.append(line.strip())
                
        #key = ''
        #sto = []
        #NEW 2014.03.14 Changing metadata back to multiline. Also reject lines longer than 74 characters.
        '''
        for line in strings[first:]: #some optimazition required in the order of statements
            m=metachop.match(line)
            if m:
                if key != '':
                    self._meta[key] = self.metahand(sto)
                    self._metakeys.append(key)
                    sto = []
                tmp = line.strip()
                key = tmp[len(m.group(0)):tmp.rfind('>')]
            else:
                sto.append(line.strip())
        if key != '':
            self._meta[key] = self.metahand(sto)
            self._metakeys.append(key)
            sto = []
        '''
        #print str(first)
        add = 0
        for i, line in enumerate(strings[first:]): #some optimazition required in the order of statements
            if line == '\n': #emptyline
                if i-add < 2:
                    add = i
                    #print 'ping! empty line on line {}'.format(i)
                    continue
                newmeta = Sdfmeta(strings[first+add:first+i+1])
                add = i+1
                self._meta[newmeta.getname()] = newmeta
                self._metakeys.append(newmeta.getname())
                
    def dists(self, point1, ignores=['H']):
        ignores = [item.upper() for item in ignores]
        self.numerize()
        coord1 = coorder(point1)
        if coord1 == None:
            #coord1 =  self.getatomloc(self._meta[point1]._data[index])
            coord1 =  self.getatomloc(self.logicgetmeta(leveler(point1+'[0]'))._data[0])
        alist=[]
        for i in range(len(self._atoms)):
            #this responsibility will be taken elsewhere
            if self.gettype(i+1).strip() in ignores:
                continue
            
            coord2 = numpy.array(self.getcoord(i+1)) #numpy.array([float(c) for c in atom[0:3]])
            dist = sum((coord1-coord2)**2)**0.5
            alist.append([float(dist),i+1])
        #alist=sorted(alist, key=lambda i: i[0])
        return alist

    def numerize(self):
        if not self._numeric:
            #self._comments = list()
            self._counts = [self._other[2][i:i+3].strip() for i in range(33)[::3]]+[self._other[2][33:39].strip()]
            blockdiv = [3,3+int(self._counts[0]),3+int(self._counts[0])+int(self._counts[1])]
            self._atoms = [[atom[atomcut[i]:atomcut[i+1]] for i in range(len(atomcut)-1)] for atom in self._other[3:blockdiv[1]]]
            self._bonds = [[bond[i:i+3].strip() for i in range(7*3)[::3]] for bond in self._other[blockdiv[1]:blockdiv[2]]]
            self._properties = self._other[blockdiv[2]:]
            del(self._other) #DANGER!!!
            self._numeric = True
        
    def getconf(self):
        ans = [None, None]
        m = confchop.search(self._name)
        if m:
            ans[0] = m.group(0)[2:-2]
        if ckey in self._meta:
            #ans[1] = str(self._meta[ckey][0]).strip()
            ans[1] = self._meta[ckey][0]
        if ans[0] and ans[1]:
            if ans[0]!=str(ans[1]):
                print(self.mes.format(ans[0],ans[1]))
        return ans
        
    def getconfn(self):
        a = self.getconf()
        return a[0] or a[1]
        #return next((item for item in a if item is not None), None)
        
    def addconf(self, conf, bolist = [True, True]):
        conf = str(conf)
        already =  self.getconf() 
        if bolist[1]:
            if already[1]:
                if conf != str(already[1]):
                    print(self.mes.format(already[0],already[1]))
            else:
                #self._meta[ckey] = [conf]
                self.addmeta(ckey, conf, literal = True)
        if bolist[0]:
            if already[0]:
                if already[0] != conf:
                    print(mes.format(ans[0],ans[1]))
            else:
                self._name = self._name + '{[' + conf + ']}'
                
    def remconf(self, bolist = [True, True]):
        already = self.getconf()
        if bolist[0] and already[0]:
            self._name = confchop.sub('',self._name)
        if bolist[1] and already[1]:
            self._metakeys.remove(ckey)
            self._meta.pop(ckey)
            
    def getname(self):
        return confchop.sub('', self._name)
    
    def getcoord(self, N):
        return numify(self._atoms[N-1][:3])
    
    def gettype(self, N):
        return self._atoms[N-1][4]

    def metacombine(self, othersdfmol,overwrite=False):
        for key in othersdfmol._metakeys:
            if not overwrite and key in self._meta:
                continue
            self.addmeta(key,othersdfmol._meta[key],overwrite=True)
            #self._metakeys.append(key)
            #self._meta[key]=list(othersdfmol._meta[key])
                
    def issame(self, othersdfmol):
        same = [False, False]
        if confchop.sub('', self._name) == confchop.sub('', othersdfmol._name):
            same[0]=True
        if self.getconfn() == othersdfmol.getconfn():
            same[1]=True
        return same
    
    def make_other(self):
        other = []
        for line in self._comment:
            other.append(line + '\n')
        other.append(listtostring(self._counts[:-1],3)+'{:>6}'.format(self._counts[-1])+'\n') # atomheader!
        for atom in self._atoms:
            other.append(atomtostring(atom)+'\n')
        for bond in self._bonds:
            other.append(listtostring(bond, 3)+'\n')
        for prop in self._properties:
            other.append(prop+'\n')
        return other
    
    def selftolist(self):
        me=[]
        me.append(self._name + '\n')
        
        if self._numeric:
            me.extend(self.make_other())
        else:
            for line in self._other:
                me.append(line + '\n')
        
        ''' #Metastructure has changed
        for key in self._metakeys:
            me.append('>  <'+key+'>\n')
            #for line in self._meta[key]:
            #    me.append(line + '\n')
            me.append(self._meta[key] + '\n')
            me.append('\n')
        '''
        for key in self._metakeys:
            mystrings = [line+'\n' for line in self._meta[key].selftolistofstrings()]
            me.extend(mystrings)
        me.append('$$$$\n')
        return me
        
    #def addmeta(self,metafield, value, overwrite = True):
    def addmeta(self,metafield, value, **dictarg): #overwrite, literal
        if 'overwrite' in dictarg and dictarg['overwrite']:
            overwrite = True
            del(dictarg['overwrite'])
        else:
            overwrite = False
        
        if metafield in self._meta and not overwrite:
            return
        '''
        if type(value) == list:
            #value = [value]
            value = ' '.join(value)
        '''
        if isinstance(value, Sdfmeta): #type(value) == Sdfmeta:
            insert = value
        else:
            insert = Sdfmeta.construct(value, name = metafield, **dictarg)
        #insert.cleandelim(True)
        #insert.setname(metafield)
        self._meta[metafield] = insert
        self._meta[metafield].cleandelim(True)
        self._meta[metafield].setname(metafield)
        if not metafield in self._metakeys:
            self._metakeys.append(metafield)
        #print metafield+' '+self._meta[metafield]
        
    def metasort(self, reverse=False):
        self._metakeys.sort(reverse, key=lambda meta: meta.lower())
        
    def changemetaname(self, oldname, newname):
        if oldname in self._meta:
            i=self._metakeys.index(oldname)
            self._metakeys[i]=newname
            self._meta[newname]=self._meta[oldname]
            self._meta[newname].setname(newname)
            del(self._meta[oldname])
            
    def atomlistdistances(self,metaatom,metalist):
        ''' oneline version
        self.numerize()
        try:
            atom=numify(self._meta[metaatom])
            manyatom=self._meta[metalist]
        except KeyError:
            return None
        atoms = [numify(one) for one in re.split('\s+',manyatom.strip())]
        dists = []
        for oneatom in atoms:
            dists.append(self.atomdist(atom,oneatom))
        return (atoms,dists)
        '''
        '''
        used to give distances from one atom in metafield to many atoms in metafield
        in future will:
        add new meta with 
        '''
        self.numerize()
        try:
            atom=numify(self._meta[metaatom])
            manyatom=self._meta[metalist]
        except KeyError:
            return None
        atoms = [numify(one) for one in re.split('\s+',manyatom.strip())]
        dists = []
        for oneatom in atoms:
            dists.append(self.atomdist(atom,oneatom))
        return (atoms,dists)
        
    def getatomloc(self,n):
        self.numerize()
        return numpy.array([numify(coord) for coord in self._atoms[n-1][0:3]])
    
    def atomdist(self, atom1, atom2):
        self.numerize()
        return sum((self.getatomloc(atom1)-self.getatomloc(atom2))**2)**0.5
    
    def pointlistdists(self, point, metalist):
        #arpoint=numpy.array(point)
        #arpoint=coorder(point)
        #self.numerize()
        alldists = self.dists(point)
        try:
            atoms = self._meta[metalist]
        except KeyError:
            return None
        #atoms = [numify(one) for one in re.split('\s+',manyatom.strip())]
        dists = []
        for oneatom in atoms:
            #dists.append(self.disttopoint(oneatom,arpoint))
            dists.append(alldists[oneatom-1])
        return (atoms,dists)
    
    '''deprecated?
    def disttopoint(self, atom1, point):
        return sum((self.getatomloc(atom1)-point)**2)**0.5
    '''
    
    def bonddists(self, atom1, intrests=None):
        #dijikstra algorithm for finding path length to all atoms, or to list of intresting atoms
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
                #if distances[target-1]==numpy.inf:
                #    distances[target-1]=1
                if distances[target-1]>distances[atom-1]+1:
                    distances[target-1]=distances[atom-1]+1
            visited[atom-1]=True
            
            if intrests:
                if atom in intrests:
                    intrests.remove(atom)
                    if len(intrests)==0:
                        #print distances
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
        
    def stringtoownmeta(self,string):
        #find parens, test if same
        leveled = leveler(string)
        pass
    
    def logicgetmeta(self, partab):
        return self.collapser(self.levopemap(partab))
        
    
    def levopemap(self, tab, par=None, length=None):
        '''
        Accepts lists made by leveler method. Maps this list for mathematical operators and metafield names. Also slices.
        '''
        if isinstance(tab, tuple):
            return (tab[0], self.levopemap(tab[1], tab[0]))
        elif isinstance(tab, list):
            for j, thing in enumerate(tab):
                if not isinstance(thing, str):
                    continue
                thing = thing.strip()
                #print thing
                '''old
                sear = re.search('[+-]', thing)
                if not sear:
                    sear = re.search('[\*/]', thing)
                if sear:
                    tab1=[]
                    tab1.extend(tab[:i])
                    if sear.start() > 0:
                        tab1.append(thing[:sear.start()])
                    if sear.end() < len(thing):
                        tab2 = [thing[sear.end():]]
                    else:
                        tab2 = []
                    tab2.extend(tab[i+1:])
                    return (sear.group(),self.levopemap(tab1),self.levopemap(tab2))
                '''
                
                ''' tryout
                def opeser(whats,fromwhat):
                    if not type(whats) in (list, tuple):
                        what = (what)
                    bi = -1
                    Lfw = len(fromwhat)
                    for what in whats:
                        if type(what) != str:
                            raise TypeError
                        Lw = len(what)
                        dup = 2*what
                        for i in range(Lfw)
                '''
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
                    #return (sear.group(),self.levopemap(tab1),self.levopemap(tab2))
                    return (thing[i:ie],self.levopemap(tab1),self.levopemap(tab2))
                #new ends
            return [self.levopemap(item, par,len(tab)) for item in tab]
        else:
            
            tab = tab.strip()
            
            if tab in self._meta:
                return self._meta[tab]
            
            if length==1 and par!=None:
                return str(tab)
                ''' slicing happens later on :)
                if par in ('[','{'):
                    try:
                        return slice(*[{True: lambda n: None, False: int}[x == ''](x) for x in (tab.split(':') + ['', '', ''])[:3]])
                    except ValueError:
                        print tab
                        raise ValueError('Your logic makes no sense...')
                if par == '(':
                    return 
                '''
            else:
                tab = numify(tab)
                if not isinstance(tab, str):
                    return Sdfmeta.construct(tab)
                else:
                    #return tab
                    #print tab
                    raise ValueError('Your logic makes no sense: '+tab)
                        
    
    def collapser(self, tab, para=None):
        #print (tab, para)
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
                #multiple parentheses in series. First must be meta.
                metaus = self.collapser(tab[0])
                for item in tab[1:]:
                    #sli = self.collapser(item[1],item[0])
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
                #print tab
                raise ValueError('Your logic makes no sense...')
        elif isinstance(tab, Sdfmeta):
            return tab #OR THIS
        else:
            #print tab
            raise TypeError('Unknown logic')
            
    
    def stripbutmeta(self, levemeta):
        self.numerize()
        atoms = self.logicgetmeta(levemeta)
        newatoms = []
        for n in atoms._data:
            newatoms.append(self._atoms[n-1])
        self._atoms = newatoms
        self._bonds = []
        self._counts[0] = len(newatoms)
        self._counts[1] = 0

#end of Sdfmole

class Sdfmeta(object):
    #name
    #datatype: string list dict ?
    #the structure
    #                 |------ singles ? -----|
    #self._datypes = ['string', 'int', 'float' 'vlist', 'hlist', 'contlist', 'dict','unknown']
    #single = 'single'
    
    def __init__(self, listofstrings=None):
        
        self._name = [None,None] #Metafield name
        self._datatype = None #int, float, str
        self._datastruct = None #list, dict, single
        self._data = None #The actual data
        self._delims = []
        if listofstrings:
            self.initialize(listofstrings)
        
    def __getitem__(self, ind):
        if self._datastruct == OrDi:
            mylist = self._data.values()
        else:
            mylist = self._data
        if isinstance( ind, (slice, int) ):
            return mylist[ind]
        else:
            raise TypeError, "Invalid argument type."
            
    def __iter__(self):
        if type(self._data) == OrDi:
            return iter(self._data.values())
        else:
            return iter(self._data)
            
    def __len__(self):
        return len(self._data)
            
    #comparisons
    def __lt__(self, other):
        return self._compare(other,operator.lt)
        
    def __le__(self, other):
        return self._compare(other,operator.le)
        
    def __eq__(self, other):
        return self._eqcompare(other,operator.eq)
        
    def __ge__(self, other):
        return self._compare(other,operator.ge)
        
    def __gt__(self, other):
        return self._compare(other,operator.gt)
        
    def __ne__(self, other):
        return not self._eqcompare(other,operator.ne)
    
    def _compare(self, other, oper):
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
            #fulfils = True
            for item1 in li1:
                if oper(item1, item2):
                    #fulfils = False
                    #break
                    return True
            #if fulfils:
                #return True
        return False
        
    def initialize(self, listofstrings):
        '''parse the metadata'''
        #get name of metafield
        if listofstrings[0][0] != '>' :
            self._name = ['  ','']
            fi = 0
        else:
            self._name = metaname.match(listofstrings[0]).groups()
            fi = 1
            if re.match('[ ]{0,2}', self._name[0]):
                self._name = ('  ', self._name[1])
        #Remove the last empty line and linechanges
        if not len(listofstrings[-1].strip())==0:
            #Raise error?
            return
        else: #             rstrip?
            mylines = [line.strip('\n ') for line in listofstrings[fi:-1] ]
        
        #instead of commented section, merge lines with '' and do whattype
        #if no delimiter or ' ' in the end of line, add ' '
        newlines = []
        for line in mylines[:-1]:
            newlines.append(line)
            if not re.match('[ ,;\t]',line[-1]):
                newlines.append(' ')
        newlines.append(mylines[-1])
        (data,dtype,delims) = whattype(''.join(newlines))

        #if string, it's special
        if dtype == str:
            self._datatype = str
            self._data = mylines
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
    
    @staticmethod
    #def construct(name, data, delims=None):
    def construct(data, **dictarg): #name, delims, literal
        
        new = Sdfmeta()
        if 'name' in dictarg:
            name = dictarg['name']
            new.setname(name)
        else:
            name = ''
        '''
        if type(name) in (list, tuple):
            new._name = list(name[:2])
        elif type(name) == str:
            new._name = ['  ', name]
        else:
            raise TypeError('Name defined incorrectly!')
        '''
        if type(data) == str:
            
            if 'literal' in dictarg and dictarg['literal']:
                new._datastruct = 'single'
            else:
                (newdata, newtype, newdelims) = whattype(data)
                if newtype != str:
                    return Sdfmeta.construct(newdata, delims = newdelims, name = name)
            data = [data]
            #else:
            #    new._datastruct = 'single'
            #    data = [data]
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
            #print data
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
        new = Sdfmeta()
        new._name = copy.copy( self._name )
        new._datatype = copy.copy( self._datatype )
        new._datastruct = copy.copy( self._datastruct )
        new._data = copy.copy( self._data )
        new._delims = copy.copy( self._delims )
        return new
        
    def __deepcopy__(self, memo):
        new = Sdfmeta()
        new._name = copy.deepcopy( self._name,memo)
        new._datatype = copy.deepcopy( self._datatype,memo)
        new._datastruct = copy.deepcopy( self._datastruct,memo )
        new._data = copy.deepcopy( self._data,memo )
        new._delims = copy.deepcopy( self._delims,memo )
        return new
    
    def getname(self):
        return self._name[-1]
        
    def setname(self, newname):
        if type(newname) in (list, tuple):
            self._name = list(name[:2])
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
            #self._delims = [', ']*(len(self._data)-1)
            self._delims = [', ']*(len(self._data))
            self._delims.extend(other._delims)
            
        if floatflag:
            if self._datastruct == list:
                self._data = map(float, self._data)
            else:
                self._data = {k: float(v) for k, v in self._data.items()}
        
        
    def getmetastrings(self, length=float('inf')):
        dictflag = self._datastruct == OrDi
        if dictflag:
            key = self._data.keys()[0]
            tmp=[str(key)+':'+str(self._data[key])]
            del(key)
            itera = enumerate(self._data.keys()[1:])
        else:
            tmp=[str(self._data[0])]
            itera = enumerate(self._data[1:])
        l = len(tmp[-1])
        
        strings=[]
        for i, item in itera:
            if dictflag:
                stuff = str(item)+':'+str(self._data[item])
            else:
                stuff = str(item)
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
            if l +lda > length:
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
            self._delims = [', ']*(len(self._data)-1)
            return
        '''
        elif len(self._delims) < (len(self._data)-1):
            self._delims = [self._delims[0]]*(len(self._data)-1)
            self.cleandelim(unify)
            return
        '''
        
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
        workmetas = copy.deepcopy(metas)
        if oper in mathopers:
            if str in types:
                #cannot calculate strings, doh
                return None
            
            if len(set(structs))>2:
                raise TypeError('Mixed datastructs')
            
            
            #if len(set(structs))==2 and 'single' in structs:
            #singles and (list or OrDi) in structs, extend singles to lists
            '''Old implementation
            try:
                minlen = min( {len(meta._data) for meta in workmetas}-{1} ) # shortest list length, ignore singles
                ostru = iter(set(structs)-{'single'}).next() #the other structuretype
                singles = True
            except ValueError:
                #Only singles
                minlen = 1
                ostru = list
                singles = False
            '''
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
                return Sdfmeta.construct( listoper(oper, workmetas, singles) ) #Nameless meta
            elif ostru == OrDi:
                #extend singles to dicts
                keys = None
                keyorder = None
                #keys = set(workmetas[0]._data)
                #'''
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
                return Sdfmeta.construct(OrDioper(oper,workmetas))
        
        elif oper == metajoiner: #Therefore, it's joiner
            newmeta = metajoiner(metas)
            newmeta.cleandelim(True)
            return newmeta
            
        
        else:
            raise ValueError('Non-existent operator ')
        
        '''
        if oper != ' '.join:
            if len(types)>1:
                #quit if multiple types. Add support for  mixed structure, e.g., single and list
                return None
            
            if list in types:
                datas = [meta._data for meta in metas]
                minlen = min( map(len, datas) )
                datas = [data[:minlen] for data in datas]
        '''
    def pickvalues(self, value, operator):
        #MAY CREATE METAS WITH LENGT OF 0
        if type( self._data ) == OrDi:
            newdata = OrDi()
            for key in self._data:
                if operator(self._data[key], value):
                    newdata[key] = self._data[key]
        else:
            newdata = []
            for thing in self._data:
                if operator(thing, value):
                    newdata.append(thing)
        self._data = newdata
        self._delims = self._delims[:len(self._data)-1]
        
    def slicer(self, sliceorindex, paren):
        #parens = {'(':getpare,'[':getbrac,'{':getcubr}
        #slices = {Sdfmeta:pass, slice:pass}
        #parens = {'(':{Sdfmeta:pass, slice:pass},'[':{Sdfmeta:pass, slice:pass},'{':{Sdfmeta:pass, slice:pass}}
        
        def itsslice( toget, slic ):
            #print 'SLICE!'
            if paren == '(':
                return Sdfmeta.construct(  OrDi( [ (i, toget[i]) for i in self._data.keys()[slic] ] ) )
            #elif self._datastruct == 'dict'
            #    self._data.values()
            #
            else:
                #try:
                return Sdfmeta.construct(  toget[slic] )
                #except TypeError:
                #    print slic
                #    print self._data
            #return toget[sliceorindex]
        
        def itsmeta( toget, meta ):
            #print 'META!'
            #if meta._datastruct == 'dict':
            if type(meta._data) == OrDi:
                indexes = meta._data.values()
            else:
                indexes = meta._data
            if paren == '(':
                if type(self._data) != OrDi:
                    raise TypeError('"(" not applicaple for lists')
                return Sdfmeta.construct(  OrDi([(i, toget[i]) for i in indexes]) )
            return Sdfmeta.construct(  [toget[i] for i in indexes] )
        
        if paren == '{':
            toget = self._data.keys()
        #elif paren == '[' and self._datastruct == 'dict':
        elif paren == '[' and type(self._data) == OrDi:
            toget = self._data.values()
        else:
            toget = self._data
        
        if not type(toget) in (list, OrDi):
            raise TypeError('There is something wrong with our toget going to slicer') 
        slices = {Sdfmeta:itsmeta, slice:itsslice}
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
            #type == list
            self._data = sorted( self._data, reverse = not ascending )

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
                        #self._molecules.append(OrDi())
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
                stachar = sum( map(len, things[:mycol]) )
                endchar = stachar+len(things[mycol])
                injectinfo.append( (form.format(inject), stachar, endchar))
            else:
                continue
                #print line
        
        #print injectinfo
        
        sta = min( [item[1] for item in injectinfo] )
        end = max( [item[2] for item in injectinfo] )
        maxlen = max( [ len(item[0]) for item in injectinfo] )
        
        for i, line in enumerate(self['ATOM']):
            #if injectinfo[i][0]:
            #    putin = injectinfo[i][0]
            #else:
            #    putin = defaultValue
            #self['ATOM'][i] = line[:sta] + form.format(putin) + line[end+1:]
            self['ATOM'][i] = line[:sta] + (maxlen-len(injectinfo[i][0]))*' ' +  injectinfo[i][0] + line[end+1:]
        
#End of Mol2Mol

def coorder(point):
    coord = numpy.array( map(numify, goodchop.split(point.strip('{[()]}'))) )
    if type(coord) == numpy.ndarray and len(coord) > 1:
        return coord
    else:
        return None


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
    
def listtostring(tab, fill):
    return ''.join([('{:>'+str(fill)+'}').format(i) for i in tab])
    
def nestfind(alist, subindex, tofind):
    return next((i for i, sublist in  enumerate(alist) if sublist[subindex]==tofind),-1)
    
def numify(stri):
    if type(stri) == int:
        return stri
    elif type(stri) == float:
        #if int(stri)-stri == 0:
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
                raise TypeError('Wrong type')

def numitest(string):
    tststr = str(numify(string))
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
    tab= []
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
    #if len(splitted)>1:
    return (splitted,delim)
    #else:
    #    return None
        
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

'''
def stringtoownmeta(string):
    #find parens, test if same
    
    def leveler(pars, otab, string):
        print pars
        print string
        tab = copy.deepcopy(otab)
        if len(pars)<2:
            return string
        curpar = alpar[0]
        ci=0
        level = 1
        tab.append(string[:pars[0][1]])
        for i, item in enumerate(alpar[1:]):
            if curpar[0] in pair.values():
                ci=i+2
                curpar = pars[ci]
                continue
            elif item[0] == pair[curpar[0]]:
                level -= 1
                if level == 0:
                    #full parentheses here
                    tab.append([curpar[1],item[1],curpar[0],leveler(pars[ci+1:i+1],[],string[curpar[1]+1:item[1]])])
                    try:
                        ci = i+2
                        curpar = pars[ci]
                    except IndexError:
                        break
                    level = 1
                    continue
            elif item[0] == curpar[0]:
                level += 1
        endstring = string[item[1]+1:]
        if len(endstring)>0:
            tab.append(endstring)
        return tab
    
    alpar = [(item.group(), item.start()) for item in parre.finditer(string)]
    pair = {'(':')','[':']','{':'}'}
    tab=leveler(alpar,[],string)
    return tab
    '''
    
def leveler(string):
    pair = {'(':')','[':']','{':'}'}
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
        
    for i, item in enumerate(pars):
        
        if not curpar:
            curpar = item
            level = 1
            if backpar:
                addy = string[backpar[1]+1:curpar[1]]
                if len(addy) > 0:
                    tab.append(addy)
                backpar = None
        elif item[0] == curpar[0]:
            level += 1
        if item[0] == pair[curpar[0]]:
            level -= 1
            if level == 0:
                tab.append((curpar[0],leveler(string[curpar[1]+1:item[1]]))) #curpar[1],item[1],
                curpar = None
                backpar = item
    endstring = string[pars[-1][1]+1:]
    if len(endstring)>0:
        tab.append(endstring)
    return tab

'''
def levdepth(tab,i):
    if isinstance(tab, tuple):
        if len(tab)==2:
            return depth(tab[1],i+1)
    elif isinstance(tab, list):
        j=i
        mi=-1
        for ind, thing in enumerate(tab):
            k=depth(thing,i)
            if k>j:
                j=k
                mi=ind
        return j
    else:
        return i
'''



def tabjoin(taber):
    pair = {'(':')','[':']','{':'}'}
    string = []
    for item in taber:
        if isinstance(item, str):
            string.append(item)
        elif isinstance(item, list):
            string.extend([item[0],tabjoin(item[1]),pair[item[0]]])
    return ''.join(string)
    

def test(thing, i):
    '''
    Method for testing things returned by leveler and Sdfmole.levopemap
    '''
    if isinstance(thing, (list,tuple)):
        for item in thing:
            test(item, i+1)
    elif isinstance(thing, sdfconf.Sdfmeta):
        print ''.join(['  ']*i) + ' '.join(thing.selftolistofstrings())
    elif isinstance(thing, (str,slice)):
        print ''.join(['  ']*i) + str(thing)


    
def whattype(onestring):
    '''
    Tries to find out what type your string is
    Return (parsed,parsedtype,delim)
    '''
    #onestring = onestring.strip()
    #Number?
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
        #print 'onestring'
        #print onestring
        #Is it a list delimited by given delimiter ?
        splitted = listsep(onestring,delimiter)
        #is dict?
        dicti = OrDi() #OrderedDict
        notdict = False
        #print 'splitted'
        #print splitted
        for i, cell in enumerate(splitted[0]):
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
    
def singleoper(oper, metas):
    if len({meta._datastruct for meta in metas} - {single})>0:
        return None
    else:
        return [oper([meta._data[0] for meta in metas])]

if __name__ == "__main__":
    
    arger = argparse.ArgumentParser(description= 'Some bad-ass manipulation of SDF-files. Notice that documentation is not up to date.' )
    
    arger.add_argument("input", metavar = 'path', nargs='+', type = str,               help="Specify the input file")
    choicewrite = arger.add_mutually_exclusive_group()
    choicewrite.add_argument("-out", "--output", type = str, default=None,  help = "Specify output file.")
    choicewrite.add_argument("-o", "--overwrite", action='store_true',      help = "overwrite. you don't need to specify output file")
    arger.add_argument("-cf", "--tofield", action = "store_true",           help = "add conformation number to metadata from name. if number in name doesn't exist, make a new one.")
    arger.add_argument("-cn", "--toname",  action = "store_true",           help = "add conformation number to name from metadata. if number in metadata doesn't exist, make a new one.")
    arger.add_argument("-mtn", "--metatoname",  type = str,                 help = "Change the name of molecule to the data in given metafield")
    arger.add_argument("-ntm", "--nametometa",  type = str,                 help = "Copy the name of molecule into given metafield")
    arger.add_argument("-rc", "--remove",  type = int,                      help = "remove conformation number from 1=metadata, 2=name, 3=both. if number doesn't exist, do nothing")
    choicecombi = arger.add_mutually_exclusive_group()
    choicecombi.add_argument("-co", "--combine", type = str,                help = "Combine metadata from specified file to the data of original file. Confromationr_epik_Ionization_Penaltys must match.")
    choicecombi.add_argument("-aco", "--allcombine", type = str,            help = "Combine metadata from specified file to the data of original file. Names must match.")
    choicecut = arger.add_mutually_exclusive_group()
    choicecut.add_argument("-cu", "--cut", type = str,                      help = "Remove molecules in specified file from original file. Confromations must match.")
    choicecut.add_argument("-acu", "--allcut", type = str,                  help = "Remove molecules in specified file from original file. Names must match. Not tested")
    arger.add_argument("-csv", "--addcsv", type = str,                      help = "Add metadata from csv-file. File must have a 1-line header, it gives names to metafields. Names of molecules must be leftmost. If name includes confnumber, meta is only added molecules with same confnumber.")
    arger.add_argument("-ex", "--extract", type = str,                      help = "Pick or remove molecules from file by metafield info. Either with logical comparison or fraction of molecules with same name. Closest_atoms{:5}==soms, 2.5>Closest_atoms(soms)[], Closest_atoms[:3]<5.5. Takes multiple statements separated with | ")
    choiceremo = arger.add_mutually_exclusive_group()
    choiceremo.add_argument("-rm", "--removemeta", type = str,              help = "Remove metadata from molecules. Takes multiple values, separaterd by comma(,) or semicolon(;). If first is '?', means 'all but'")
    choiceremo.add_argument("-pm", "--pickmeta", type = str,                help = "Remove all nonspecified metadata from molecules. Takes multiple values, separaterd by comma(,) or semicolon(;). If first is '?', means 'all but'")
    arger.add_argument("-s", "--split", type = int,                         help = "Split the file into even pieces. Positibe means number of file while negative means number of molecules per file. 0 doesn't apply.")
    arger.add_argument("-mf", "--makefolder", action = "store_true",        help = "Put outputfile(s) into folder(s).")
    outputtype = arger.add_mutually_exclusive_group()
    outputtype.add_argument("-gc", "--getcsv", type = str ,                 help = "Writes a .csv-file istead of .sdf-file. Specify which fields you'll need, separated by ','.")
    outputtype.add_argument("-ml", "--metalist", action = "store_true",     help = "Writes a list of metafields.")
    outputtype.add_argument("-nm", "--counts", type = int,                  help = "Number of different molecules and different conformations. 0=just sums, 1=by molecule name, 2=both")
    outputtype.add_argument("-dnp", "--donotprint", action = "store_true",  help = "No output")

    arger.add_argument("-ca", "--closestatom", type = str,                  help = "Calculates the closest atoms (distances by atom number) to given point. Adds 'Closest_atoms' metafield with optional prefix. Needs either coordinates separated by ',' or or single atom number")
    #arger.add_argument("-mc", "--multiclosestatom", type = str,             help = "Calculates the distances from atom number in a given field, multiple atoms in second field. Creates a field <[secondfieldname]_closest_distance> ['Closest_atom,PRIMARY_SOM|Closest_atom,SECONDARY_SOM]")
    #arger.add_argument("-cb", "--closestbybonds", type = str,               help = "Calculates the closest atom of interest (atom number) from given atom (fieldname) to given atoms (fieldname). Adds '<input3>_atom' and '<input3>_distance' metafields. Needs three fieldnames, separated by ','. Fields are 'from', 'to' and 'newfield'. You can give multiple parameters, separate by '|'")
    arger.add_argument("-cla", "--closeratoms", type = str,                 help = "Calculates number of atoms closer to the given point, than the ones given adds metafields 'Closest_atom_from_{meta}' and 'Closer_atoms_than_{meta}'. Needs point and metafield name separated by ',', point first. Takes multiple parameters separated by '|'")
    arger.add_argument("-v", "--verbose", action = "store_true" ,           help = "More info on your run.")
    arger.add_argument("-mm", "--mergemeta", type = str,                    help = "Makes a new metafield based on old ones. newmeta=sum(meta1,meta2). operator are sum, max, min, avg, prod, div, power and join. Takes�multiple arguments separated by |")
    arger.add_argument("-mnm", "--makenewmeta", type = str,                 help = "Makes a new metafield based on logical statement and value picking inside the metafield. newmeta = meta1 + meta2 < 50. Takes�multiple arguments separated by |")
    arger.add_argument("-cm", "--changemeta", type = str,                   help = "Changes names of metafields. [olname1>newname1|oldname2>newname2]. ")
    arger.add_argument("-sm", "--sortmeta", type = str,                     help = "Sorts cells of a metafield in ascending [<metaname] or descending [>metaname] order. Additional '+' as the last character sorts by key in dictionary type metas. Takes multiple values separated by |")
    arger.add_argument("-so", "--sortorder", type = str,                    help = "Sorts molecules of a file in order of metafield. <MolecularWeight|>Id Sorts molecules first by highest weight first, then by smallest name first")
    arger.add_argument("-hg", "--histogram", type = str,                    help = "Plots a 1D or 2D histogram, multiple plots with '|'. 'Xname,Yname,Xtitle=x-akseli,Ytitle=y-akseli,bins=[30,30]'")
    
    arger.add_argument("-gm2", "--getmol2", type = str,                     help = "Reads atom block column data from mol2-file and adds it to sdf-file as metadata. pathto.mol2, column, metaname.")#meta column path
    arger.add_argument("-pm2", "--putmol2", type = str,                     help = "Injects meta-data from sdf-file and adds it to mol2-file as atom block column data. input.mol2,output.mol2, column, metaname, default, precision.")#metaname, column, path, defaultValue, precision, outpath
    
    arger.add_argument("-sbm", "--stripbutmeta", type = str,                help = "Removes all atoms from molecules, except for those in given logical statement. Takes multiple parameters separated by '|'")
    
    args = arger.parse_args()
    manyfiles = args.input #glob.glob(args.input)
    
    def main(inputfile):
        times = [time.time()]
        sdf1 = Sdffile(inputfile)
        times.append(time.time())
        
        if args.verbose:
            print 'Reading file done. It took {} seconds.'.format(times[-1]-times[-2])
        
        #check
        if args.tofield:
            sdf1.addconfs([False,True])
            times.append(time.time())
            if args.verbose:
                print 'Conformation numbers added to metadata. It took {} seconds.'.format(times[-1]-times[-2])
        
        #check
        if args.toname:
            sdf1.addconfs([True,False])
            times.append(time.time())
            if args.verbose:
                print 'Conformation numbers added to names. It took {} seconds.'.format(times[-1]-times[-2])
        
        #check
        if args.nametometa:
            sdf1.nametometa(args.nametometa)
            times.append(time.time())
            if args.verbose:
                print 'Name written to metafieldfield '+args.nametometa+'. It took {} seconds.'.format(times[-1]-times[-2])
        
        #check
        if args.remove==2:
            sdf1.remconfs([True, False])
            times.append(time.time())
            if args.verbose:
                print 'Conformation numbers removed from names. It took {} seconds.'.format(times[-1]-times[-2])
        
        #check
        elif args.remove==1:
            sdf1.remconfs([False, True])
            times.append(time.time())
            if args.verbose:
                print 'Conformation numbers removed from metafields. It took {} seconds.'.format(times[-1]-times[-2])
        
        #check
        elif args.remove==3:
            sdf1.remconfs([True, True])
            times.append(time.time())
            if args.verbose:
                print 'Conformation numbers from names and metafields. It took {} seconds.'.format(times[-1]-times[-2])
        
        
        if args.cut:
            for statement in args.cut.split('|'):
                sdf2 = Sdffile(statement.strip())
                sdf1.sdflistremove(sdf2,True) #sdf1.metacombi(sdf2)
                del(sdf2)
                times.append(time.time())
                if args.verbose:
                    print 'Removing molecules complete. It took {} seconds.'.format(times[-1]-times[-2])
        
        elif args.allcut:
            for statement in args.allcut.split('|'):
                sdf2 = Sdffile(statement.strip())
                sdf1.sdflistremove(sdf2,False) #sdf1.metacombi(sdf2)
                del(sdf2)
                times.append(time.time())
                if args.verbose:
                    print 'Removing molecules complete. It took {} seconds.'.format(times[-1]-times[-2])
        
        #should work, allcombine tested
        if args.combine:
            sdf2 = Sdffile(args.combine)
            sdf1.sdfmetacombi(sdf2,[True, True]) #sdf1.metacombi(sdf2)
            del(sdf2)
            times.append(time.time())
            if args.verbose:
                print 'sdf-file metadata combination complete. It took {} seconds.'.format(times[-1]-times[-2])
        
        #check
        elif args.allcombine:
            sdf2 = Sdffile(args.allcombine)
            sdf1.sdfmetacombi(sdf2,[True, False]) #metacombi(sdf2)
            del(sdf2)
            times.append(time.time())
            if args.verbose:
                print 'sdf-file metadata combination complete. It took {} seconds.'.format(times[-1]-times[-2])
        
        #check
        if args.addcsv:
            for statement in args.addcsv.split('|'):
                sdf1.addcsvmeta(statement.strip(), args.verbose)
                times.append(time.time())
                if args.verbose:
                    print 'Metadata from csv-file added. It took {} seconds.'.format(times[-1]-times[-2])
        
        if args.getmol2:
            for statement in args.getmol2.split('|'):
                try:
                    path, column, metaname = re.split('\s*,\s*', statement.strip())
                except ValueError:
                    raise ValueError('Give: path, column, metaname')
                sdf1.getMol2Data(metaname, int(column), path)
                times.append(time.time())
                if args.verbose:
                    print 'Metadata from mol2-file added. It took {} seconds.'.format(times[-1]-times[-2])
        
        
        #should work
        if args.closestatom:
            for statement in args.closestatom.split('|'):
                argus = splitter(statement)
                if len(argus)==1:
                    sdf1.closest(argus[0])
                elif len(argus)==2:
                    sdf1.closest(argus[0],name=argus[1])
                times.append(time.time())
                if args.verbose:
                    print 'Closest atoms to point calculated. It took {} seconds.'.format(times[-1]-times[-2])
        
        #should work
        if args.closeratoms:
            for statement in args.closeratoms.split('|'):
                params =  splitter(statement) #re.split('\s*,|;\s*',statement)
                sdf1.closer(params[0],params[1])
                times.append(time.time())
                if args.verbose:
                    print 'Atoms closer than {} to point calculated. It took {} seconds.'.format(params[1],times[-1]-times[-2])
        
        '''deprecation?
        if args.multiclosestatom:
            for statement in args.multiclosestatom.split('|'):
                argus = statement.split(',')
                sdf1.closestatoms(argus[0].strip(),argus[1].strip())
                times.append(time.time())
                if args.verbose:
                    print 'Closest atoms to given atom calculated. It took {} seconds.'.format(times[-1]-times[-2])
        '''
        '''
        #doesn't work
        if args.closestbybonds:
            for statement in args.closestbybonds.split('|'):
                argus = statement.split(',')
                sdf1.closestbybonds(argus[0].strip(),argus[1].strip(),argus[2].strip())
                times.append(time.time())
                if args.verbose:
                    print 'Closest atom to given atom by bond calculated. It took {} seconds.'.format(times[-1]-times[-2])
        '''
        
        #check
        if args.changemeta:
            for statement in args.changemeta.split('|'):
                argtmp = statement.split('>')
                sdf1.changemetaname(argtmp[0].strip(),argtmp[1].strip())
                times.append(time.time())
                if args.verbose:
                    print 'Changing metafield name {} to {} done. It took {} seconds.'.format(argtmp[0].strip(),argtmp[1].strip(),times[-1]-times[-2])
        
        #check
        if args.mergemeta:
            for statement in args.mergemeta.split('|'):
                sdf1.metamerger(statement.strip())
                if args.verbose:
                    print 'New metafield merged.'
            times.append(time.time())
            if args.verbose:
                print 'Merging metafields done. It took {} seconds.'.format(times[-1]-times[-2])
        
        if args.makenewmeta:
            for statement in args.makenewmeta.split('|'):
                sdf1.makenewmetastr(statement.strip())
                if args.verbose:
                    print 'New metafield made.'
            times.append(time.time())
            if args.verbose:
                print 'Making new metafields done. It took {} seconds.'.format(times[-1]-times[-2])
        
        #check
        if args.sortorder:
            sortsies = args.sortorder.split('|')
            for sorty in sortsies:
                sdf1._orderlist=sdf1.sorter(sorty)
                if args.verbose:
                    print 'Sort {} done.'.format(sorty)
                pass
            times.append(time.time())
            if args.verbose:
                print 'Sorting done. It took {} seconds.'.format(times[-1]-times[-2])
        
        if args.sortmeta:
            sortsies = args.sortmeta.split('|')
            for sorty in sortsies:
                sdf1.sortmetas(sorty)
                if args.verbose:
                    print 'Sort {} done.'.format(sorty)
                pass
            times.append(time.time())
            if args.verbose:
                print 'Meta sorting done. It took {} seconds.'.format(times[-1]-times[-2])
        
        if args.stripbutmeta:
            for statement in args.makenewmeta.split('|'):
                sdf1.stripbutmeta(statement)
                if args.verbose:
                    print 'All atoms, execpt for those in statement '+statement+' removed!'
            times.append(time.time())
            if args.verbose:
                print 'Stripping atoms done. It took {} seconds.'.format(times[-1]-times[-2])
        
        #check :)
        if args.extract:
            logics = args.extract.split('|')
            for logic in logics:
                sdf1.logicparse(logic.strip())
                if args.verbose:
                    print 'After logical chop {}, sdf-file has {} molecules and {} conformations left.'.format(logic,len(sdf1._dictomoles),len(sdf1))
            times.append(time.time())
            if args.verbose:
                print 'Logical chopping done. It took {} seconds.'.format(times[-1]-times[-2])
        
        #check
        if args.removemeta:
            sdf1.removemeta(args.removemeta, False)
            times.append(time.time())
            if args.verbose:
                print 'Given metafields removed. It took {} seconds.'.format(times[-1]-times[-2])
        #check
        elif args.pickmeta!=None:
            sdf1.removemeta(args.pickmeta, True)
            times.append(time.time())
            if args.verbose:
                print 'Non-specified metafields removed. It took {} seconds.'.format(times[-1]-times[-2])
        times.append(time.time())
        
        if args.putmol2:
            for statement in args.putmol2.split('|'):
                try:
                    input, output, column, metaname, default, precision = re.split('\s*,\s*', statement)
                except ValueError:
                    raise ValueError('Give: input, output, column, metaname, default, precision')
                sdf1.injectMol2Data(metaname, int(column), input, numify(default), int(precision), output)
                if args.verbose:
                    print 'New mol2-file with applied meta-data created.'
            times.append(time.time())
            if args.verbose:
                print 'Creating mol2-files done. It took {} seconds.'.format(times[-1]-times[-2])
        
        #check
        if args.histogram:
            showflag=True
            plots = args.histogram.split('|')
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
                sdf1.mehist(*larg,**darg)
                if path:
                    pylab.savefig(path, bbox_inches='tight')
            if showflag:
                sdf1.show()
            times.append(time.time())
            if args.verbose:
                print 'Plotting histograms done. It took {} seconds.'.format(times[-1]-times[-2])
        times.append(time.time())        
        
        wriarg=dict()
        
        #check
        if args.overwrite:
            wriarg['path']=inputfile
        elif args.output:
            wriarg['path']=args.output
        if args.getcsv:
            writetype='getcsv'
            wriarg['csv']=args.getcsv
        elif args.metalist:
            writetype='metalist'
        elif args.counts!=None:
            writetype='counts'
            wriarg['counts']=args.counts
        elif args.donotprint:
            writetype='none'
        else:
            writetype='sdf'
        
        if args.split:
            wriarg['split']=args.split
        if args.makefolder:
            wriarg['makefolder']=args.makefolder #True #makefolder = True
        sdf1.writer(writetype,**wriarg)
        times.append(time.time())
        
        if args.verbose:
                print 'File writing done. It took {} seconds.'.format(times[-1]-times[-2])
        times.append(time.time())
    
    for onefile in manyfiles:
        main(onefile)
