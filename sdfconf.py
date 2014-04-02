#!/usr/bin/env python
# -*- coding: latin-1 -*-


import os
import sys
import re
import math
import argparse
import numpy
import time
import pylab
import copy
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

ckey = "confnum"
atomcut=[0, 10, 20, 30, 31, 34, 36, 39, 42, 45, 48, 51, 54, 57, 60, 63, 66, 69]



class Sdffile:
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
        for i, molord in enumerate(self._orderlist):
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
        #Method to merge multiple metafields into one. Comes with a few operations. Frontend to mergenewmeta
        [newmeta, task] = re.split('\s*=\s*',string)
        funki = task.find('(')
        funk = task[:funki]
        metas = goodchop.split(task[funki:].strip('{[()]}'))
        
        operators = {'sum':sum, 'max':max, 'min':min, 'avg':avg, 'join':' '.join, 'prod':numpy.prod, 'sub':sub}
        if funk in operators:
            self.mergenewmeta(newmeta, metas, operators[funk.lower().strip()])
        else:
            if args.verbose:
                print 'Command {} not found'.format(funk)
    
    
    def mergenewmeta(self, newmeta, metas, oper):
        for mol in self:
            tab=[]
            for meta in metas:
                if meta in mol._meta:
                    tab.append(numify(mol._meta[meta]))
            try:
                mol.addmeta(newmeta, [str(oper(tab))])
            except TypeError:
                mol.addmeta(newmeta, [str(oper(map(str,tab)))])
                
    def metaoper(self, oper, metas):
        '''
        method to make operations for Sdfmeta objects
        return resulting metavalue of operation.
        '''
        structs = [meta._datastruct for meta in metas]
        types   = [meta._datatype   for meta in metas]
        mathopers = (sum, sub, numpy.prod, max, min, avg)
        workmetas = copy.deepcopy(metas)
        if oper in mathopers:
            if str in types:
                #cannot calculate strings, doh
                return None
            if len(set(structs))==2 and 'single' in structs:
                #singles and (list or OrDi) in structs, extend singles to lists
                ostru = iter(set(structs)-{'single'}).next() #the other structuretype

                minlen = min( map(len, {meta._data for meta in workmetas}-{1}) ) # shortest list length, ignore singles
                if ostru == list:
                    #make lists have the same length
                    for i, meta in enumerate(workmetas): 
                        if meta._datastruct == 'single':
                            workmeta[i]._data = meta._data * minlen
                        else:
                            workmeta[i]._data = meta._data[:minlen]
                    
                elif ostru == OrDi:
                    keys = set()
                    for i, meta in enumerate(workmetas):
                        if meta._datastruct == OrDi:
                            keys = keys.union(meta._data)
                    for meta in workmetas:
                        if meta._datastruct == 'single':
                            pass
        
        if oper != ' '.join:
            if len(types)>1:
                #quit if multiple types. Add support for  mixed structure, e.g., single and list
                return None
            
            if list in types:
                datas = [meta._data for meta in metas]
                minlen = min( map(len, datas) )
                datas = [data[:minlen] for data in datas]
    
    def nametoid(self, meta):
        #Adds a metafield including the molecule name
        for mol in self:
            mol.addmeta(meta, mol._name)
            
    def changemetaname(self, oldname, newname):
        #Change the name of a metafield
        for mol in self:
            mol.changemetaname(oldname, newname)
    
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
        if len(listofmeta)==1 and listofmeta[0]=='?':
            listofmeta = self.listmetas()
        csv = [separator.join(listofmeta)]
        for info in self._orderlist:
            mol = self._dictomoles[info[0]][info[1]]
            line = []#[mol._name]
            for meta in listofmeta:
                if meta in mol._meta:
                    line.append(mol._meta[meta])
                else:
                    line.append('_')
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
    

    def closest(self, point, name=None):
        coord1 = numpy.array(coorder(point)) #map(numify, goodchop.split(point.strip('{[()]}')))
        #coord1 = numpy.array(coord1)
        if not name:
            name = ''
        for mol in self:
            #is coord1 meta or not?
            alist = mol.dists(coord1)
            alist=sorted(alist, key=lambda i: i[0])
            mol.addmeta(name+'Closest_atom',str(alist[0][1]))
            mol.addmeta(name+'Closest_distance',str(alist[0][0]))


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
        #kwargs: title, Xtitle, Ytitle, bins
        newargs = dict()
        titargs = ['Xtitle','Ytitle','title']
        for key in kwargs:
            if key in titargs:
                newargs[key]=kwargs[key]
                #del(kwargs[key])
        for key in newargs:
            del(kwargs[key])
        
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


    def sorter(self,sortstring):
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
        return sorted(self._orderlist, key=lambda avain: numify(self._dictomoles[avain[0]][avain[1]]._meta[sortstring[1:]]),reverse=rever)
        
    
    def addcsvmeta(self, path):
        f = open(path)
        chop = re.compile('\s*,|\t\s*')
        csvdata = [[cell.strip('"|\n*') for cell in chop.split(line)] for line in f.readlines()] #lukee tiedoston, splittaa pilkuista ja poistaa alkioista hipsut/rivinvaihdot
        f.close()
        header = csvdata[0]
        for csvmol in csvdata[1:]:
            m = confchop.search(csvmol[0])
            if m:
                n=m.group(0)[2:-2]
                #l=len(n)+4
                #seeker = csvmol[0][:-l]
                seeker = csvmol[0][:-(len(n)+4)]
            else:
                #l=0
                seeker = csvmol[0]
            #try:
            #    moldic = self._dictomoles[seeker]
            #except KeyError:
            if not seeker in self._dictomoles:
                #print csvmol[0]
                #print csvmol[0][:-l] in self._dictomoles
                #print l
                if args.verbose:
                    print seeker+' not found!'
                continue
            if m:
                keys = [n]
            else:
                keys = self._dictomoles[seeker].keys()
            #moldic = {n:moldic[n]}
            for key in keys:
                if key in self._dictomoles[seeker]:
                    mol = self._dictomoles[seeker][key]
                else:
                    #print self._dictomoles[seeker].keys()
                    continue
                #mol = moldic[key]
                for i, newmeta in enumerate(csvmol[1:]):
                    if len(newmeta)!=0:
                        mol.addmeta(header[i+1],newmeta,True)
                        #print 'added '+header[i+1]+' to '+key
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
        
#end of Sdffile

class Sdfmole:
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
            new._counts     = copy.deepcopy(self._counts)
            new._atoms      = copy.deepcopy(self._atoms)
            new._bonds      = copy.deepcopy(self._bonds)
            new._properties = copy.deepcopy(self._properties)
        else:
            new._other      = copy.deepcopy(self._other,memo)
        return new
    
    def metahand(self, strilist):
            return strilist[0]
            
    
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
                
    def dists(self, coord1):
        self.numerize()
        alist=[]
        for i in range(len(self._atoms)):
            #if self.gettype(i+1).lower()=='h':
            #if self.gettype(i+1) in ['H','F','Br','I']:
            if self.gettype(i+1).strip() in ['H','F','Br','I']:
                continue
            coord2 = numpy.array(self.getcoord(i+1)) #numpy.array([float(c) for c in atom[0:3]])
            dist = sum((coord1-coord2)**2)**0.5
            alist.append([dist,i+1])
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
            ans[1] = self._meta[ckey].strip()
        if ans[0] and ans[1]:
            if ans[0]!=ans[1]:
                print(mes.format(ans[0],ans[1]))
        return ans
        
    def getconfn(self):
        a = self.getconf()
        return a[0] or a[1]
        #return next((item for item in a if item is not None), None)
        
    def addconf(self, conf, bolist = [True, True]):
        already = self.getconf()
        if bolist[1]:
            if already[1]:
                if conf != already[1]:
                    print(mes.format(ans[0],ans[1]))
            else:
                #self._meta[ckey] = [conf]
                self.addmeta(ckey,conf)
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
            if not overwrite:
                if key in self._meta:
                    continue
            self.addmeta(key,othersdfmol._meta[key],overwrite)
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
        
    def addmeta(self,metafield, value, overwrite = True):
        if metafield in self._meta and not overwrite:
            return
        if type(value) == list:
            #value = [value]
            value = ' '.join(value)
        if not metafield in self._meta:
            self._metakeys.append(metafield)
        self._meta[metafield] = value
        #print metafield+' '+self._meta[metafield]
        
    def changemetaname(self, oldname, newname):
        if oldname in self._meta:
            i=self._metakeys.index(oldname)
            self._metakeys[i]=newname
            self._meta[newname]=self._meta[oldname]
            del(self._meta[oldname])
            
    def atomlistdistances(self,metaatom,metalist):
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
    
    def pointlistdists(self, point,metalist):
        arpoint=numpy.array(point)
        self.numerize()
        try:
            manyatom=self._meta[metalist]
        except KeyError:
            return None
        atoms = [numify(one) for one in re.split('\s+',manyatom.strip())]
        dists = []
        for oneatom in atoms:
            dists.append(self.disttopoint(oneatom,arpoint))
        return (atoms,dists)
    
    def disttopoint(self, atom1, point):
        return sum((self.getatomloc(atom1)-point)**2)**0.5
        
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
    
#end of Sdfmole

class Sdfmeta:
    #name
    #datatype: string list dict ?
    #the structure
    #                 |------ singles ? -----|
    #self._datypes = ['string', 'int', 'float' 'vlist', 'hlist', 'contlist', 'dict','unknown']
    
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
        if isinstance( ind, slice ) or isinstance( ind, int ) :
            return mylist[ind]
        else:
            raise TypeError, "Invalid argument type."

        
    def initialize(self, listofstrings):
        '''parse the metadata'''
        #get name of metafield
        if listofstrings[0][0] != '>' :
            #Raise error
            print listofstrings[0]
            #print 'No meta'
            return
        else:
            #print listofstrings[0]
            self._name = metaname.match(listofstrings[0]).groups()
            #print self._name
        #Remove the last empty line and linechanges
        if not len(listofstrings[-1].strip())==0:
            #Raise error?
            return
        else: #             rstrip?
            mylines = [line.strip('\n ') for line in listofstrings[1:-1] ]
        
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
            
    def construct(self, data, delims=None):
        new = Sdfmeta()
        if type(data) in (Ordi, dict, list):
            if type(data)==dict:
                data = OrDi(data)
                new._datastruct = OrDi
            elif type(data)==list and len(data)==1:
                new._datastruct = 'single'
            else:
                new._datastruct = type(data)
        elif type(data) in (str, int, float):
            new._datastruct = 'single'
            data = [data]
        else:
            raise TypeError('Datastructure type not list, dict, OrderedDict, str, int or float.')
        
        if type(new._data) == list:
            data = map(numify, data)
            types = set(map(type, data))
        else:
            for key in data:
                newkey = numify(key)
                data[newkey] = numify(data[key])
                if key != newkey:
                    del(data[key])
            types = {type(data[key]) for key in data}
            
        if not all({typ in (int, float, str) for typ in types}):
            raise TypeError('Data type not str, int or float.')
        if len(types)>1 and str in types:
            raise TypeError('Mixed datatypes. Strings and numbers.')
        elif len(types)==1:
            newtype = iter(types).next()
        else:
            newtype = float
            if type(new._data) == list:
                
            
            
    def __copy__(self):
        new = Sdfmeta()
        new._name = copy.copy( self._name )
        new._datatype = copy.copy( self._datatype )
        new._datastruct = copy.copy( self._datastruct )
        new._data = copy.copy( self._data )
        new._delims = copy.copy( self._delims )
        return new
        
    def __deepcopy__(self):
        new = Sdfmeta()
        new._name = copy.deepcopy( self._name )
        new._datatype = copy.deepcopy( self._datatype )
        new._datastruct = copy.deepcopy( self._datastruct )
        new._data = copy.deepcopy( self._data )
        new._delims = copy.deepcopy( self._delims )
        return new
    
    def getname(self):
        return self._name[1]
        
    def dtype(self):
        return self._datatype
        
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
            lde = len(self._delims[i])
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
        
    def cleandelim(self):
        '''Clean excess whitespaces from delimiters'''
        for i in range(len(delim)):
            #new = 
            #if len(new)==0:
            #    new = ' '
            self._delim[i] = self._delim[i].strip() + ' '
            
#End of Sdfmeta

def coorder(point):
    return map(numify, goodchop.split(point.strip('{[()]}')))


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
    try:
        return int(stri)
    except ValueError:
        try:
            return float(stri)
        except ValueError:
            return stri
    except TypeError:
        if type(stri)==list:
            return map(numify,stri)
        else:
            raise TypeError('Wrong type')
            
def avg(num):
    return float(sum(num))/len(num)
    
def sub(num):
	if len(num)<2:
		return num
	else:
		return num[0]-sum(num[1:])

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
    for i, item in enumerate(original):
        sta=list(stapa.finditer(item))
        end=list(endpa.finditer(item))
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
    if string.find('[')==0 and (string.rfind(']')-len(string))==-1:
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
        
def whattype(onestring):
    '''
    Tries to find out what type your string is
    Return (parsed,parsedtype,delim)
    '''
    #onestring = onestring.strip()
    #Number?
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
        tab.append(oper(caltab))
        i+=1
    return tab
    
def OrDioper(oper, metas):
    keys = set.intersection(*[set(meta._data) for meta in metas if meta._datastruct != 'single'])
    dic = dict()
    for key in keys:
        caltab = []
        for meta in metas:
            if meta._datastruct == 'single':
                caltab.append(meta._data[0])
            else:
                caltab.append(meta._data[key])
        dic[key] = oper(caltab)
    
    for meta in metas:
        refmeta = meta
        if refmeta._datastruct == OrDi:
            break
    if refmeta._datastruct != OrDi:
        return None
    outordi = OrDi()
    for key in refmeta._datastruct:
        if key in dic:
            outordi[key]=dic[key]
    return outordi
    
def singleoper(oper, metas):
    if len({meta._datastruct for meta in metas} - {single})>0:
        return None
    else:
        return [oper([meta._data[0] for meta in metas])]

if __name__ == "__main__":
    
    arger = argparse.ArgumentParser(description= 'Some bad-ass manipulation of SDF-files.' )
    
    arger.add_argument("input", metavar = 'path', nargs='+', type = str,               help="Specify the input file")
    choicewrite = arger.add_mutually_exclusive_group()
    choicewrite.add_argument("-out", "--output", type = str, default=None,  help = "Specify output file.")
    choicewrite.add_argument("-o", "--overwrite", action='store_true',      help = "overwrite. you don't need to specify output file")
    arger.add_argument("-cf", "--tofield", action = "store_true",           help = "add conformation number to metadata from name. if number in name doesn't exist, make a new one.")
    arger.add_argument("-cn", "--toname",  action = "store_true",           help = "add conformation number to name from metadata. if number in metadata doesn't exist, make a new one.")
    arger.add_argument("-mn", "--metatoname",  type = str,                  help = "Change the name of molecule to the data in given metafield")
    arger.add_argument("-ni", "--nametoid",  type = str,                    help = "Copy the name of molecule into given metafield")
    arger.add_argument("-rc", "--remove",  type = int,                      help = "remove conformation number from 1=metadata, 2=name, 3=both. if number doesn't exist, do nothing")
    choicecombi = arger.add_mutually_exclusive_group()
    choicecombi.add_argument("-co", "--combine", type = str,                help = "Combine metadata from specified file to the data of original file. Confromationr_epik_Ionization_Penaltys must match.")
    choicecombi.add_argument("-aco", "--allcombine", type = str,            help = "Combine metadata from specified file to the data of original file. Names must match.")
    choicecut = arger.add_mutually_exclusive_group()
    choicecut.add_argument("-cu", "--cut", type = str,                      help = "Remove molecules in specified file from original file. Confromations must match.")
    choicecut.add_argument("-acu", "--allcut", type = str,                  help = "Remove molecules in specified file from original file. Names must match. Not tested")
    arger.add_argument("-csv", "--addcsv", type = str,                      help = "Add metadata from csv-file. File must have a 1-line header, it gives names to metafields. Names of molecules must be leftmost. If name includes confnumber, meta is only added molecules with same confnumber.")
    arger.add_argument("-ex", "--extract", type = str,                      help = "Give logical statement [metafield]logicaloperator[number]. For example [molweight<550]. [Shaepresult>90%%T] gives those who have 90%% value of top-molecule, works always per molecule name. A=average, T=top, B=bottom. M gives %% of molecules, while N gives number of molecules. ie. [dataresult>5%%N] gives 5 molecules per name with highest dataresult value. ~ makes comparison between metafield values, = makes a search for fixed value on datafield. Takes multiple statements separated with | ")
    choiceremo = arger.add_mutually_exclusive_group()
    choiceremo.add_argument("-rm", "--removemeta", type = str,              help = "Remove metadata from molecules. Takes multiple values, separaterd by comma(,) or semicolon(;)")
    choiceremo.add_argument("-pm", "--pickmeta", type = str,                help = "Remove all nonspecified metadata from molecules. Takes multiple values, separaterd by comma(,) or semicolon(;)")
    arger.add_argument("-s", "--split", type = int,                         help = "Split the file into even pieces. Positibe means number of file while negative means number of molecules per file. 0 doesn't apply.")
    arger.add_argument("-mf", "--makefolder", action = "store_true",        help = "Put outputfile(s) into folder(s).")
    outputtype = arger.add_mutually_exclusive_group()
    outputtype.add_argument("-gc", "--getcsv", type = str ,                 help = "Writes a .csv-file istead of .sdf-file. Specify which fields you'll need")
    outputtype.add_argument("-ml", "--metalist", action = "store_true",     help = "Writes a list of metafields.")
    outputtype.add_argument("-nm", "--counts", type = int,                  help = "Number of different molecules and different conformations. 0=just sums, 1=by molecule name, 2=both")
    outputtype.add_argument("-dnp", "--donotprint", action = "store_true",  help = "No output")

    arger.add_argument("-ca", "--closestatom", type = str,                  help = "Calculates the closest atom (atom number) to given point. Adds 'Closest_atom' and 'Closest_distance' metafields. Needs coordinates, separated by , or ;")
    arger.add_argument("-mc", "--multiclosestatom", type = str,             help = "Calculates the distances from atom number in a given field, multiple atoms in second field. Creates a field <[secondfieldname]_closest_distance> ['Closest_atom,PRIMARY_SOM|Closest_atom,SECONDARY_SOM]")
    arger.add_argument("-cb", "--closestbybonds", type = str,               help = "Calculates the closest atom of interest (atom number) from given atom (fieldname) to given atoms (fieldname). Adds '<input3>_atom' and '<input3>_distance' metafields. Needs three fieldnames, separated by ','. Fields are 'from', 'to' and 'newfield'. You can give multiple parameters, separate by '|'")
    arger.add_argument("-cla", "--closeratoms", type = str,                 help = "Calculates number of atoms closer to the given point, than the ones given adds metafields 'Closest_atom_from_{meta}' and 'Closer_atoms_than_{meta}'. Needs point and metafield name separated by ',', point first. Takes multiple parameters separated by '|'")
    arger.add_argument("-v", "--verbose", action = "store_true" ,           help = "More info on your run.")
    arger.add_argument("-mm", "--mergemeta", type = str,                    help = "Makes a new metafield based on old ones. newmeta=sum(meta1,meta2). operator are sum, max, min, avg, prod, join. Takes multiple arguments separated by |")
    arger.add_argument("-cm", "--changemeta", type = str,                   help = "Changes names of metafields. [olname1>newname1|oldname2>newname2]. ")
    arger.add_argument("-so", "--sortorder", type = str,                    help = "Sorts molecules of a file in order of metafield. <MolecularWeight|>Id Sorts molecules first by highest weight first, then by smallest name first")
    arger.add_argument("-hg", "--histogram", type = str,                    help = "Plots a 1D or 2D histogram, multiple plots with '|'. 'Xname,Yname,Xtitle=x-akseli,Ytitle=y-akseli,bins=[30,30]'")
    
    
    args = arger.parse_args()
    manyfiles = args.input #glob.glob(args.input)
    
    def main(inputfile):
        times = [time.time()]
        sdf1 = Sdffile(inputfile)
        times.append(time.time())
        
        if args.verbose:
            print 'Reading file done. It took {} seconds.'.format(times[-1]-times[-2])
        
        if args.tofield:
            sdf1.addconfs([False,True])
            times.append(time.time())
            if args.verbose:
                print 'Conformation numbers added to metadata. It took {} seconds.'.format(times[-1]-times[-2])
        
        if args.toname:
            sdf1.addconfs([True,False])
            times.append(time.time())
            if args.verbose:
                print 'Conformation numbers added to names. It took {} seconds.'.format(times[-1]-times[-2])
        
        if args.nametoid:
            sdf1.nametoid(args.nametoid)
            times.append(time.time())
            if args.verbose:
                print 'Name written to ID-field. It took {} seconds.'.format(times[-1]-times[-2])
        
        if args.remove==2:
            sdf1.remconfs([True, False])
            times.append(time.time())
            if args.verbose:
                print 'Conformation numbers removed from names. It took {} seconds.'.format(times[-1]-times[-2])
        
        elif args.remove==1:
            sdf1.remconfs([False, True])
            times.append(time.time())
            if args.verbose:
                print 'Conformation numbers removed from metafields. It took {} seconds.'.format(times[-1]-times[-2])
        
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
        
        if args.combine:
            sdf2 = Sdffile(args.combine)
            sdf1.sdfmetacombi(sdf2,[True, True]) #sdf1.metacombi(sdf2)
            del(sdf2)
            times.append(time.time())
            if args.verbose:
                print 'sdf-file metadata combination complete. It took {} seconds.'.format(times[-1]-times[-2])
        
        elif args.allcombine:
            sdf2 = Sdffile(args.allcombine)
            sdf1.sdfmetacombi(sdf2,[True, False]) #metacombi(sdf2)
            del(sdf2)
            times.append(time.time())
            if args.verbose:
                print 'sdf-file metadata combination complete. It took {} seconds.'.format(times[-1]-times[-2])
        
        if args.addcsv:
            for statement in args.addcsv.split('|'):
                sdf1.addcsvmeta(statement.strip())
                times.append(time.time())
                if args.verbose:
                    print 'Metadata from csv-file added. It took {} seconds.'.format(times[-1]-times[-2])
        
        if args.closestatom:
            for statement in args.closestatom.split('|'):
                argus = splitter(statement)
                sdf1.closest(*argus)
                times.append(time.time())
                if args.verbose:
                    print 'Closest atoms to point calculated. It took {} seconds.'.format(times[-1]-times[-2])

        if args.closeratoms:
            for statement in args.closeratoms.split('|'):
                params =  splitter(statement)#re.split('\s*,|;\s*',statement)
                sdf1.closer(params[0],params[1])
                times.append(time.time())
                if args.verbose:
                    print 'Atoms closer than {} to point calculated. It took {} seconds.'.format(params[1],times[-1]-times[-2])

        if args.multiclosestatom:
            for statement in args.multiclosestatom.split('|'):
                argus = statement.split(',')
                sdf1.closestatoms(argus[0].strip(),argus[1].strip())
                times.append(time.time())
                if args.verbose:
                    print 'Closest atoms to given atom calculated. It took {} seconds.'.format(times[-1]-times[-2])
        
        if args.closestbybonds:
            for statement in args.closestbybonds.split('|'):
                argus = statement.split(',')
                sdf1.closestbybonds(argus[0].strip(),argus[1].strip(),argus[2].strip())
                times.append(time.time())
                if args.verbose:
                    print 'Closest atom to given atom by bond calculated. It took {} seconds.'.format(times[-1]-times[-2])

        if args.changemeta:
            for statement in args.changemeta.split('|'):
                argtmp = statement.split('>')
                sdf1.changemetaname(argtmp[0].strip(),argtmp[1].strip())
                times.append(time.time())
                if args.verbose:
                    print 'Changing metafield name {} to {} done. It took {} seconds.'.format(argtmp[0].strip(),argtmp[1].strip(),times[-1]-times[-2])
        
        if args.mergemeta:
            for statement in args.mergemeta.split('|'):
                sdf1.metamerger(statement.strip())
                if args.verbose:
                    print 'New metafield merged.'
            times.append(time.time())
            if args.verbose:
                print 'Merging metafields done. It took {} seconds.'.format(times[-1]-times[-2])
        
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
        
        if args.extract:
            logics = args.extract.split('|')
            for logic in logics:
                sdf1.logicchop(logic.strip())
                if args.verbose:
                    print 'After logical chop {}, sdf-file has {} molecules left.'.format(logic,len(sdf1))
            times.append(time.time())
            if args.verbose:
                print 'Logical chopping done. It took {} seconds.'.format(times[-1]-times[-2])
        
        if args.removemeta:
            sdf1.removemeta(args.removemeta, False)
            times.append(time.time())
            if args.verbose:
                print 'Given metafields removed. It took {} seconds.'.format(times[-1]-times[-2])
        elif args.pickmeta!=None:
            sdf1.removemeta(args.pickmeta, True)
            times.append(time.time())
            if args.verbose:
                print 'Non-specified metafields removed. It took {} seconds.'.format(times[-1]-times[-2])
        times.append(time.time())
        
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
