#!/usr/bin/env python
# -*- coding: latin-1 -*-

myversion='0.801'

#import os
import sys
#import re
#import math
import argparse
#import numpy
import time
#import copy
#import operator
from collections import OrderedDict as OrDi
#import warnings
#import bisect as bi

try:
    from future.utils import lmap
except ImportError:
    lmap = map

if sys.version_info[0]==2 and sys.version_info[1]>=7:
    pass
else:
    raise SystemError('Python version must be 2.7. or later, but not 3.x.')

from sdfconf.sdf import Sdffile
#import functions
from sdfconf import functions

class Runner(object):

    order = OrDi( (
             ('verbose','v'), 
             ('input','in'), 
             ('tofield','cf'), 
             ('toname','cn'), 
             ('nametometa','ntm'), 
             ('removeconfname','rcn'), 
             ('removeconfmeta','rcm'), 
             ('ignores', 'ig'),
             ('cut','cu'), 
             ('allcut','acu'), 
             ('combine','co'), 
             ('allcombine','aco'), 
             ('addcsv','csv'), 
             ('getmol2','gm2'), 
             ('addescape','aesc'), 
             ('addinside','ains'), 
             ('config','con'), 
             ('closestatoms','ca'), 
             ('closeratoms','cla'), 
             ('changemeta','cm'), 
             ('makenewmeta','mnm'), 
             ('sortmeta','sm'), 
             ('stripbutmeta','sbm'), 
             ('proportion','pro'), 
             ('sortorder','so'), 
             ('proportion','pro'), 
             ('extract','ex'), 
             ('metatoname','mtn'), 
             ('removemeta','rm'), 
             ('pickmeta','pm'), 
             ('putmol2','pm2'), 
             ('histogram','hg'), 
             ('getcsv','gc'), 
             ('getatomcsv','gac'), 
             ('metalist','ml'), 
             ('counts','nm'), 
             ('donotprint','dnp'), 
             ('split','s'), 
             ('makefolder','mf'), 
             ('output','out'), 
             ('overwrite','o'), 
             ) )
    
    
    simplelist =    ('tofield', 'toname', 'removeconfname', 'removeconfmeta', 
                     'ignores', 'nametometa', 'metatoname', 'removemeta', 
                     'pickmeta', 'input', 'proportion' , 'histogram', 'verbose', 
                     )
    
    simpleloops =   ('getmol2', 'closestatoms', 'closeratoms', 'changemeta', #'mergemeta',
                     'sortorder', 'sortmeta', 
                     'stripbutmeta', 'extract','makenewmeta','config',
                     'cut', 'allcut', 'combine', 'allcombine', 'addcsv', 
                     'putmol2', 'addescape', 'addinside', 
                     )
    
    listbatch = ()
    
    writers =       ('overwrite', 'output', 'stdout') #default stdout
    writetypes =    {'getcsv':True,'getatomcsv':True,'metalist':True,'counts':True,'donotprint':True,'sdf':True,'split':False,'makefolder':False} #default 'sdf'
    #graphers =      ('histogram')
    
    def __init__(self, options=dict()):
        self.options = options if options else dict()
        
        
        for option, value in (('verbose', False),('proportion', False),('ignores',[])):
            if not option in self.options:
                self.options[option] = value
        
        self.times = [time.time()]
        self.inpath = None
        self.wriarg = dict()
        self.writetype = 'sdf'
        self.writer = None
        
        self.setVerbose()
        self.setIgnores()
        
        self.sdf = Sdffile()
        
    def setPropor(self,parameter=None):
        if not parameter or parameter is True:
            self.propor = False
        else:
            self.propor = parameter.strip()
    
    def setVerbose(self, verbose=False):
        if isinstance(verbose, (list, tuple)):
            self.verbose = verbose[0]
        else:
            self.verbose = verbose
    
    def setIgnores(self,ignore=['H']):
        self.ignore = ignore
        try:
            self.sdf.setIgnores(self.ignore)
        except AttributeError:
            pass
    
    
    def runOptions(self):
        #after config, do the rest anyway.
        #print self.options
        for option in Runner.order:
            if option in self.options:
                self.times.append(time.time())
                self.funcselector(option, self.options[option])
        
        if self.verbose:
            print('Run complete, it took {:.3f} seconds.'.format(self.times[-1]-self.times[0]))
        
    def runConfig(self,confpath):
        
        def parsecon(oneline):
            
            params = [para.strip() for para in oneline.partition('#')[0].partition('::')]
            if params[0] == '':
                return None
            if params[0] not in Runner.order:
                for longname in Runner.order:
                    if params[0] == Runner.order[longname]:
                        params[0] = longname
                        break
            if params[1] == '::' and params[2] != '':
                return (params[0], [functions.numify(cell.strip()) for cell in params[2].split(';;')])
            else:
                return (params[0], (True,))
        with open(confpath, 'r') as confile:
            configures = (parsecon(confline)  for confline in confile)
            for option in configures:
                if option is None:
                    continue
                elif option[0] in Runner.order: #self.options:
                    self.times.append(time.time())
                    self.funcselector(option[0], option[1], True)
        
    
    def taskLib(self,task,option,param=None,**kwargs):
        steps = kwargs.get('steps',1)
        timedif = lambda : '{:.3f}'.format(self.times[-1]-self.times[-(1+steps)])
        
        NoLamb = lambda : None
        def extMes(n):
            mypropo = self.sdf.propomes(self.propor)
            mypropo = mypropo.replace('{','{{')
            mypropo = mypropo.replace('}','}}')
            messes = (
                      lambda : ('\n'.join(('After logical chop {}, sdf-file has {} molecules and {} conformations left.',mypropo)).strip() ,(param,len(self.sdf._dictomoles),len(self.sdf))), 
                      lambda : ('\n'.join(('Initially sdf-file has {} molecules and {} conformations.',mypropo)).strip(),(len(self.sdf._dictomoles),len(self.sdf))), 
                      lambda : ('Chopping complete, it took {} seconds.',(timedif(),)) 
                      )
            return messes[n]
            
        tasks =   {
                    'tofield'        :lambda i : (self.sdf.addconfs, lambda : ((False,True),),      NoLamb,NoLamb,lambda : ('Conformation numbers added to metadata. It took {} seconds.', (timedif(),)))[i], 
                    'toname'         :lambda i : (self.sdf.addconfs, lambda : ((True,False),),      NoLamb,NoLamb,lambda : ('Conformation numbers added to names. It took {} seconds.', (timedif(),)))[i], 
                    'removeconfname' :lambda i : (self.sdf.remconfs, lambda : ((True,False),),      NoLamb,NoLamb,lambda : ('Conformation numbers removed from names. It took {} seconds.',(timedif(),)))[i], 
                    'removeconfmeta' :lambda i : (self.sdf.remconfs, lambda : ((True,False),),      NoLamb,NoLamb,lambda : ('Conformation numbers removed from metafield \'confnum\'. It took {} seconds.',(timedif(),)))[i], 
                    'nametometa'     :lambda i : (self.sdf.nametometa, lambda : (param,),           NoLamb,NoLamb,lambda : ('Name written to metafieldfield {}. It took {} seconds.',(param,timedif())))[i], 
                    'metatoname'     :lambda i : (self.sdf.metatoname, lambda : (param,),           NoLamb,NoLamb,lambda : ('Metafield written to name {}. It took {} seconds.',(param,timedif())))[i],
                    'removemeta'     :lambda i : (self.sdf.removemeta, lambda : (param, False),      lambda : ('Metafieldfield(s) in ({}) Removed. It took {} seconds.',(param,timedif())),NoLamb,NoLamb)[i], 
                    'pickmeta'       :lambda i : (self.sdf.removemeta, lambda : (param, True),       lambda : ('All metafields except those in ({}) Removed. It took {} seconds.',(param,timedif())),NoLamb,NoLamb)[i], 
                    'getmol2'        :lambda i : (self.sdf.getMol2DataStr, lambda : (param,),       NoLamb,NoLamb,NoLamb)[i], 
                    'closestatoms'   :lambda i : (self.sdf.closestStr, lambda : (param,),                  NoLamb,NoLamb,NoLamb)[i], 
                    'closeratoms'    :lambda i : (lambda x: self.sdf.closer(*functions.splitter(x)), lambda : (param,),     NoLamb,NoLamb,NoLamb)[i], 
                    'changemeta'     :lambda i : (lambda x: self.sdf.changemetaname([item.strip() for item in x.split('>')]), lambda : (param,), NoLamb, NoLamb, NoLamb)[i], 
                    'sortorder'      :lambda i : (self.sdf.sortme, lambda : (param,),               lambda : ('Sort {} done.',(param,)),NoLamb,NoLamb)[i], 
                    'stripbutmeta'   :lambda i : (self.sdf.stripbutmeta, lambda : (param,),         lambda : ('All atoms, execpt for those in statement {} removed!',(param,)),NoLamb,NoLamb)[i], 
                    'extract'        :lambda i : (self.sdf.mollogicparse, lambda : (param,),        extMes(0),extMes(1),extMes(2))[i],
                    'makenewmeta'    :lambda i : (self.sdf.makenewmetastr, lambda : (param,),       lambda : ('New metafield {} made.',(param,)), NoLamb, lambda : ('Making new metafields done. It took {} seconds.',(timedif(),)))[i], 
                    'addcsv'         :lambda i : (self.sdf.addcsvmeta, lambda : (param, self.verbose,), lambda : ('Metadata from csv-file {} added. It took {} seconds.',(param,timedif(),)),NoLamb,NoLamb)[i], 
                    'input'          :lambda i : (self.sdf.xreadself, lambda : (param,),            NoLamb, lambda : ('Starting to read file {}',(param,)), lambda : ('Reading file done. It took {} seconds.', (timedif(),)))[i],
                    'verbose'        :lambda i : (self.setVerbose, lambda : (param,) ,NoLamb, NoLamb, lambda : ('Verbose enabled.' if param else 'Verbose disabled.',(param,)))[i], 
                    'ignores'        :lambda i : (self.setIgnores, lambda : (param,) ,NoLamb, NoLamb, lambda : ('Ignores set to {}.',(param,)))[i], 
                    'config'         :lambda i : (self.runConfig, lambda : (param,),                lambda : ('Config-file {} done.',(param,)), lambda : ('Run config-files.',()), lambda : ('Running config-files done',()))[i], 
                    'histogram'      :lambda i : (self.sdf.histogramFromList, lambda : (param,),    NoLamb, lambda : ('Start plotting',()), lambda : ('Plotting done',()))[i], 
                    'cut'            :lambda i : (lambda x :self.sdf.listremove(Sdffile(x.strip()),True), lambda : (param,), lambda : ('Removing all molecules (matching name) present in {} complete. It took {} seconds.', (param, timedif())), NoLamb, NoLamb)[i], 
                    'allcut'         :lambda i : (lambda x :self.sdf.listremove(Sdffile(x.strip()),False), lambda : (param,), lambda : ('Removing all molecules (matching confnum) present in {} complete. It took {} seconds.', (param, timedif())), NoLamb, NoLamb)[i], 
                    'combine'        :lambda i : (lambda x :self.sdf.sdfmetacombi(Sdffile(x.partition(';')[0].strip()),(True,True),x.partition(';')[2].strip().lower() in ('o','over','overwrite')), lambda : (param,), lambda : ('Combining metadata from {} (matching confnum) complete. It took {} seconds.', (param, timedif())), NoLamb, NoLamb)[i], 
                    'allcombine'     :lambda i : (lambda x :self.sdf.sdfmetacombi(Sdffile(x.partition(';')[0].strip()),(True,False),x.partition(';')[2].strip().lower() in ('o','over','overwrite')), lambda : (param,), lambda : ('Combining metadata from {} (matching name) complete. It took {} seconds.', (param, timedif())), NoLamb, NoLamb)[i], 
                    'proportion'     :lambda i : (self.setPropor, lambda : (param,) ,NoLamb, NoLamb, (lambda : ('Propor set to {}.',(param,))) if param else NoLamb )[i], 
                    'putmol2'        :lambda i : (lambda x :self.sdf.injectMol2DataStr(x.strip()), lambda : (param,), NoLamb, NoLamb, NoLamb)[i],
                    'addescape'      :lambda i : (self.sdf.escapeStr, lambda : (param, False), lambda : ('Escape number from {} added. It took {} seconds.',(param, timedif())), lambda : ('Calculating escape numbers...', ()), NoLamb)[i],
                    'addinside'      :lambda i : (self.sdf.escapeStr, lambda : (param, True), lambda : ('Inside number from {} added. It took {} seconds.',(param, timedif())), lambda : ('Calculating inside numbers...', ()), NoLamb)[i],
                   }
        selector = {'func':0,'loop':2,'initial':3,'final':4}
        
        try:
            if task == 'func':
                mytask=tasks.get(option)
                fu = mytask(0)
                para = mytask(1)()
            else:
                fu=tasks.get(option)(selector.get(task))
                para=()
        except KeyError:
            raise KeyError('Wrong task={} or option={}'.format(task,option))
        retu = fu(*para)
        if task == 'func':
            self.times.append(time.time())
        return retu
        
    def funcselector(self,option,params, fromconfig=False):
        
        def messenger(task, mypara=params,**kwargs):
            if self.verbose:
                mes = self.taskLib(task, option, mypara, **kwargs)
                if mes:
                    print(mes[0].format(*mes[1]))
        
        if option == 'input':
            messenger('initial')
            self.inpath = params
            self.taskLib('func', option, params)
            messenger('final')
        
        elif option in Runner.simpleloops or option in Runner.simplelist:
            if not fromconfig and option in Runner.simplelist:
                
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
            self.times.append(time.time())
            if self.verbose:
                print('File writing done, it took {} seconds.'.format(self.times[-1]-self.times[-2]))
    
#End of Runner

if __name__ == "__main__":
    
    #main should only collect arguments...
    
    arger = argparse.ArgumentParser(description='Some bad-ass manipulation of SD-files. Also data retrieval/injection for .mol2-files. Notice that documentation is not completely up to date. Version '+myversion)
    
    arger.add_argument("input", metavar = 'input.sdf', nargs='*', type = str, default = None,  help="Specify the input  sd-file")
    
    choicewrite = arger.add_mutually_exclusive_group()
    choicewrite.add_argument("-out", "--output", metavar='output.file', type = str, default=None,  help = "Specify output file. It may be sdf or csv, depending on other arguments.")
    choicewrite.add_argument("-o", "--overwrite", action='store_true',      help = "Overwrite to input file. You don't need to specify output file")
    
    arger.add_argument("-v", "--verbose", action = "store_true" ,      help = "More info on your run.")
    
    arger.add_argument("-con","--config", metavar = 'config.txt', nargs='+', type = str, help="Specify the config file. Config file includes lines of argument name, followed by '::' and argument value. Separate multiple values with ';;'.")
    arger.add_argument("-cf", "--tofield", action = "store_true",           help = "Add conformation number to metafielf 'confnum'. If number in name doesn't exist, makes a new one.")
    
    arger.add_argument("-cn", "--toname",  action = "store_true",           help = "add conformation number to name from metafield 'confnum'. If number in metafield doesn't exist, make a new one.")
    arger.add_argument("-mtn", "--metatoname", metavar='meta' , type = str,                 help = "Change the name of molecule to the data in given metafield.")
    arger.add_argument("-ntm", "--nametometa",  metavar='meta', type = str,                 help = "Copy the name of molecule into given metafield.")
    arger.add_argument("-rcn", "--removeconfname",  action="store_true",    help = "remove conformation number from name.")
    arger.add_argument("-rcm", "--removeconfmeta",  action="store_true",    help = "remove conformation number from metafield 'confnum'.")
    
    arger.add_argument("-aesc", "--addescape",  type=str, nargs='+', metavar="File,mol_num,range,name[,max]", help = "Add metafield name_escapenum which is a list of atoms not in range of atoms in other molecule. mol_num specifies which molecule in a file you want to compare to, starts from 0. With option max quits if max atoms are found, 0 means all.")
    arger.add_argument("-ains", "--addinside",  type=str, nargs='+', metavar="File,mol_num,range,name[,max]", help = "Add metafield name_insidenum which is a list of atoms in range of atoms in other molecule. mol_num specifies which molecule in a file you want to compare to, starts from 0. With option max quits if max atoms are found, 0 means all.")
    arger.add_argument("-co", "--combine", metavar='addition.sdf', type = str, nargs='+',     help = "Combine metadata from specified file to the data of original file. Confromation numbers must match.")
    arger.add_argument("-aco", "--allcombine", metavar='addition.sdf', type = str, nargs='+', help = "Combine metadata from specified file to the data of original file. Names must match.")
    
    arger.add_argument("-cu", "--cut", metavar='unwanted.sdf', type = str, nargs='+',           help = "Remove molecules in specified file from original file. Confromations must match.")
    arger.add_argument("-acu", "--allcut", metavar='unwanted.sdf', type = str, nargs='+',       help = "Remove molecules in specified file from original file. Names must match. Not tested.")
    
    arger.add_argument("-csv", "--addcsv", metavar='data.csv', type = str, nargs='+',           help = "Add metadata from csv-file. File must have a 1-line header, it gives names to metafields. Names of molecules must be leftmost. If name includes confnumber, meta is only added molecules with same confnumber.")
    arger.add_argument("-ex", "--extract", metavar='statement', type = str, nargs='+',           help = "Pick or remove molecules from file by metafield info. Either with logical comparison or fraction of molecules with same name. Closest_atoms{:5}==soms, 2.5>Closest_atoms(soms)[], Closest_atoms[:3]<5.5, ID='benzene'.")
    arger.add_argument("-pro", "--proportion", metavar='statement', type = str,                  help = "Takes one exctract-like metastatement and prints proportion of molecules and conformations fulfilling it after every chop, if you are on verbose mode.")
    
    choiceremo = arger.add_mutually_exclusive_group()
    choiceremo.add_argument("-rm", "--removemeta", metavar='unwanted1,unwanted2,...', type = str,              help = "Remove metadata from molecules. Takes multiple values, separaterd by comma(,) or semicolon(;). If first is '?', means 'all but'.")
    choiceremo.add_argument("-pm", "--pickmeta", metavar='wanted1,wanted2,...', type = str,                help = "Remove all nonspecified metadata from molecules. Takes multiple values, separaterd by comma(,) or semicolon(;). If first is '?', means 'all but'.")
    
    arger.add_argument("-s", "--split", type = int,                         help = "Split the file into even pieces. Positive means number of files while negative means number of molecules per file. 0 doesn't apply.")
    arger.add_argument("-mf", "--makefolder", action = "store_true",        help = "Put outputfile(s) into folder(s).")
    
    outputtype = arger.add_mutually_exclusive_group()
    outputtype.add_argument("-gc", "--getcsv", type = str ,                 help = "Writes a .csv-file istead of .sdf-file. Specify which fields you'll need, separated by ','.")
    outputtype.add_argument("-gac", "--getatomcsv", type = str ,            help = "Writes a .csv-file istead of .sdf-file. Specify which fields you'll need, separated by ','. Writes dictionaries to separate lines.")
    outputtype.add_argument("-ml", "--metalist", action = "store_true",     help = "Writes a list of metafields.")
    outputtype.add_argument("-nm", "--counts", nargs='?', type = int, const=0, choices=(0,1,2),  help = "Number of different molecules and different conformations. 0=just sums, 1=by molecule name, 2=both.")
    outputtype.add_argument("-dnp", "--donotprint", action = "store_true",  help = "No output")
    
    arger.add_argument("-ca", "--closestatoms", type = str, nargs='+',  metavar='(x, y, z) [,name] [,interests=value]', help = "Calculates the closest atoms (distances by atom number) to given point. Creates a metafield with given name, if no name is given 'Closest_atoms' is created. (xx, yy, zz) may be replaced by metastatement describing single atom number.")
    arger.add_argument("-cla", "--closeratoms", type = str, nargs='+', metavar= "(x, y, z),meta", help = "Calculates number of atoms closer to the given point, than the ones given in meta. Adds metafields 'Closest_atom_from_{meta}' and 'Closer_atoms_than_{meta}'.")
    
    arger.add_argument("-mnm", "--makenewmeta", type = str, nargs='+', metavar='newmeta=metastatement',     help = "Makes a new metafield based on metastatement.")
    arger.add_argument("-cm", "--changemeta", type = str, nargs='+',   metavar='olname1>newname1' , help = "Changes names of metafields. [olname1>newname1|oldname2>newname2].")
    arger.add_argument("-so", "--sortorder", type = str, nargs='+',    metavar='meta', help = "Sorts molecules of a file in order of metafield. <MolecularWeight|>Id Sorts molecules first ascending by weight, then descenting by name.")
    arger.add_argument("-hg", "--histogram", type = str, nargs='+',    metavar="X-metastatement [,Y-metastatement] [,title=figtitle] [,Xtitle=x-axel [,Ytitle=y-axel]] [,args]",        help = "Plots a 1D or 2D histogram, multiple plots with '|'.")
    
    arger.add_argument("-gm2", "--getmol2", type = str, nargs='+', metavar='pathto.mol2,column,metaname',         help = "Reads atom block column data from mol2-file and adds it to sdf-file as metadata.")#meta column path
    arger.add_argument("-pm2", "--putmol2", type = str, nargs='+',     metavar='input.mol2,output.mol2, column, metastatement, default, precision',        help = "Injects meta-data from sdf-file and adds it to mol2-file as atom block column data.")#metaname, column, path, defaultValue, precision, outpath
    
    arger.add_argument("-sbm", "--stripbutmeta", type = str, nargs='+', metavar='statement', help = "Removes all atoms from molecules, except for those in given logical statement.")
    arger.add_argument("-ig", "--ignores", type = str, nargs='*', metavar='Element', default=['H'], help = "Ignores given atoms in distance calculations, etc. Default is H.")
    
    args = arger.parse_args()
    
    if not (args.input or args.config):
        arger.print_help()
        sys.exit(1)
    
    manyfiles = args.input
    
    
    def main(inputfile, options):
        options['input'] = inputfile
        onerun = Runner(options)
        onerun.runOptions()
        
        if onerun.writer==None:
            onerun.sdf.writer(onerun.writetype,**onerun.wriarg)
            
        
    
    for onefile in manyfiles:
        options = vars(args)
        deli = []
        for key in options.keys():
            if not options[key] and type(options[key])!=int:
                deli.append(key)
        options = dict([(key, options[key]) for key in options if not key in deli])
        del(deli)
        main(onefile, dict(options))
    
    if 'plt' in globals():
        onefile.plt.show()