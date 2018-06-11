#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import time
#from kitchen.iterutils import isiterable
#import os
#from collections import OrderedDict as OrDi


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

try:
    import functions
    from sdf import Sdffile
    from _version import __version__
except ImportError:
    from sdfconf import functions
    from sdfconf.sdf import Sdffile
    from sdfconf._version import __version__



#__version__ = ""
#execfile(os.path.dirname(os.path.realpath(__file__))+'/src/sdfconf/version.py')

#__version__ = open(os.path.join(os.path.dirname(__file__), 'VERSION')).read()

class Runner(object):
    
    #Order of excecution. Also includes both long and short names.
    order = ((
             ('verbose',        'v'     ), 
             ('nospam',         'ns'    ), 
             ('input',          'in'    ), 
             
             ('molgrouper',     'mog'   ), 
             ('confgrouper',    'cog'   ), 
             
             ('conftometa',     'ctf'   ), 
             ('conftoname',     'ctn'   ), 
             ('nametometa',     'ntm'   ), 
             ('removeconfname', 'rcn'   ), 
             ('removeconfmeta', 'rcm'   ), 
             ('ignores',        'ig'    ), 
             ('cut',            'cu'    ), 
             ('allcut',         'acu'   ), 
             ('combine',        'co'    ), 
             ('allcombine',     'aco'   ), 
             ('addcsv','csv'), 
             ('addatomiccsv',   'acsv'  ), 
             ('getmol2',        'gm2'   ), 
             ('addescape',      'aesc'  ), 
             ('addinside',      'ains'  ), 
             ('config',         'con'   ), 
             ('closestatoms',   'ca'    ), 
             ('closeratoms',    'cla'   ), 
             ('changemeta',     'cm'    ), 
             ('makenewmeta',    'mnm'   ), 
             ('stripbutmeta',   'sbm'   ), 
             ('proportion',     'pro'   ), 
             ('sortorder',      'so'    ), 
             #('proportion','pro'), 
             ('extract',        'ex'    ), 
             ('metatoname',     'mtn'   ), 
             ('removemeta',     'rm'    ), 
             ('pickmeta',       'pm'    ), 
             ('putmol2',        'pm2'   ), 
             ('histogram',      'hg'    ), 
             ('scatter',        'sca'   ), 
             ('getcsv',         'gc'    ), 
             ('getatomcsv',     'gac'   ), 
             ('metalist',       'ml'    ), 
             ('counts',         'nm'    ), 
             ('donotprint',     'dnp'   ), 
             ('donotplot',      'dnpl'  ), 
             ('split',          's'     ), 
             ('makefolder',     'mf'    ), 
             ('output',         'out'   ), 
             ('overwrite',      'o'     ), 
             
             ('sdf',            'sdf'   ), #not in argparse. Should convert writetype back to sdf
             
             ) )
    
    '''
    simplelist =    ('conftometa', 'conftoname', 'removeconfname', 'removeconfmeta', 
                     'ignores', 'nametometa', 'metatoname', 'removemeta', 
                     #'pickmeta', 'input', 'proportion' , 'histogram', 'verbose', 
                     'pickmeta', 'input', 'proportion' , 'verbose', ##CHANGED: removed histogram
                     )
    
    simpleloops =   ('getmol2', 'closestatoms', 'closeratoms', 'changemeta', #'mergemeta',
                     'sortorder', 'stripbutmeta', 'extract','makenewmeta', 
                     'config', 'cut', 'allcut', 'combine', 'allcombine', 
                     'addcsv', 'addatomiccsv', 'putmol2', 'addescape', 
                     'addinside',
                     #'histogram' ##CHANGED: added histogram 
                     )
    '''
    
    
    singulars =    ('overwrite', 'verbose', 'conftometa', 'conftoname', 'removeconfname',
                     #'removeconfmeta', 'makefolder', 'metalist', 'donotprint',
                     'removeconfmeta', 'metalist', 'donotprint',
                     'molgrouper', 'confgrouper', #new arrivals 24.2.2017
                     'ignores', 'donotplot', 'nospam'
                     )
    
    
    #listbatch = ()
    
    writers =       ('overwrite', 'output', 'stdout') #default stdout
    writetypes =    {'getcsv':True,'getatomcsv':True,'metalist':True,'counts':True,'donotprint':True,'sdf':True,'split':False,'makefolder':False} #default 'sdf'
    parametricwrite =   ('getcsv', 'getatomcsv', 
                         'split', 'makefolder') #NEW 4.4.2017
    graphers =      ('histogram', 'scatter')
    
    
    initials = singulars + graphers +  ('extract', 'addescape', 'addinside', 'closestatoms', 'closeratoms', ) 
    
    def __init__(self, options=dict()):
        
        self.options = options if options else dict()
        
        for option, value in (('verbose', False),('proportion', []),('ignores',[])):
            if not option in self.options:
                self.options[option] = value
        
        self.times = [time.time()]
        self.inpath = None
        self.wriarg = dict()
        self.writetype = 'sdf'
        self.writer = None
        self.propor = []
        
        self.setVerbose()
        self.setIgnores()
        
        self.sdf = Sdffile()
        
        self._noPlot = False
        
    def setPropor(self,parameter=None):
        if not parameter or parameter is True:
            self.propor = []
        else:
            #self.propor = parameter.strip()
            self.propor.append(parameter)
    
    def setVerbose(self, verbose=False):
        if isinstance(verbose, (list, tuple)):
            self.verbose = verbose[0]
        else:
            self.verbose = verbose
    
    def setNoPlot(self, noPlot=False):
        if isinstance(noPlot, (list, tuple)):
            self._noPlot = noPlot[0]
        else:
            self._noPlot = noPlot
    
    def setIgnores(self,ignore=['H']):
        self.ignore = ignore
        try:
            self.sdf.setIgnores(self.ignore)
        except AttributeError:
            pass
    
    
    def runOptions(self):
        #after config, do the rest anyway.
        #print self.options
        for optionpair in Runner.order:
            if optionpair[0] in self.options or optionpair[1] in self.options:
                self.times.append(time.time())
                self.funcselector(optionpair[0], self.options[optionpair[0]], False)
        
        if self.verbose:
            print('Run complete, it took {:.3f} seconds.'.format(self.times[-1]-self.times[0]))
        
    def runConfig(self,confpath):
        with open(confpath, 'r') as confile:
            configures = (Runner.parsecon(confline)  for confline in confile)
            #print(list(configures))
            for option in configures:
                #print option
                if option is None:
                    continue
                #elif option[0] in Runner.order: #self.options:
                elif option[0] in zip(*Runner.order)[0]: #self.options:
                    self.times.append(time.time())
                    self.funcselector(option[0], option[1], True)
        
    
    def taskLib(self,task,option,param=None,**kwargs):
        steps = kwargs.get('steps',1)
        timedif = lambda : '{:.3f}'.format(self.times[-1]-self.times[-(1+steps)])
        
        NoLamb = lambda : None
        def extMes(n):
            '''
            mypropo = ' ' + self.sdf.propomes(self.propor)
            mypropo = mypropo.replace('{','{{')
            mypropo = mypropo.replace('}','}}')
            '''
            mypropo = []
            for pro in self.propor:
                txt  = '  ' + self.sdf.propomes(pro)
                txt = txt.replace('{','{{')
                txt = txt.replace('}','}}')
                mypropo.append(txt)
            
            messes = (
                      #lambda : ('\n'.join((' After logical chop {}, sdf-file has {} molecules and {} conformations left.', mypropo)).strip('\n') ,(param,len(self.sdf._dictomoles),len(self.sdf))), 
                      #lambda : ('\n'.join((' After logical chop {}, sdf-file has {} molecules and {} conformations left.', mypropo)).rstrip() ,(param,len(self.sdf._dictomoles),len(self.sdf))), 
                      #lambda : ('\n'.join(('Initially sdf-file has {} molecules and {} conformations.',mypropo)).strip(),(len(self.sdf._dictomoles),len(self.sdf))), 
                      lambda : ('\n'.join([' After logical chop {}, sdf-file has {} molecules and {} conformations left.',] + mypropo).rstrip() ,(param,len(self.sdf._dictomoles),len(self.sdf))), 
                      lambda : ('\n'.join(['Initially sdf-file has {} molecules and {} conformations.',] + mypropo).strip(),(len(self.sdf._dictomoles),len(self.sdf))), 
                      lambda : (' Chopping complete, it took {} seconds.',(timedif(),)) 
                      )
            return messes[n]
            
        tasks =   {
                    'conftometa'    :lambda i : (
                        self.sdf.addConfs,
                        #lambda : ((False,True),),
                        lambda : (False,True,),
                        NoLamb,
                        NoLamb,
                        lambda : ('Conformation numbers added to metadata. It took {} seconds.', (timedif(),))
                        )[i], 
                    'conftoname'     :lambda i : (
                        self.sdf.addConfs,
                        #lambda : ((True,False),),
                        lambda : (True,False,),
                        NoLamb,
                        NoLamb,
                        lambda : ('Conformation numbers added to names. It took {} seconds.', (timedif(),))
                        )[i], 
                    'removeconfname' :lambda i : (
                        self.sdf.removeConfs,
                        #lambda : ((True,False),),
                        lambda : (True,False,),
                        NoLamb,
                        NoLamb,
                        lambda : ('Conformation numbers removed from names. It took {} seconds.',(timedif(),))
                        )[i], 
                    'removeconfmeta' :lambda i : (
                        self.sdf.removeConfs,
                        #lambda : ((True,False),),
                        lambda : (True,False,),
                        NoLamb,
                        NoLamb,
                        lambda : ('Conformation numbers removed from metafield \'confnum\'. It took {} seconds.',(timedif(),))
                        )[i], 
                    'nametometa'     :lambda i : (
                        self.sdf.nameToMeta,
                        lambda : (param,),
                        NoLamb,
                        NoLamb,
                        lambda : ('Name written to metafieldfield {}. It took {} seconds.',(param,timedif()))
                        )[i], 
                    'metatoname'     :lambda i : (
                        self.sdf.metaToName,
                        lambda : (param,),
                        NoLamb,
                        NoLamb,
                        lambda : ('Metafield written to name {}. It took {} seconds.',(param,timedif()))
                        )[i],
                    'removemeta'     :lambda i : (
                        self.sdf.removemeta,
                        lambda : (param, False),
                        lambda : ('Metafieldfield(s) in ({}) Removed. It took {} seconds.',(param,timedif())),
                        NoLamb,
                        NoLamb
                        )[i], 
                    'pickmeta'       :lambda i : (
                        self.sdf.removemeta,
                        lambda : (param, True),
                        lambda : ('All metafields except those in ({}) Removed. It took {} seconds.',(param,timedif())),
                        NoLamb,
                        NoLamb
                        )[i], 
                    'nospam'    :lambda i : (
                        self.sdf.setSpam,
                        #lambda : ((False,True),),
                        lambda : ('',),
                        NoLamb,
                        NoLamb,
                        #NoLamb,
                        lambda : ('Spam removed from comments. Please remeber to cite sdfconf in your work!', () )
                        )[i], 
                    'getmol2'        :lambda i : (
                        self.sdf.getMol2DataStr,
                        lambda : (param,),
                        NoLamb,
                        NoLamb,
                        NoLamb
                        )[i], 
                    'closestatoms'   :lambda i : (
                        self.sdf.closestStr,
                        lambda : (param,),
                        lambda : (' Calculated distances to atoms from point {}.',(param, )),
                        lambda : ('Calculate point of interest distances.',()), 
                        lambda : (' Distances calculated, it took {} seconds.',(timedif(),))
                        )[i], 
                    'closeratoms'    :lambda i : (
                        lambda x: self.sdf.closer( *functions.splitter(x)[:2] ),
                        lambda : (param,),
                        lambda : (' Calculated number of atoms closer: {}.',(param, )),
                        lambda : ('Calculate number of atoms closer to point of interest.',()), 
                        lambda : (' Numbers of atoms calculated, it took {} seconds.',(timedif(),))
                        )[i], 
                    'changemeta'     :lambda i : (
                        lambda x: self.sdf.changeMetaName(*[item.strip() for item in x.split('>')]),
                        lambda : (param,),
                        lambda : ('Change {} done.',(param,)),
                        NoLamb,
                        NoLamb
                        )[i], 
                    'sortorder'      :lambda i : (
                        self.sdf.sortMe,
                        lambda : (param,),
                        lambda : ('Sort {} done.',(param,)),
                        NoLamb,
                        NoLamb
                        )[i], 
                    'stripbutmeta'   :lambda i : (
                        self.sdf.stripbutmeta,
                        lambda : (param,),
                        lambda : ('All atoms, execpt for those in statement {} removed!',(param,)),
                        NoLamb,
                        NoLamb
                        )[i], 
                    'extract'        :lambda i : (
                        self.sdf.mollogicparse,
                        lambda : (param,),
                        extMes(0),
                        extMes(1),
                        extMes(2)
                        )[i],
                    'makenewmeta'    :lambda i : (
                        self.sdf.makeNewMetaStr,
                        lambda : (param,),
                        lambda : (' New metafield {} made.',(param,)),
                        lambda : ('Make new metafields.',()),
                        lambda : (' Making new metafields done. It took {} seconds.',(timedif(),))
                        )[i], 
                    'addcsv'         :lambda i : (
                        self.sdf.addCsvMeta,
                        lambda : (param, ),
                        lambda : ('Metadata from csv-file {} added. It took {} seconds.',(param,timedif(),)),
                        NoLamb,
                        NoLamb
                        )[i],
                    'addatomiccsv'   :lambda i : (
                        self.sdf.addAtomicCsvMeta,
                        lambda : (param, ),
                        lambda : ('Dictionary metadata from atomic csv-file {} added. It took {} seconds.',(param,timedif(),)),
                        NoLamb,
                        NoLamb
                        )[i], 
                    #'addcsv'         :lambda i : (self.sdf.addCsvMeta, lambda : param, lambda : ('Metadata from csv-file {} added. It took {} seconds.',(param,timedif(),)),NoLamb,NoLamb)[i],
                    #'addatomiccsv'   :lambda i : (self.sdf.addAtomicCsvMeta, lambda : param, lambda : ('Dictionary metadata from atomic csv-file {} added. It took {} seconds.',(param,timedif(),)),NoLamb,NoLamb)[i], 
                    'input'          :lambda i : (
                        self.sdf.xreadself,
                        lambda : (param,),
                        NoLamb,
                        lambda : ('Starting to read file {}',(param,)),
                        lambda : (' Reading file done. It took {} seconds.', (timedif(),))
                        )[i],
                    'verbose'        :lambda i : (
                        self.setVerbose,
                        lambda : (param,),
                        NoLamb,
                        NoLamb,
                        lambda : ('Verbose enabled.' if param else 'Verbose disabled.',())
                        )[i], 
                    'ignores'        :lambda i : (
                        self.setIgnores,
                        lambda : (param,),
                        NoLamb,
                        NoLamb,
                        #lambda : ('Ignores set to {}.',(param,))
                        lambda : ('Ignores set to {}.', param)
                        )[i], 
                    'config'         :lambda i : (
                        self.runConfig,
                        lambda : (param,),
                        lambda : (' Config-file {} done.',(param,)),
                        lambda : ('Run config-files.',()),
                        lambda : (' Running config-files done',())
                        )[i], 
                    'histogram'      :lambda i : (
                        self.sdf.histogramFromListOfStrings,
                        lambda : (param,),
                        NoLamb,
                        lambda : ('Start plotting histograms',()),
                        lambda : (' Plotting done',())
                        )[i], 
                    'scatter'        :lambda i : (
                        self.sdf.scatterFromListOfStrings,
                        lambda : (param,),
                        NoLamb,
                        lambda : ('Start plotting scatter plots',()),
                        lambda : (' Plotting done',())
                        )[i], #FIXME
                    'donotplot'       :lambda i : (
                        self.setNoPlot,
                        lambda : (param,),
                        NoLamb,
                        NoLamb,
                        lambda : ('No plotting enabled.' if param else 'No plotting disabled.',())
                        )[i], #FIXME
                    'cut'            :lambda i : (
                        #lambda x :self.sdf.listremove(Sdffile(x.strip()),True),
                        lambda x :self.sdf.listremove(x.strip(),True),
                        lambda : (param,),
                        lambda : ('Removing all molecules (matching name and confnum) present in {} complete. It took {} seconds.', (param, timedif())),
                        NoLamb,
                        NoLamb
                        )[i], 
                    'allcut'         :lambda i : (
                        #lambda x :self.sdf.listremove(Sdffile(x.strip()),False),
                        lambda x :self.sdf.listremove(x.strip(),False),
                        lambda : (param,),
                        lambda : ('Removing all molecules (matching name) present in {} complete. It took {} seconds.', (param, timedif())),
                        NoLamb,
                        NoLamb
                        )[i], 
                    'combine'        :lambda i : (
                        #lambda x :self.sdf.sdfMetaCombi(Sdffile(x.partition(';')[0].strip()),(True,True),x.partition(';')[2].strip().lower() in ('o','over','overwrite')), 
                        lambda x :self.sdf.sdfMetaCombi(Sdffile(x.partition(';')[0].strip()), True, True, x.partition(';')[2].strip().lower() in ('o','over','overwrite')), 
                        lambda : (param,), 
                        lambda : ('Combining metadata from {} (matching name and confnum) complete. It took {} seconds.', (param, timedif())), 
                        NoLamb, 
                        NoLamb
                        )[i], 
                    'allcombine'     :lambda i : (
                        #lambda x :self.sdf.sdfMetaCombi(Sdffile(x.partition(';')[0].strip()),(True,False),x.partition(';')[2].strip().lower() in ('o','over','overwrite')), 
                        lambda x :self.sdf.sdfMetaCombi(Sdffile(x.partition(';')[0].strip()), True, False, x.partition(';')[2].strip().lower() in ('o','over','overwrite')), 
                        lambda : (param,), 
                        lambda : ('Combining metadata from {} (matching name) complete. It took {} seconds.', (param, timedif())), 
                        NoLamb, 
                        NoLamb
                        )[i], 
                    'proportion'     :lambda i : (
                        self.setPropor,
                        lambda : (param,),
                        NoLamb,
                        NoLamb, 
                        (lambda : ('Propor set to {}.',(param,))) if param else NoLamb 
                        )[i], 
                    #'putmol2'        :lambda i : (lambda x :self.sdf.injectMol2DataStr(x.strip()), lambda : (param,), lambda : ('Injection "{}" done.', (param,)) , NoLamb, NoLamb)[i],
                    'putmol2'        :lambda i : (
                        self.sdf.injectMol2DataStr,
                        lambda : (param,),
                        lambda : ('Injection "{}" done.', (param,)) ,
                        NoLamb,
                        NoLamb
                        )[i],
                    'addescape'      :lambda i : (
                        self.sdf.escapeStr,
                        lambda : (param, False),
                        lambda : (' Escape number from {} added. It took {} seconds.',(param, timedif())),
                        lambda : ('Calculating escape numbers...', ()),
                        NoLamb
                        )[i],
                    'addinside'      :lambda i : (
                        self.sdf.escapeStr,
                        lambda : (param, True),
                        lambda : (' Inside number from {} added. It took {} seconds.',(param, timedif())),
                        lambda : ('Calculating inside numbers...', ()),
                        NoLamb
                        )[i],
                    'molgrouper'      :lambda i : (
                        self.sdf.setGrouper,
                        lambda : (param, None),
                        NoLamb,
                        NoLamb,
                        lambda : ('Molgrouper changed to {}', (param,),()),
                        )[i],
                    'confgrouper'      :lambda i : (
                        self.sdf.setGrouper,
                        lambda : (self.sdf.molgrouper, param),
                        NoLamb,
                        NoLamb,
                        lambda : ('Confgrouper set to {}.',(param,)),
                        )[i],
                   }
        selector = {'func':0,'loop':2,'initial':3,'final':4}
        
        try:
            if task == 'func':
                mytask=tasks.get(option)
                fu = mytask(0)
                para = mytask(1)()
            elif option in tasks:#else:
                #print('**task: {}, option: {}'.format(task, option))
                fu=tasks.get(option)(selector.get(task))
                para=()
            else:
                fu = None
        except KeyError:
            raise KeyError('Wrong task={} or option={}'.format(task,option))
        if fu:
            retu = fu(*para)
        else:
            retu = None
        if task == 'func':
            self.times.append(time.time())
        return retu
    
    #def messenger(self, task, mypara=params,**kwargs):
    def messenger(self, task, option, params,**kwargs):
        if self.verbose:
            #mes = self.taskLib(task, option, mypara, **kwargs)
            mes = self.taskLib(task, option, params, **kwargs)
            if mes:
                print(mes[0].format(*mes[1]))
    
        
    
    def funcselector(self,option,params, fromconfig=False):
        '''
        Takes an option handle and its parameters and executes them accordingly.
        Parameters may include list of parameters.
        '''
        
        #print('option:{}; params:{}, fromconfig:{}'.format(option,params,fromconfig))
        
        #elif option in Runner.simpleloops or option in Runner.simplelist:
        if option in Runner.singulars:
            #if not fromconfig and option in Runner.simplelist:
            #if  option in Runner.simplelist: ##TODO
                
            params = (params,)
                
            #self.messenger('initial', option, params)
            #print('option:{}; params:{}, fromconfig:{}'.format(option,params,fromconfig))
        
        if option in Runner.initials:
            self.messenger('initial', option, params)
        
        if option == 'input':
            self.messenger('initial', option, params)
            self.inpath = params
            self.taskLib('func', option, params)
            self.messenger('final', option, params)
        
        elif option in Runner.writetypes:
            #if Runner.writetypes[option]:
            if Runner.writetypes[option]:
                self.writetype = option
            if option in Runner.parametricwrite:
                self.wriarg[option] = params[0]
                #if option == 'output':
                
            else:
                self.wriarg[option] = params
            #self.wriarg[option]=params
            
        elif option in Runner.writers:
            writers = {'overwrite':lambda x : (self.inpath,), 'output': lambda x : x, 'stdout': lambda x: None} #default stdout
            #try:
            self.wriarg['path'] = writers.get(option, None)(params)[0]
            
            #except TypeError:
                #pass
            self.writer = option
            self.sdf.writer(self.writetype,**self.wriarg)
            self.times.append(time.time())
            if self.verbose:
                print('File writing done, it took {} seconds.'.format(self.times[-1]-self.times[-2]))
        elif option in Runner.graphers:
            self.taskLib('func', option, params)
            self.messenger('loop', option, params)
            self.messenger('final', option, params)
            
        else:
            #elif option in Runner.taskLi  : #tasks
            for oneparam in params:
                self.taskLib('func', option, oneparam)
                self.messenger('loop', option, oneparam)
            #self.messenger('final', option, oneparam, steps=len(params))
            self.messenger('final', option, params, steps=len(params))
            
        
    @staticmethod
    def parsecon(oneline):
        '''
        Parse a line from config file into a proper option-parameter-pair
        '''
        params = [para.strip() for para in oneline.partition('#')[0].partition('::')]
        if params[0] == '':
            return None
        '''
        if params[0] not in zip(*Runner.order):
            for longname in Runner.order:
                if params[0] == Runner.order[longname]:
                    params[0] = longname
                    break
        '''
        for optionpair in Runner.order:
            if params[0] in optionpair:
            #if params[0] == Runner.order[longname]:
                params[0] = optionpair[0]
                break
        if params[1] == '::' and params[2] != '':
            return (params[0], [functions.numify(cell.strip()) for cell in params[2].split(';;')])
        else:
            return (params[0], (True,))

    
#End of Runner

    
def main(arguments=None):
    #main should only collect arguments...
    
    #arger = argparse.ArgumentParser(description='******************\nSome bad-ass manipulation of SDfiles. Also data retrieval/injection for .mol2-files. \n\nVersion {0} \n\nIf you publish work using sdfconf, please cite:\nmanuscript\n\nNotice that documentation is not completely up to date. \n******************\n'.format(__version__), formatter_class=argparse.RawTextHelpFormatter)
    arger = argparse.ArgumentParser(description='******************\nSome bad-ass manipulation of SDfiles. Also data retrieval/injection for .mol2-files. \n\nVersion {0} \n\nIf you publish work using sdfconf, please cite:\nmanuscript\n\nNotice that documentation is not completely up to date. \n******************\n'.format(__version__),  formatter_class=argparse.RawDescriptionHelpFormatter)
    
    arger.add_argument("input", metavar = 'input.sdf', nargs='*', type = str, default = None,  help="Specify the input  SDfile")
    
    choicewrite = arger.add_mutually_exclusive_group()
    choicewrite.add_argument("-out", "--output", nargs=1, metavar='output.file', type = str, default=None,  help = "Specify output file. It may be sdf or csv, depending on other arguments.")
    choicewrite.add_argument("-o", "--overwrite", action='store_true',      help = "Overwrite to input file. You don't need to specify output file")
    
    arger.add_argument("-v", "--verbose", action = "store_true" ,      help = "More info on your run.")
    
    arger.add_argument("-con","--config", metavar = 'config.txt', nargs='+', type = str, help="Specify the config file. Config file includes lines of argument name, followed by '::' and argument value. Separate multiple values with ';;'.")
    arger.add_argument("-ctm", "--conftometa", action = "store_true",           help = "Add conformation number to metafielf 'confnum'. If number in name doesn't exist, makes a new one.")
    
    arger.add_argument("-ctn", "--conftoname",  action = "store_true",           help = "add conformation number to name from metafield 'confnum'. If number in metafield doesn't exist, make a new one.")
    arger.add_argument("-mtn", "--metatoname", nargs=1, metavar='meta' , type = str,                 help = "Change the name of molecules to the data in given metafield.")
    arger.add_argument("-ntm", "--nametometa", nargs=1,  metavar='meta', type = str,                 help = "Copy the name of molecules into given metafield.")
    arger.add_argument("-rcn", "--removeconfname",  action="store_true",    help = "remove conformation number from name.")
    arger.add_argument("-rcm", "--removeconfmeta",  action="store_true",    help = "remove conformation number from metafield 'confnum'.")
    
    arger.add_argument("-aesc", "--addescape",  type=str, nargs='+', metavar="file,r[,optional args]", help = "Add metafield escapenum which is a list of atoms not in range of atoms in other molecule (file). Optional args include mol, max, and name. mol specifies which molecule in a file you want to compare to, starts from 0. With option max quits if max atoms are found, 0 means all. N and M default to 0.")
    arger.add_argument("-ains", "--addinside",  type=str, nargs='+', metavar="file,r[,optional args]", help = "Add metafield insidenum which is a list of atoms in range of atoms in other molecule (file). Optional args include mol, max, and name. mol specifies which molecule in a file you want to compare to, starts from 0. With option max quits if max atoms are found, 0 means all. N and M default to 0.")
    #"File,range,name[,name=some][,mol=N][,max=M]"
    arger.add_argument("-co", "--combine", metavar='addition.sdf', type = str, nargs='+',     help = "Combine metadata from specified file to the data of original file. Confromation numbers must match.")
    arger.add_argument("-aco", "--allcombine", metavar='addition.sdf', type = str, nargs='+', help = "Combine metadata from specified file to the data of original file. Names must match.")
    
    arger.add_argument("-cu", "--cut", metavar='unwanted.sdf', type = str, nargs='+',           help = "Remove molecules in specified file from original file. Confromations must match.")
    arger.add_argument("-acu", "--allcut", metavar='unwanted.sdf', type = str, nargs='+',       help = "Remove molecules in specified file from original file. Names must match. Not tested.")
    
    arger.add_argument("-csv", "--addcsv", metavar='path [,molcol=<column>] [,confkey=<column>]', type = str, nargs='+',           help = "Add metadata from csv-file. File must have a 1-line header, it gives names to metafields. By default reads molecule names from column 0. If name includes confnumber, meta is only added molecules with same confnumber. Columns including molecule names and conformation numbers may also be specified with either column number or header name.")
    arger.add_argument("-acsv", "--addatomiccsv", metavar='path[,molcol=<column>] [,confkey=<column>] [,atomnnumber=<column>]', type = str, nargs='+',    help = "Add atomic metadata from csv-file. File must have a 1-line header, it gives names to metafields. By default reads molecule names from column 0 and atom numbers from column 'atomn_number'. If name includes confnumber, meta is only added molecules with same confnumber. Columns including molecule names, conformation numbers and atom numbers may also be specified with either column number or header name.")
    
    arger.add_argument("-ex", "--extract", metavar='statement', type = str, nargs='+',          help = "Pick or remove molecules from file by metafield info. Either with logical comparison or fraction of molecules with same name. Closest_atoms{:5}==soms, 2.5>Closest_atoms(soms)[], Closest_atoms[:3]<5.5, ID=='benzene'.")
    arger.add_argument("-pro", "--proportion", metavar='statement', nargs='+', type = str,                 help = "Takes one exctract-like metastatement and prints proportion of molecules and conformations fulfilling it after every chop, if you are on verbose mode.")
    
    choiceremo = arger.add_mutually_exclusive_group()
    choiceremo.add_argument("-rm", "--removemeta", nargs='+', metavar='unwanted', type = str,              help = "Remove metadata from molecules. Takes multiple values, separaterd by comma(,) or semicolon(;). If first is '?', means 'all but'.")
    choiceremo.add_argument("-pm", "--pickmeta", nargs='+', metavar='wanted', type = str,                help = "Remove all nonspecified metadata from molecules. Takes multiple values, separaterd by comma(,) or semicolon(;). If first is '?', means 'all but'.")
    
    arger.add_argument("-s", "--split", type = str, nargs=1,                         help = "Split the file into even pieces. Positive means number of files while negative means number of molecules per file. 0 doesn't apply.")
    arger.add_argument("-mf", "--makefolder", type = str, nargs=1,        help = "Put outputfile(s) into folder.")
    #arger.add_argument("-mf", "--makefolder", action = "store_true",        help = "Put outputfile(s) into folder(s).")
    
    outputtype = arger.add_mutually_exclusive_group()
    outputtype.add_argument("-gc",  "--getcsv",     type = str, nargs=1 ,                 help = "Writes a .csv-file istead of .sdf-file. Specify which fields you'll need, separated by ','.")
    outputtype.add_argument("-gac", "--getatomcsv", type = str, nargs=1 ,            help = "Writes a .csv-file istead of .sdf-file. Specify which fields you'll need, separated by ','. Writes dictionaries to separate lines.")
    outputtype.add_argument("-ml",  "--metalist",   action = "store_true",     help = "Writes a list of metafields.")
    outputtype.add_argument("-nm",  "--counts",     type = int, nargs='?', const=0, choices=(0,1,2),  help = "Number of different molecules and different conformations. 0=just sums, 1=by molecule name, 2=both.")
    outputtype.add_argument("-dnp", "--donotprint", action = "store_true",  help = "No output")
    arger.add_argument("-dnpl", "--donotplot", action = "store_true",  help = "No do print plots on screen. Saving still works")
    
    arger.add_argument("-ca",  "--closestatoms",    type = str, nargs='+',  metavar='(x, y, z) [,name] [,interests=value]', help = "Calculates the closest atoms (distances by atom number) to given point. Creates a metafield with given name, if no name is given 'Closest_atoms' is created. (xx, yy, zz) may be replaced by metastatement describing single atom number.")
    arger.add_argument("-cla", "--closeratoms",     type = str, nargs='+', metavar= "(x, y, z),meta", help = "Calculates number of atoms closer to the given point, than the ones given in meta. Adds metafields 'Closest_atom_from_{meta}' and 'Closer_atoms_than_{meta}'.")
    
    arger.add_argument("-mnm", "--makenewmeta",     type = str, nargs='+', metavar='newmeta=metastatement',     help = "Makes a new metafield based on metastatement.")
    arger.add_argument("-cm",  "--changemeta",      type = str, nargs='+',   metavar='olname1>newname1' , help = "Changes names of metafields. [olname1>newname1|oldname2>newname2].")
    arger.add_argument("-so",  "--sortorder",       type = str, nargs='+',    metavar='meta', help = "Sorts molecules of a file in order of metafield. <MolecularWeight|>Id Sorts molecules first ascending by weight, then descenting by name.")
    
    arger.add_argument("-hg",  "--histogram",       type = str, nargs='+',    metavar="X-metastatement [,Y-metastatement] [,title=figtitle] [,Xtitle=x-axel [,Ytitle=y-axel]] [,args]",        help = "Plots a 1D or 2D histogram. Multiple plots as separate strings.")
    arger.add_argument("-sca", "--scatter",         type = str, nargs='+',    metavar="X-metastatement ,Y-metastatement [,group-metastatement] [,title=figtitle] [,Xtitle=x-axel [,Ytitle=y-axel]] [,args]",        help = "Plots a scatter plot. Multiple plots as separate strings.")
    
    arger.add_argument("-gm2", "--getmol2",         type = str, nargs='+', metavar='pathto.mol2,column,metaname',         help = "Reads atom block column data from mol2-file and adds it to sdf-file as metadata.")#meta column path
    arger.add_argument("-pm2", "--putmol2",         type = str, nargs='+',     metavar='input.mol2,output.mol2, column, metastatement[, default[, precision]]',        help = "Injects meta-data from sdf-file and adds it to mol2-file as atom block column data.")#metaname, column, path, defaultValue, precision, outpath
    
    arger.add_argument("-sbm", "--stripbutmeta",    type = str, nargs='+', metavar='statement', help = "Removes all atoms from molecules, except for those in given logical statement.")
    arger.add_argument("-ig",  "--ignores",         type = str, nargs='*', metavar='Element', default=['H'], help = "Ignores given atoms in distance calculations, etc. Default is H.")
    
    #arger.add_argument("-mog", "--molgrouper",      type = str, nargs='?', metavar='metaname', const='', default='', help="Give existing metaname, by which molecules should be grouped. If you give empty, molecule name is used.")
    #arger.add_argument("-cog", "--confgrouper",      type = str, nargs='?', metavar='metaname', const='confnum', default='confnum', help="Give existing metaname, by which conformations will be identified. If you give empty, 'confnum' is used.")
    arger.add_argument("-mog", "--molgrouper",      type = str, nargs='?', metavar='metaname', const='', help="Give existing metaname, by which molecules should be grouped. If you give empty, molecule name is used.")
    arger.add_argument("-cog", "--confgrouper",     type = str, nargs='?', metavar='metaname', const='confnum', help="Give existing metaname, by which conformations will be identified. If you give empty, 'confnum' is used.")
    
    arger.add_argument("-ns" , "--nospam",          action = "store_true", help = "Don't modify the comment line of files. Please, remember to cite sdfconf in your work! Works only for new entries in conf files.")
        
    args = arger.parse_args(arguments) if arguments else arger.parse_args() 
    
    if not (args.input or args.config):
        print('You must specify either an sdf or config file.')
        arger.print_usage()
        sys.exit(1)
    
    manyfiles = args.input
    
    global plots 
    plots = None
    
    def run(inputfile, options):
        options['input'] = inputfile
        onerun = Runner(options)
        onerun.runOptions()
        
        if onerun.writer==None:
            onerun.sdf.writer(onerun.writetype,**onerun.wriarg)
        
        global plots
        if (not onerun._noPlot) and (not plots) and (onerun.sdf._plt):
            plots = onerun.sdf._plt
        
        #print onerun.sdf._ignores
        
    
    for onefile in manyfiles:
        options = vars(args)
        deli = []
        for key in options.keys(): 
            #if not options[key] and type(options[key])!=int: ##FIXME: what's this?
            #if not options[key] and not isinstance(type(options[key]), int): #parameter is None or False but not 0 =>> do not excecute handle. 
            #print('{}: {}'.format(key, options[key]))
            '''
            if (not options[key]) and (not isinstance(options[key], int)): #parameter is None or False but not 0 =>> do not excecute handle. 
                print('- {}: {}'.format(key, options[key]))
                deli.append(key)
            else:
                print('+ {}: {}'.format(key, options[key]))
            '''
            item = options[key]
            if (isinstance(item, bool) and not item) or (item is None):
                deli.append(key)
            '''
                print('- {}: {}'.format(key, options[key]))
            else:
                print('+ {}: {}'.format(key, options[key]))
            '''
            
        options = dict([(key, options[key]) for key in options if not key in deli])
        del(deli)
        #run(onefile, dict(options))
        run(onefile, options) #no need for cast
    
    '''
    if 'plt' in globals():
        onefile.plt.show()
    '''
    if plots:
        plots.show()
        
'''def main():
    run(sys.argv[1:])'''

if __name__ == "__main__":
    #print(os.getcwd())
    main()