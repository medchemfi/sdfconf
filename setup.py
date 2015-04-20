# -*- coding: latin-1 -*-

from setuptools import setup, find_packages
import re
import os

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

#verstrline = open("src/sdfconf/version.py", "rt").read()
with open("src/sdfconf/version.py", "rt") as vf:
    VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
    for line in vf:
        mo = re.search(VSRE, line, re.M)
        if mo:
            verstr = mo.group(1)
            break
    if not mo:
        raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

setup(name =                'sdfconf',
      version = verstr, 
      description =         ("Diverse manipulation and alysis tool for .sdf files."),
      long_description =    read('README'),
      install_requires =    ['numpy>=1.7.1','matplotlib>=1.4.2'],
      author =              'Sakari Lätti',
      author_email =        'sakari.latti@jyu.fi',
      maintainer =          'Sakari Lätti',
      maintainer_email =    'sakari.latti@jyu.fi',
      packages =            ['sdfconf'],
      package_dir =         {'sdfconf':'src/sdfconf'},
      keywords =            'sdf mol2 conformation analyze histogram',
      url =                 'http://users.jyu.fi/~pentikai/',
      license =             'MIT/expat', 
      entry_points = 
                            {'console_scripts': ['sdfconf = sdfconf.runner:main'], 
                             'setuptools.installation': ['eggsecutable = sdfconf.runner:main',], 
                             }, 
      classifiers=           ['Development Status :: 4 - Beta',
                              'Environment :: Console',
                              'Intended Audience :: Science/Research',
                              'License :: OSI Approved :: MIT License',
                              'Natural Language :: English',
                              'Operating System :: OS Independent',
                              'Programming Language :: Python :: 2.7', 
                              'Programming Language :: Python :: 2 :: Only',
                              'Topic :: Scientific/Engineering :: Bio-Informatics', 
                              'Topic :: Scientific/Engineering :: Chemistry' ,
                              'Topic :: Software Development :: Libraries',
                              
                              ], 
       package_data =         {
                               '': ['VERSION'],
    }
      )