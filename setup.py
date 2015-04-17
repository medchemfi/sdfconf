# -*- coding: latin-1 -*-

from setuptools import setup, find_packages
#from distutils.core import setup

setup(name =                'sdfconf',
      #version='0.804',
      version =             '0.8.1.0',
      description =         ("Diverse manipulation and alysis tool for .sdf files."),
      long_description =    read('README'),
      #requires = ['numpy','matplotlib','argparse'],
      install_requires =    ['numpy>=1.7.1','matplotlib>=1.4.2'],
      author =              'Sakari Lätti',
      author_email =        'sakari.latti@jyu.fi',
      maintainer =          'Sakari Lätti',
      maintainer_email =    'sakari.latti@jyu.fi',
      packages =            ['sdfconf'],
      package_dir =         {'sdfconf':'src/sdfconf'},
      #scripts=['bin/sdfconf'],
      keywords =            'sdf mol2 conformation analyze histogram',
      url =                 'http://users.jyu.fi/~pentikai/',
      licence =             'MIT/expat', 
      entry_points = 
                            {'console_scripts': ['sdfconf = sdfconf.runner'], 
                             'setuptools.installation': ['eggsecutable = sdfconf.runner',], 
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
      )