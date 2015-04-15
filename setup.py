# -*- coding: latin-1 -*-
from distutils.core import setup
setup(name='sdfconf',
      version='0.804',
      requires = ['numpy','matplotlib','argparse'],
      author = 'Sakari Lätti',
      author_email = 'sakari.latti@jyu.fi',
      packages = ['sdfconf'],
      package_dir = {'sdfconf':'src/sdfconf'},
      #package_dir = { '':'src/sdfconf'},
      #py_modules=['sdfconf.sdf', 'sdfconf.mol2', 'sdfconf.functions', 'sdfconf.runner', 'sdfconf.findable'],
      #py_modules=['sdfconf.sdfconf', 'sdfconf.runner'],
      scripts=['bin/sdfconf'],
      #scripts=['src/sdfconf/runner.py'],
      )
