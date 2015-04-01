# -*- coding: latin-1 -*-
from distutils.core import setup
setup(name='sdfconf.py',
      version='0.753x',
      install_requires = ['numpy','matplotlib','argparse'],
      author = 'Sakari Lätti',
      author_email = 'sakari.latti@jyu.fi',
      package_dir = {'':'src/sdfconf'},
      #py_modules=['./src/sdfconf/sdfconf.py'],
      py_modules=['sdfconf'],
      #scripts=['./src/sdfconf/sdfconf.py'],
      scripts=['sdfconf'],
      )