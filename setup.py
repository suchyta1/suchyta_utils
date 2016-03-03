#from distutils.core import setup
from setuptools import setup
import os

bindir = 'bin'
files = os.listdir(bindir)
binfiles = [os.path.join(bindir,f) for f in files]

setup(name="suchyta_utils", 
      version="0.1.01",
      description="Utilities written by Eric Suchyta",
      license = "GPL",
      packages=['suchyta_utils'],
      package_data={'suchyta_utils':['custom-sytle.mpl']},
      #scripts=['bin/.screenrc-nest_inner','bin/.screenrc-nest_outer','bin/screen-split','bin/screen-retach','bin/screen-detach','bin/catpdfs'],
      scripts=binfiles,
      author="Eric Suchyta",
      author_email="eric.d.suchyta@gmail.com")
