#from distutils.core import setup
from setuptools import setup

setup(name="suchyta_utils", 
      version="0.1.0",
      description="Utilities written by Eric Suchyta",
      license = "GPL",
      packages=['suchyta_utils'],
      package_data={'suchyta_utils':['custom-sytle.mpl']},
      scripts=['bin/.screenrc-nest_inner','bin/.screenrc-nest_outer','bin/screen-split','bin/screen-retach','bin/screen-detach','bin/catpdfs'],
      author="Eric Suchyta",
      author_email="eric.d.suchyta@gmail.com")
