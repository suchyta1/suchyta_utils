#from distutils.core import setup
from setuptools import setup, find_packages

setup(name="suchyta_utils", 
      version="0.1.0",
      description="Utilities written by Eric Suchyta",
      license = "GPL",
      author="Eric Suchyta",
      author_email="eric.d.suchyta@gmail.com",
      packages=['suchyta_utils'],
      package_data={'plot':'custom-style.mpl'})
