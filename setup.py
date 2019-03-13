from distutils.core import setup
from distutils.extension import Extension

setup(name='rawor',
      version='0.1',
      description='Analytic calculation of the DDR, DRR, and RRR terms in the galaxy three point correlation function for uniform periodic cubes.',
      url='https://github.com/dpearson1983/rawor',
      author='David W. Pearson',
      author_email='dpearson1983@gmail.com',
      license='MIT',
      ext_modules=[
          Extension("rawor", ["rawor/pyRawor.cpp"],
          libraries=["boost_python3", "boost_numpy3"])
      ],
      classifiers=[
          "Programming Language :: Python :: 3",
          "Operating System :: Linux"
      ],)
                   
