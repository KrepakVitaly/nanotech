from distutils.core import setup
import py2exe
import matplotlib

setup(
    name='nonlinear_cap',
    version='',
    packages=[''],
    url='',
    license='',
    author='user',
    author_email='',
    description='',
    console=['main.py'],
    data_files=matplotlib.get_py2exe_datafiles(),
    options={'py2exe': {'excludes': ['_gtkagg', '_tkagg'], }}
)
