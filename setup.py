from setuptools import setup, find_packages

# The version lives in TeloBP/_version.py so the package and its run metadata report the same one
version = {}
with open('TeloBP/_version.py') as f:
    exec(f.read(), version)

setup(
    name='TeloBP',
    version=version['__version__'],
    packages=find_packages(),
    install_requires=[
        'biopython',
        'matplotlib',
        'matplotlib-inline',
        'numpy',
        'pandas',
        'pyparsing',
        'zipp'
    ],
    # install_requires=[
    #     'biopython==1.81',
    #     'matplotlib==3.7.3',
    #     'matplotlib-inline==0.1.6',
    #     'numpy==1.26.0',
    #     'pandas==2.2.0',
    #     'pyparsing==3.1.1',
    #     'zipp==3.16.2'
    # ],
)
