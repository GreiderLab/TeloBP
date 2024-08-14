from setuptools import setup, find_packages

setup(
    name='TeloBP',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'biopython==1.81',
        'matplotlib==3.7.3',
        'matplotlib-inline==0.1.6',
        'numpy==1.26.0',
        'pandas==2.2.0',
        'pyparsing==3.1.1',
        'zipp==3.16.2'
    ],
)
