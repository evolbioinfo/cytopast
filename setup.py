import os
from distutils.core import setup

setup(
    name='cytopast',
    packages=['cytopast'],
    package_data={'cytopast': [os.path.join('templates', '*.html'), os.path.join('templates', '*.js'),
                               os.path.join('examples', 'Albania', 'data',  'Albanian.tree.152tax.tre'),
                               os.path.join('examples', 'Albania', 'data',  'Annotation.Albanian.5chars.txt'),
                               os.path.join('examples', 'Albania', 'main.py'),
                               os.path.join('examples', 'C', 'data',  'HIV-1C.fasta'),
                               os.path.join('examples', 'C', 'data',  'HIV-1C.tab'),
                               os.path.join('examples', 'C', 'snakemake',  'Snakefile'),
                               os.path.join('examples', 'C', 'snakemake',  'config.yaml'),
                               os.path.join('..', 'README.md')]},
    long_description=open('README.md').read(),
    include_package_data=True,
    platform=['Linux'],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Software Development',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    version='0.1.1',
    description='Visualisation of PASTML trees.',
    author='Anna Zhukova',
    author_email='anna.zhukova@pasteur.fr',
    url='https://gitlab.pasteur.fr/phylo/cytopast',
    download_url='https://gitlab.pasteur.fr/phylo/cytopast/archive/0.1.zip',
    keywords=['PASTML', 'visualisation', 'phylogeny'],
    install_requires=['ete3', 'pandas', 'numpy', 'jinja2']
)
