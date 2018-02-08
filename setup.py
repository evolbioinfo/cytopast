import os
from distutils.core import setup

setup(
    name='cytopast',
    packages=['cytopast'],
    package_data={'cytopast': [os.path.join('templates', '*.html'), os.path.join('templates', '*.js'),
                               os.path.join('..', 'README.md')]},
    long_description=open('README.md').read(),
    include_package_data=True,
    platform=['Linux', 'Windows', 'Mac OS'],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Software Development',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    version='0.3.5',
    description='Visualisation of PASTML trees.',
    author='Anna Zhukova',
    author_email='anna.zhukova@pasteur.fr',
    url='https://gitlab.pasteur.fr/phylo/cytopast',
    download_url='https://gitlab.pasteur.fr/phylo/cytopast/archive/0.1.zip',
    keywords=['PASTML', 'visualisation', 'phylogeny'],
    install_requires=['ete3', 'pandas', 'numpy', 'jinja2', 'pastml'],
    requires=['ete3', 'pandas', 'numpy', 'jinja2', 'pastml']
)
