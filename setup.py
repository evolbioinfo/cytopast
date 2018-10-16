import os
from setuptools import setup

setup(
    name='cytopast',
    packages=['cytopast'],
    include_package_data=True,
    package_data={'cytopast': [os.path.join('templates', '*.html'), os.path.join('templates', '*.js'),
                               os.path.join('..', 'README.md')]},
    long_description=open('README.md').read(),
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Software Development',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    version='1.0.9',
    description='Visualisation of PastML trees.',
    author='Anna Zhukova',
    author_email='anna.zhukova@pasteur.fr',
    url='https://github.com/evolbioinfo/cytopast',
    keywords=['PASTML', 'visualisation', 'phylogeny'],
    install_requires=['ete3', 'pandas', 'numpy', 'jinja2', 'pastml>=1.0.8'],
    entry_points={
            'console_scripts': [
                'cytopast = cytopast.pastml_analyser:main',
            ]
    },
)
