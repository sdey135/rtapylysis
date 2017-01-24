from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='rtapylysis',
    version='0.3',
    description='Synchrotorn Rapid Thermal Annealing (RTA) X-ray Diffraction (XRD) Data Analysis Project',
    long_description=long_description,
    url='https://github.com/sdey135/rtapylysis',
    download_url = 'https://github.com/sdey135/rtapylysis/tarball/0.3',
    author='SD',
    author_email='sdey135@users.noreply.github.com',
    license='MIT',

    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Manufacturing',
        'Intended Audience :: Education',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: POSIX :: Linux',
        'Operating System :: Unix',
        'Operating System :: Microsoft :: Windows :: Windows 7',
        'Operating System :: MacOS',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Software Development :: Build Tools',
        'Topic :: Education',        
    ],

    keywords=[ 'Synchrotron', 'Rapid Thermal Annealing', 'RTA', 'X-ray Diffraction', 'XRD', 'Big Data' ],
    packages=['rtapylysis'],
    platforms='any',
    install_requires=['peppercorn'],

    package_data={
        'rtapylysis': ['example/*', 'dhkl_dir/*'],
    },

    entry_points={},
)
