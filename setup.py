
from os import path
from setuptools import setup
import platon


# Get the long description from the README file
setup_dir = path.abspath(path.dirname(__file__))
with open(path.join(setup_dir, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


setup(
    name='cb-platon',
    version=platon.__version__,
    description='Platon: identification and characterization of bacterial plasmid contigs from short-read draft assemblies.',
    keywords=['bioinformatics', 'plasmids', 'wgs'],
    long_description=long_description,
    long_description_content_type='text/markdown',
    license='GPLv3',
    author='Oliver Schwengers',
    author_email='oliver.schwengers@computational.bio.uni-giessen.de',
    url='https://github.com/oschwengers/platon',
    packages=['platon'],
    python_requires='>=3.5',
    include_package_data=True,
    zip_safe=False,
    install_requires=[
        'biopython >= 1.71'
    ],
    entry_points={
        'console_scripts': [
            'platon=platon.platon:main'
        ]
    },
    classifiers=[
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3 :: Only',
        'Development Status :: 5 - Production/Stable',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Natural Language :: English'
    ],
    project_urls={
        'Bug Reports': 'https://github.com/oschwengers/platon/issues',
        'Source': 'https://github.com/oschwengers/platon'
    },
)
