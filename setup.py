"""
Koeken setup

To run: python setup.py install

"""
try:
    import setuptools
except ImportError:
    sys.exit("Please install setuptools.")

# Pkg info
setuptools.setup(
    name="koeken",
    version="0.3.0",
    url="https://github.com/twbattaglia/koeken",
    author="Thomas W. Battaglia",
    author_email="tb1280@nyu.edu",
    description="A Linear Discriminant Analysis effect size (LEfSe) wrapper.",
    keywords=['microbial', 'microbiome', 'bioinformatics', 'LEfSe',
              'metagenomic', 'QIIME', 'koeken', 'biology'],
    platforms=['Linux', 'MacOS'],
    classifiers=[
        "Programming Language :: Python",
        'Development Status :: 2 - Pre-Alpha',
        "Environment :: Console",
        "Operating System :: MacOS",
        "Operating System :: Unix",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.4",
        "Topic :: Scientific/Engineering :: Bio-Informatics"],
    entry_points={
        'console_scripts': [
            'koeken = koeken.koeken:main',
            'pretty_lefse = humann2.tools.humann2_databases:main']},
    packages=setuptools.find_packages(),
    install_requires=['rpy2', 'argparse', 'pandas', 'biopython']
)
