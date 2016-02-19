import setuptools
from setuptools.command.install import install


# Install Necessary R packages
class CustomInstallPackages(install):
    """Customized setuptools install command - runs R package install one-liner."""
    def run(self):
        import subprocess
        import shlex
        print "Attempting to install R packages...Please wait."
        cmd =''' R -e "install.packages(c('optparse', 'gtools', 'klaR','survival', 'mvtnorm', 'modeltools', 'coin', 'MASS'), repos = 'http://cran.stat.ucla.edu')" '''
        try:
            subprocess.call(shlex.split(cmd))
            print "Necessary R packages were sucessfully installed"
        except:
            print "Error installing R dependecies! Check to see if R is properly installed or see online documentation for more answers."
        install.run(self)



# Pkg info
setuptools.setup(
    name="koeken",
    version="0.2.4",
    url="https://github.com/twbattaglia/koeken",

    author="Thomas W. Battaglia",
    author_email="tb1280@nyu.edu",

    description="Linear Discriminant Analysis (LEfSe) on a Longitudinal Microbial Dataset.",
    long_description=open('README.rst').read(),
    keywords="Biology Microbiome LEFSE QIIME Formatting Diversity Python Bioinformatics",


    scripts=['koeken/koeken.py', 'koeken/lefse_src/format_input.py', 'koeken/lefse_src/run_lefse.py', 'koeken/lefse_src/lefse.py', 'koeken/lefse_src/plot_cladogram.py', 'koeken/lefse_src/export2graphlan.py', 'koeken/lefse_src/hclust2/hclust2.py', 'koeken/pretty_lefse.py'],

    cmdclass={'install': CustomInstallPackages},

    packages=setuptools.find_packages(),

    install_requires=['rpy2', 'numpy', 'matplotlib', 'argparse', 'pandas', 'biopython', 'qiime'],

    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4'
    ]
)
