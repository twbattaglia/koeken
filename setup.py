import setuptools

setuptools.setup(
    name = "koeken",
    version = "0.2.1",
    url = "[https://github.com/twbattaglia/koeken",

    author = "Thomas W. Battaglia",
    author_email = "tb1280@nyu.edu",

    description = "Linear Discriminant Analysis (LEfSe) on A Longitudinal Microbial Dataset.",
    long_description = open('README.rst').read(),
    keywords = "Biology Microbiome LEFSE QIIME Formatting Diversity Python Bioinformatics",


    scripts = ['koeken/koeken.py','koeken/lefse_src/format_input.py', 'koeken/lefse_src/run_lefse.py', 'koeken/lefse_src/lefse.py', 'koeken/lefse_src/plot_cladogram.py', 'koeken/lefse_src/export2graphlan.py', 'koeken/lefse_src/hclust2/hclust2.py'],

    #entry_points={
        #'console_scripts':
            #['koeken = koeken.koeken:main']
        #},


    packages=setuptools.find_packages(),

    install_requires = ['rpy2', 'numpy', 'matplotlib', 'argparse', 'pandas', 'biopython', 'qiime'],

    classifiers = [
        'Development Status :: 2 - Pre-Alpha',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
    ],
)
