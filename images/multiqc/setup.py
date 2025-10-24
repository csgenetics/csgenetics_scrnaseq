import os
from setuptools import setup, find_packages

setup(
    name='multiqc-csgenetics',
    version='0.1.0',
    author='CS Genetics',
    author_email='support@csgenetics.com',
    description='Custom MultiQC plugins for CS Genetics scRNA-seq pipeline',
    long_description=open('README.md').read() if os.path.exists('README.md') else '',
    long_description_content_type='text/markdown',
    url='https://github.com/csgenetics/multiqc-plugins',
    packages=find_packages(),
    package_data={
        'unified_qc': ['*.yaml'],
    },
    entry_points={
        'multiqc.modules.v1': [
            'unified_qc = unified_qc:MultiqcModule',
        ],
        'multiqc.hooks.v1': [
            'before_config = unified_qc:register_plugin_search_patterns',
        ],
    },
    install_requires=[
        'multiqc>=1.14',
    ],
    python_requires='>=3.7',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)
