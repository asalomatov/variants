from setuptools import setup, find_packages

if __name__ == '__main__':
    setup(
        name='variants',
        packages=find_packages(),
        version='0.1.0',
        description='Genomic variants manipulation and filtering.',
        long_description='',
        url='https://github.com/asalomatov/variants',
        author='Andrei Salomatov',
        license='MIT',
        classifiers=[
            'Development Status :: 3 - Alpha',
            'Environment :: Console',
            'Operating System :: OS Independent',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: MIT License', 
            'Programming Language :: Python',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            ],
        install_requires=[
            'pandas>=0.17.1',
            'numpy>=1.10.2',
            'keras>=0.3.3'
        ],
)
