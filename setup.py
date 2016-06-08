from setuptools import setup, find_packages
from glob import glob

if __name__ == '__main__':
    setup(
        name='variants',
        packages=['variants'],
        version='0.1.0',
        package_data={'variants': ['denovo_classifier_model_SNP/*', 'denovo_classifier_config/*']},
        include_package_data=True,
        scripts=['variants/scripts/call_de_novo.py'],
        description='Genomic variants manipulation and filtering',
        long_description='',
        url='https://github.com/asalomatov/variants',
        author='Andrei Salomatov',
        author_email='Andrei.Salomatov@gmail.com',
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
            'sklearn',
            'pandas>=0.17.1',
            'numpy>=1.10.2',
            'keras>=0.3.3'
        ],
)
