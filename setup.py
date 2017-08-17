from setuptools import setup, find_packages
from glob import glob

if __name__ == '__main__':
    setup(
        name='variants',
        packages=['variants'],
        version='0.1.0',
        package_data={'variants': ['denovo_classifier_model_SNP/*',
                                   'denovo_classifier_model_INDEL/*',
                                   'denovo_classifier_config/*', '*.sh']},
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
            'scikit-learn==0.18.1',
            'pandas==0.19.2',
            'numpy==1.11.3',
            'keras==1.1.1',
            'seaborn==0.6.0',
            'xgboost==0.6a2',
            'matplotlib==2.0.0',
            'pysam==0.8.3',
            'PyVcf'
        ],
    )
