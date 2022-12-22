from setuptools import setup

long_description = open('README.md', encoding='utf-8').read()

setup(name='FastORF',
      version='0.0',
      description='FastORF: A fast and inaccurate open reading frame finder',
      long_description = long_description,
      long_description_content_type = 'text/markdown',
      classifiers=[
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Artificial Intelligence',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
      ],
      url='https://github.com/luispedro/fastorf',
      author='Luis Pedro Coelho',
      author_email='luispedro@big-data-biology.org',
      license='MIT',
      packages = ['fastorf'],
      install_requires=[],
      package_data={},
      zip_safe=False,
      entry_points={}
)

