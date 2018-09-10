import setuptools


def readme():
    with open('README.rst') as f:
        return f.read()

setuptools.setup(name='signature-invariants-py',
      #python_requires='>=3.5.2',
      version='0.1.1',
      packages=['signature_invariants'],
      description='Implements certain invariants for the iterated-integrals signature.',
      long_description=readme(),
      author='Joscha Diehl',
      #author_email='',
      url='https://github.com/diehlj/signature-invariants-py',
      license='Eclipse Public License',
      install_requires=['numpy', 'scipy', 'sympy', 'linear-combination-py'],
      #setup_requires=['setuptools_git >= 0.3', ],
      test_suite='tests'
      )
