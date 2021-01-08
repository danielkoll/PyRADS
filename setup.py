from setuptools import setup

dependencies = ["numpy",
                "scipy"]

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='PyRADS',
      version='0.1.0',
      description='PyRADS is the "Python line-by-line RADiation model for planetary atmosphereS"',
      long_description=readme(),
      url='',
      author='Daniel B. Koll',
      author_email='dbkoll@mit.edu',
      license='MIT',
      packages=['pyrads'],
      install_requires=dependencies,
      zip_safe=False)
