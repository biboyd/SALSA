from setuptools import setup, find_packages
setup(name="salsa",
      version="0.0.0",
      author="Brendan Boyd",
      packages=find_packages(),
      install_requires=[
        'numpy',
        'yt',
        'trident',
        'spectacle',
        'matplotlib',
		'pandas'
      ])
