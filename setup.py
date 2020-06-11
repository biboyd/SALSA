from setuptools import setup
setup(name="SALS",
      version="0.0.0",
      author="Brendan Boyd",
      packages=["SALS"],
      install_requires=[
        'numpy',
        'yt',
        'trident',
        'spectacle',
        'matplotlib',
		'pandas'
      ])
