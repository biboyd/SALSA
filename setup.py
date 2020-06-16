from setuptools import setup
setup(name="salsa",
      version="0.0.0",
      author="Brendan Boyd",
      packages=["salsa"],
      install_requires=[
        'numpy',
        'yt',
        'trident',
        'spectacle',
        'matplotlib',
		'pandas'
      ])
