from setuptools import setup
setup(name="CGM",
      version="0.0.0",
      author="Brendan Boyd",
      packages=["CGM"],
      install_requires=[
        'numpy',
        'yt',
        'trident',
        'spectacle',
        'matplotlib'
      ])
