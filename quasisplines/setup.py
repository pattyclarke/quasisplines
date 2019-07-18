from setuptools import setup

setup(
   name='quasi-splines',
   version='0.1.0',
   description='Symbolic Spline Computations in Sage',
   author='Patrick Clarke',
   author_email='pattyclarke@gmail.com',
   packages=['quasi-splines'],  #same as name
   install_requires=['sage'], #external packages as dependencies
)
