#!/usr/bin/env sh

# website
# https://bitbucket.org/gutenkunstlab/dadi/overview

# install dependencies
sudo apt-get install python-numpy python-scipy python-matplotlib ipython ipython-notebook python-pandas python-sympy python-nose

# install dadi
cd ~/Software
git clone https://bitbucket.org/gutenkunstlab/dadi.git
cd dadi
sudo python setup.py install

# check it installed correctly
cd tests
python run_tests.py

# ALL GOOD :)