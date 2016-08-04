#!/usr/bin/env sh
# script to install fastStructure
# works on Ubuntu Linux 16.04 (July 2016)


# install SciPy, Numpy and Cython

sudo apt-get install python-scipy
sudo apt-get install python-numpy
sudo apt-get install Cython

# build dependencies of matplotlib, and install matplotlib

sudo apt-get build-dep python-matplotlib
sudo apt-get install python-matplotlib

# put matplotlib libraries into the compile path
# open the ".bashrc" file in the home directory and insert at the end of the file these directives
# the "export LD_LIBRARY_PATH" should be there by default already
# the ".bashrc" file is hiden and needs to be visualized in GUI via CRTL+h

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
export CFLAGS="-I/usr/local/include"
export LDFLAGS="-L/usr/local/lib"

# get fastStructure via git (do this from a directory where you keep programs - e.g. ~/bin)

git clone https://github.com/rajanil/fastStructure

# change into fastStructure subdirectory var, and build (assumes fastStructure is in ~/bin directory)

cd ~/bin/fastStructure/vars/
python setup.py build_ext --inplace

# change to fastStructure directory, and compile

cd ~/bin/fastStructure/
python setup.py build_ext --inplace

# test compilation by running

python ~/bin/fastStructure/structure.py