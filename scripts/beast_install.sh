#!/usr/bin/env sh

# installing BEAST
cd Downloads #change to whatever
wget https://github.com/beast-dev/beast-mcmc/releases/download/v1.8.3/BEASTv1.8.3.tgz
tar -zxvf BEASTv1.8.3.tgz
export PATH=~/Downloads/BEASTv1.8.3/bin:$PATH
echo 'export PATH=~/Downloads/BEASTv1.8.3/bin:$PATH' >> ~/.bashrc

# installing BEAGLE lib without CUDA
# installs to '/usr/local/bin'
git clone https://github.com/beagle-dev/beagle-lib.git
cd beagle-lib
git tag
git checkout beagle_release_2_1_2
./autogen.sh
./configure
make
make check
sudo make install
make clean
make distclean
echo 'export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH' >> ~/.bashrc
git checkout master

# test it's working
beast -h

# run BEAST on an UNPARTITIONED dataset from ANY directory ...
beast -threads 4 -beagle_instances 4 -beagle -beagle_SSE -beagle_CPU -beagle_double -beagle_scaling none -seed 300316 -overwrite ddRAD.gblocks.xml
