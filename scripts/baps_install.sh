#!usr/bin/env sh
# To install BAPS v6 on Ubuntu 15.10
# Rupert A. Collins 12/03/16

# this puts BAPS in your Downloads folder, but it can be moved anywhere
cd ~/Downloads
wget http://www.helsinki.fi/bsg/software/BAPS/linux/BAPS6.0_Linux_64bit.zip
unzip BAPS6.0_Linux_64bit.zip
cd BAPS_package

# download the MCR
wget http://www.helsinki.fi/bsg/software/BEBaC/BEBaC_64bit/MCRInstaller_unix_2010a_64bit.bin
chmod 755 MCRInstaller_unix_2010a_64bit.bin

# install missing libs
wget http://fr.archive.ubuntu.com/ubuntu/pool/main/libx/libxp/libxp6_1.0.2-2_amd64.deb
# when software centre opens, click on 'install'
software-center libxp6_1.0.2-2_amd64.deb

# run the MCR auto installer
sudo ./MCRInstaller_unix_2010a_64bit.bin

# follow instructions on screen and install into 
'/opt/MATLAB/MATLAB_Compiler_Runtime'

# finally run BAPS, including path to MCR
./run_baps6.sh /opt/MATLAB/MATLAB_Compiler_Runtime/v713


