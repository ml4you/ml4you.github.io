#!/bin/sh
echo "
*************************************************
*                                               *
*           GROMACS 2019.4 INSTALLER            *
*                  BY YASSIR B.                 *
*                                               *
*************************************************
                ";
while true; do
    read -p "Do you wish to install this program? " yn
    case $yn in
        [Yy]* ) make install; break;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac
done
read yesno
pwd
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install cmake
cmake --version
sudo apt-get install build-essential
mkdir gromacs/
cd gromacs/
wget http://ftp.gromacs.org/pub/gromacs/gromacs-2019.4.tar.gz
tar xvzf gromacs-2019.4.tar.gz
sudo apt-get install libfftw3-dev
wget http://gerrit.gromacs.org/download/regressiontests-2019.4.tar.gz
tar xvzf regressiontests-2019.4.tar.gz
cd gromacs-2019.4/
mkdir build
cd build
sudo cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=OFF -DCMAKE_C_COMPILER=gcc -DREGRESSIONTEST_PATH=./../../regressiontests-2019.4
sudo make check
sudo make install
source /usr/local/gromacs/bin/GMXRC
echo "Done installing"
gmx pdb2gmx --version
