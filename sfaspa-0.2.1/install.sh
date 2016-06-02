#/bin/bash

# Install perl libs
echo -e "\nInstalling perl libraries ..."
./perl_libs.sh


export SPA_HOME=`pwd`

echo -e "\nInstalling FragGeneScan ..."
cd $SPA_HOME/3rd
if [ -d FragGeneScan1.16 ]; then 
	rm -rf FragGeneScan1.16
fi
tar zxvf FragGeneScan1.16.tar.gz
cd FragGeneScan1.16
make clean
make fgs
[ $? -ne 0 ] && exit $?

echo -e "\nInstalling MetaGeneAnnotator ..."
cd $SPA_HOME/3rd
if [ -d mga_ia64 ]; then 
	rm -rf mga_ia64
fi
tar zxvf mga_ia64.tar.gz
ln -s $SPA_HOME/3rd/mga_ia64/mga_linux_ia64 $SPA_HOME/3rd/mga_ia64/mga_linux

echo -e "\nInstalling libdivsufsort ..."
cd $SPA_HOME/3rd
if [ -d libdivsufsort-2.0.1 ]; then 
	rm -rf libdivsufsort-2.0.1
fi
tar xvf libdivsufsort-2.0.1.tar
cd libdivsufsort-2.0.1
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE="Release" -DCMAKE_INSTALL_PREFIX=`pwd` ..
[ $? -ne 0 ] && exit $?
make
make install
[ $? -ne 0 ] && exit $?

## SPA installation
echo -e "\nInstalling SPA ..."
if [ ! -e $SPA_HOME/bin ]; then 
	mkdir $SPA_HOME/bin 
fi
cd $SPA_HOME/src
make all
[ $? -ne 0 ] && exit $?

cd $SPA_HOME

echo -e "\nDone"
