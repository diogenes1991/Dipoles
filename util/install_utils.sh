#!/bin/bash

BUILDDIR="${PWD}/built"
LOGFILE=install_utils.log

GSL=gsl-2.6
CUBA=Cuba-4.2
LHAPDF=LHAPDF-6.3.0

echo "--------------------------------------------" >> $LOGFILE
echo "------- Installing all dependancies --------" >> $LOGFILE
echo "--------------------------------------------" >> $LOGFILE

echo "" >> $LOGFILE
echo "" >> $LOGFILE
echo "" >> $LOGFILE

mkdir built

echo "--------------------------------------------" >> $LOGFILE
echo "------------ Now working on GSL ------------" >> $LOGFILE
echo "--------------------------------------------" >> $LOGFILE

echo "" >> $LOGFILE

echo "Extracting GSL"
tar xvzf "${GSL}.tar.gz" >> $LOGFILE
mv $LOGFILE $GSL
cd $GSL
echo "Configuring GSL"
./configure >> $LOGFILE
echo "Compiling GSL (this may take a while)"
make -j8 >> $LOGFILE 2>&1
mv $LOGFILE ../
cd ..

sed 's/####GSL PATH####/${PWD}/built/${GSL}/g' "../tpl/makefile"

LIBGSL=${GSL}/.libs/libgsl.a

if [ ! -f "$LIBGSL" ]
then
    echo ""
    echo "Error: ${LIBGSL} not found."
    echo "       Please check utilities installation."
    echo ""
    echo ""
    echo "" >> $LOGFILE
    echo "" >> $LOGFILE
    exit 1
else
    echo "GSL installation successful"
    mv $GSL built
fi

echo "" >> $LOGFILE
echo "" >> $LOGFILE
echo "" >> $LOGFILE

echo "--------------------------------------------" >> $LOGFILE
echo "----------- Now working on CUBA ------------" >> $LOGFILE
echo "--------------------------------------------" >> $LOGFILE

echo "" >> $LOGFILE

echo "Extracting CUBA"
tar xvzf ${CUBA}.tar.gz >> $LOGFILE
mv $LOGFILE $CUBA
cd $CUBA
echo "Configuring CUBA"
./configure >> $LOGFILE
echo "Compiling CUBA (this may take a while)"
make -j8 >> $LOGFILE 2>&1
mv $LOGFILE ../
cd ..

LIBCUBA=${CUBA}/libcuba.a

if [ ! -f "$LIBCUBA" ]
then
    echo ""
    echo "Error: ${LIBCUBA} not found."
    echo "       Please check utilities installation."
    echo ""
    echo ""
    echo "" >> $LOGFILE
    echo "" >> $LOGFILE
    exit 1
else
    echo "CUBA installation successful"
    mv $CUBA built
fi

echo "" >> $LOGFILE
echo "" >> $LOGFILE
echo "" >> $LOGFILE

echo "--------------------------------------------" >> $LOGFILE
echo "---------- Now working on LHAPDF -----------" >> $LOGFILE
echo "--------------------------------------------" >> $LOGFILE

echo "" >> $LOGFILE

echo "Extracting LHAPDF"
tar xvzf ${LHAPDF}.tar.gz >> $LOGFILE
mv $LOGFILE $LHAPDF
cd $LHAPDF
echo "Configuring LHAPDF"
./configure --prefix="${PWD}/built" >> $LOGFILE
echo "Compiling LHAPDF (this may take a while)"
make -j8 >> $LOGFILE 2>&1
make install >> $LOGFILE 2>&1
mv $LOGFILE ../
cd ..

LIBLHAPDF=${LHAPDF}/built/lib/libLHAPDF.a

if [ ! -f "$LIBLHAPDF" ]
then
    echo ""
    echo "Error: ${LIBLHAPDF} not found."
    echo "       Please check utilities installation."
    echo ""
    echo ""
    echo "" >> $LOGFILE
    echo "" >> $LOGFILE
    exit 1
else
    echo "LHAPDF installation successful"
    mv $LHAPDF built
fi

