#! /bin/sh 

# Gemini GMOS cookbook
# Reduction of Long-Slit Spectra with PyRAF
# https://noirlab.edu/science/programs/csdc/usngo/gmos-cookbook/Processing/PyrafProcLS.html#retrieve-organize-data

# create directory
if [ ! -e raw ]
  then
    mkdir raw 
fi
if [ ! -e work ]
  then
    mkdir work 
fi

# copy tar file into ./raw dir
cp raw_tar/*.tar raw/

cd raw
rm *.fits

# extract the tarball
cat *.tar | tar -xvf - -i
rm *.tar
bunzip2 *.bz2 
rm *.bz2 *.txt

# Remove the exisitng combine bias file

rm *bias.fits

# download and build observing log database
if [ ! -e obslog.py ]
  then
  #wget http://ast.noao.edu/sites/default/files/GMOS_Cookbook/_downloads/obslog.py
  wget https://noirlab.edu/science/programs/csdc/usngo/gmos-cookbook/_downloads/a305aea0feb0194376c731fd9cb6c783/obslog.py
fi
rm obsLog.sqlite3
python obslog.py obsLog.sqlite3

# download the python file selection module
cd ../work
if [ ! -e fileSelect.py ]
  then
  #wget http://ast.noao.edu/sites/default/files/GMOS_Cookbook/_downloads/fileSelect.py
  wget https://noirlab.edu/science/programs/csdc/usngo/gmos-cookbook/_downloads/24dc1fedce8c8b74d809e47526826207/fileSelect.py
fi



