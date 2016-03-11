if [ "$CC" == "" ]; then CC=gcc; fi
if [ "$FC" == "" ]; then FC=gfortran; fi
if [ "$CXX" == "" ]; then CXX=g++; fi


DIR="$CONDA_DEFAULT_ENV"
mkdir -p "$DIR/lib/"
mkdir -p "$DIR/bin/"

printenv

cd ACEMD
make TCL="-I$SYS_PREFIX/include -L$SYS_PREFIX/lib" CC=$CC
cp libplumed1plugin.so libplumed1plugin-rex.so "$DIR/lib"

cd ../utilities/sum_hills
make FC=$FC
cp sum_hills "$DIR/bin"





