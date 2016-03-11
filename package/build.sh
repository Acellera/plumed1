CC=gcc-4.4
FC=gfortran-4.4

which $CC
if [ "$?" != "0" ]; then
  CC=gcc
fi
which $FC
if [ "$?" != "0" ]; then
  FC=gfortran
fi



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





