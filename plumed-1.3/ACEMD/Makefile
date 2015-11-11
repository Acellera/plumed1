CPPFLAGS=-I./ 
SRC=*.c

default: libplumed1plugin.so libplumed1plugin-rex.so

libplumed1plugin.so: ../common_files/*.c plug-in-metad.c
	gcc  $(CPPFLAGS) -O3 -DACEMD -fPIC --shared -o libplumed1plugin.so ../common_files/*.c plug-in-metad.c -I../common_files

libplumed1plugin-rex.so: ../common_files/*.c plug-in-metad-rex.c
	gcc  $(CPPFLAGS) -O3 -DACEMD -fPIC --shared -o libplumed1plugin-rex.so ../common_files/*.c plug-in-metad-rex.c -I../common_files

clean:
	rm -f lib*.so

