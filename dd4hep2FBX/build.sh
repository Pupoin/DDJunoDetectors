source /home/wln/DD4hep/bin/thisdd4hep.sh
g++ dd4hep2FBXWriter.cc $(root-config --cflags --libs) -I$ROOT_INCLUDE_PATH -L$DD4HEP_LIBRARY_PATH -o run.exe