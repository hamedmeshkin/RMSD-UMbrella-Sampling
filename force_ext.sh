# on big red2 before compiling the code, gcc should be load first.
g++ -shared -o Rms_ext.so  -DUSE_TCL_STUBS -DUSE_TK_STUBS -I$TCLINC -L$TCLLIB -ltclstub8.5 -fPIC frc_ext.c
g++ -shared -o RmsF_ext.so -DUSE_TCL_STUBS -DUSE_TK_STUBS -I$TCLINC -L$TCLLIB -ltclstub8.5 -fPIC frc_ext.c