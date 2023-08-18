gfortran -c kindset.f90
gfortran -c mod_splines.f90
gfortran -c tools.f90 -ffree-line-length-none
gfortran -o main main.f90 tools.o mod_splines.o kindset.o -llapack -I/opt/homebrew/include -I/usr/local/include -ffree-line-length-none
./main
rm *.mod *.o
