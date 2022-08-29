
# Compilation script used to create executable file


# chmod u+x cmp.sh
# Uses Mersenne-Twister module for random number generation

 gfortran precision.f90 \
          mtmod.f90 \
          modulo_ti.f90 \
          main.f90 -l lapack -o main2.e -fcheck=all

 echo
