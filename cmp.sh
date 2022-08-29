
# chmod u+x cmp.sh
# Si se usa el generador de Marsaglia reemplzar mtmod.f90 por mzran.f90


 gfortran precision.f90 \
          mtmod.f90 \
          modulo_ti.f90 \
          main.f90 -l lapack -o main2.e -fcheck=all

 echo
