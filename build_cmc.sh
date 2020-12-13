module purge all
module load python/anaconda
module load cuda
module unload gcc
module unload mpi
module load gcc/4.8.3
module load mpi/openmpi-1.10.5-gcc-4.8.3
module load gsl/1.16-gcc4.8.3
module load CfitsIO
module load intel

autoreconf -ivf

F77=mpifort FC=mpifort CXX=mpic++ CC=mpicc ./configure --prefix=/projects/b1095/syr904/cmc/cmc-mpi-tidalcapture/ --enable-mpi --enable-experimental #--enable-tausworth


#module purge all
#module load CfitsIO
#module unload gcc/4.8.3
#module load gsl/2.5-gcc-6.4.0
#module load mpi/openmpi-3.1.3-gcc-6.4.0
