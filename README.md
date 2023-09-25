# Installation 
## QUEST
```
module purge
module load cmake/3.12.0 mpi/openmpi-4.1.1-gcc.10.2.0 gcc/9.2.0 hdf5/1.8.10 fftw/3.3.8-openmpi-4.0.5-gcc-10.2.0 kim-api/2.2.1-intel-19.0.5.281 blas-lapack/3.5.0_gcc hwloc/2.1.0-gcc-10.2.0 #glibc/2.25

cmake -D PKG_ASPHERE=on -D PKG_BODY=on -D PKG_CLASS2=on -D PKG_COLLOID=on -D PKG_COLVARS=on -D PKG_CORESHELL=-D PKG_DIPOLE=on -D PKG_EXTRA-COMPUTE=on -D PKG_EXTRA-DUMP=on -D PKG_EXTRA-FIX=on -D PKG_EXTRA-MOLECULE=on -D PKG_EXTRA-PAIR=on -D PKG_GRANULAR=on -D PKG_KIM=off -D PKG_KSPACE=on -D PKG_MANYBODY=on -D PKG_MC=on -D PKG_MEAM=on -D PKG_MISC=on -D PKG_ML-HDNNP=on -D PKG_ML-QUIP=on -D PKG_ML-SNAP=on -D PKG_MOLECULE=on -D PKG_MPIIO=on -D PKG_OPT=on -D PKG_PERI=on -D PKG_PHONON=on -D PKG_REAXFF=on -D PKG_REPLICA=on -D PKG_RIGID=on -D PKG_SHOCK=on -D PKG_SRD=on -D PKG_USER-MLIP=on -D PKG_USER-VCSGC=onG_ASPHERE=on -D PKG_BODY=on -D PKG_CLASS2=on -D PKG_COLLOID=on -D PKG_COLVARS=on -D PKG_CORESHELL=-D PKG_DIPOLE=on -D PKG_EXTRA-COMPUTE=on -D PKG_EXTRA-DUMP=on -D PKG_EXTRA-FIX=on -D PKG_EXTRA-MOLECULE=on -D PKG_EXTRA-PAIR=on -D PKG_GRANULAR=on -D PKG_KIM=on -D PKG_KSPACE=on -D PKG_MANYBODY=on -D PKG_MC=on -D PKG_MEAM=on -D PKG_MISC=on -D PKG_ML-HDNNP=on -D PKG_ML-QUIP=on -D PKG_ML-SNAP=on -D PKG_MOLECULE=on -D PKG_MPIIO=on -D PKG_OPT=on -D PKG_PERI=on -D PKG_PHONON=on -D PKG_REAXFF=on -D PKG_REPLICA=on -D PKG_RIGID=on -D PKG_SHOCK=on -D PKG_SRD=on -D PKG_USER-MLIP=on -D PKG_USER-VCSGC=on -D BUILD_MPI=yes DOWNLOAD_KIM=off ../cmake
make -j29
```

## RadonPy
```
conda create -n <envname> python==3.10
conda activate <envname>
conda install -c psi4 psi4 resp
conda install -c conda-forge rdkit matplotlib mdtraj
pip install radonpy-pypi
```

