Prerequisities
==============
You need
* CMake
* compiler (gcc / visual studio)
* BLAS, LAPACK
* UMFPACK (if you want to use the code without MPI)
* P4EST (if you want to use the code with MPI)
    * installs typically with ./p4est-setup.sh p4est-1.1.tar.gz /.../p4est-build
* OpenMPI (if you want to use the code with MPI)

Trilinos
========
Get the commit

  **b260ac0579d90dbb07820525c5ce393a7be3e193**

Use this CMake command (which you can fine-tune, e.g. add TBB support) if you **want to** use MPI::
```
mkdir build
cd build
cmake -D TPL_ENABLE_Pthread=OFF -D Trilinos_ENABLE_TESTS=OFF -DCMAKE_BUILD_TYPE:STRING=RELEASE -DTrilinos_ENABLE_Amesos=ON -DTrilinos_ENABLE_Epetra=ON -DTrilinos_ENABLE_EpetraExt=ON -DTrilinos_ENABLE_Ifpack=ON -DTrilinos_ENABLE_AztecOO=ON -DTrilinos_ENABLE_Sacado=ON -DTrilinos_ENABLE_Teuchos=ON -DTrilinos_ENABLE_MueLu=ON -DTrilinos_ENABLE_ML=ON -DTrilinos_VERBOSE_CONFIGURE=OFF -DTPL_ENABLE_MPI=ON -DBUILD_SHARED_LIBS=OFF -DCMAKE_VERBOSE_MAKEFILE=OFF -D Trilinos_ENABLE_ALL_PACKAGES=OFF -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF -DTrilinos_ENABLE_Gtest=OFF -DTrilinos_ENABLE_Kokkos=OFF -D MPI_BASE_DIR=/usr -D TPL_ENABLE_TBB:BOOL=OFF -DCMAKE_CXX_FLAGS:STRING=-fPIC -DCMAKE_LINKER_FLAGS=-fPIC ..
```
Use this CMake command if you do **NOT want to** use MPI::

```
mkdir build
cd build
cmake -D TPL_ENABLE_Pthread=OFF -D Trilinos_ENABLE_TESTS=OFF -DCMAKE_BUILD_TYPE:STRING=RELEASE -DTrilinos_ENABLE_Amesos=ON -DTrilinos_ENABLE_Epetra=ON -DTrilinos_ENABLE_EpetraExt=ON -DTrilinos_ENABLE_Ifpack=ON -DTrilinos_ENABLE_AztecOO=ON -DTrilinos_ENABLE_Sacado=ON -DTrilinos_ENABLE_Teuchos=ON -DTrilinos_ENABLE_MueLu=ON -DTrilinos_ENABLE_ML=ON -DTrilinos_VERBOSE_CONFIGURE=OFF -DTPL_ENABLE_MPI=OFF -DBUILD_SHARED_LIBS=OFF -DCMAKE_VERBOSE_MAKEFILE=OFF -D Trilinos_ENABLE_ALL_PACKAGES=OFF -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF -DTrilinos_ENABLE_Gtest=OFF -DTrilinos_ENABLE_Kokkos=OFF -D TPL_ENABLE_TBB:BOOL=OFF -DCMAKE_CXX_FLAGS:STRING=-fPIC -DCMAKE_LINKER_FLAGS=-fPIC ..
```

Deal.II
=======
Get the commit::

  **d7c645108ddca62905c510f5d33a40c0c87f952e**
  
Use this CMake command if you **want to** use MPI

```
  cmake -DCMAKE_BUILD_TYPE=Release -DDEAL_II_WITH_MPI:BOOL=ON -DDEAL_II_CXX_FLAGS:STRING=-fPIC -DDEAL_II_LINKER_FLAGS:STRING=-fPIC -DDEAL_II_WITH_P4EST:BOOL=ON -DDEAL_II_WITH_TRILINOS:BOOL=ON .
```  

Use this CMake command if you do **NOT want** to use MPI::
```
  cmake -DCMAKE_BUILD_TYPE=Release -DDEAL_II_WITH_MPI:BOOL=OFF -DDEAL_II_CXX_FLAGS:STRING=-fPIC -DDEAL_II_LINKER_FLAGS:STRING=-fPIC -DDEAL_II_WITH_TRILINOS:BOOL=ON .
```

MHDeal
======
Use this CMake command:
```
cmake -DDEAL_II_DIR=... ,
```

using the root directory of your deal.II repository or installation directory.
