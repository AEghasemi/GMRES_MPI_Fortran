# Distributed parallel solution of a linear system using GMRES
## Submitted to Professor Jack Dongarra
### Final report: Scientific Computing for Engineers, CS594, Spring 2022
#### By: Amirehsan Ghasmei
##### Bredesen Center for Interdisciplinary Research and Graduate Education, University of Tennessee, Knoxville, TN 37996, USA. Email: aghasemi@vols.utk.edu
---
Codes in the src folder: 
1)  **laplace_par.f90  => Explicit Solution** \
Suggest compile with something like:\
mpif90 -Og -g -fimplicit-none -fcheck=all -fbacktrace -pedantic -fbounds-check -Wall -Wextra -Wconversion -Wunderflow laplace_par.f90

2)  **custom_gmres.f90 => Implicit Solution**
---






<video src="https://github.com/AEghasemi/GMRES_MPI_Fortran/assets/93692630/7cb6a772-0233-44da-a6e0-9703a252d2d2"/>


