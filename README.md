## $$\text{Distributed parallel solution of a linear system using GMRES}$$
## $$\text{{Submitted to Professor Jack Dongarra}}$$
#### $$\text{Final Report: Scientific Computing for Engineers, CS594, Spring 2022}$$
### $$\text{By: Amirehsan Ghasmei*}$$
##### *Bredesen Center for Interdisciplinary Research and Graduate Education, University of Tennessee, Knoxville, TN 37996, USA. Email: aghasemi@vols.utk.edu
---
***Codes in the src folder:*** 
1)  **laplace_par.f90  => Explicit Solution** \
Suggest compile with something like:\
mpif90 -Og -g -fimplicit-none -fcheck=all -fbacktrace -pedantic -fbounds-check -Wall -Wextra -Wconversion -Wunderflow laplace_par.f90

2)  **custom_gmres.f90 => Implicit Solution**
---

Recommended Citation & Final report:\
DOI: https://doi.org/10.13140/RG.2.2.17501.42727


