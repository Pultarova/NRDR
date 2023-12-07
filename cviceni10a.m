
function abc
% CG a GMRES
clc;
load("Matice_MKP") % nacteni matice ze souboru (z programu 
% cond(A) % cislo podminenosti - kappa(A) - NAROCNY VYPOCET
N = size(A) % veliost matice A
A = sparse(A); % prevedeni na ridky format
tic; % zacatek mereni casu
u = A\B; % GEA
cas_GEA = toc % konec mereni casu
tic; % zacatek mereni casu
u = pcg(A,B,1e-6,500); % CG
cas_CG = toc % konec mereni casu

load("Matice_MS") % nacteni matice ze souboru
N = size(A) % veliost matice A
A = sparse(A); % prevedeni na ridky format
tic; % zacatek mereni casu
u = A\B; % GEA
cas_GEA = toc % konec mereni casu
tic; % zacatek mereni casu
u = gmres(A,B,20,1e-6,500); % GMRES
cas_GMRES = toc % konec mereni casu