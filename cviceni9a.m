% MKP pro 2d ulohu
% -div.(a grad u)+pu=f na (0,L1)x(0,L2), 
% Na hranici je Dirichletova OP (u=gD) nebo homog.Neumann (a*grad(u)*n=0) podle
% podminky NeumannOP=0 nebo 1.
% Oblast je L-shape (L_shape=1 nebo 0).

function abc
clc;
L1 = 2; % delka oblasti ve smeru x1
L2 = 1; % delka oblasti ve smeru x2
N1 = 18; % SUDY pocet intervalu ve smeru x1
N2 = 12; % SUDY pocet intervalu ve smeru x2
x1 = linspace(0,L1,N1+1); % sit bodu ve smeru x1
x2 = linspace(0,L2,N2+1); % sit bodu ve smeru x2
NeumannOP = 1; % 0 nebo 1 (ano nebo ne pro homog. Neumann. OP pro x1=0)
L_shape = 1;
Q = 0.8; % jak velky gradient vzhledem k nejvetsimu je pripustny

[xx1,xx2] = meshgrid(x1,x2); % vsechny uzly
xx1 = xx1(:); % x1 souradnice srovnane do vektoru
xx2 = xx2(:); % x2 souradnice srovnane do vektoru

Nuzel = size(xx1,1); % pocet uzlu

% uprava oblasti do L-shape:
if (L_shape==1)
    p = Nuzel;
    for k = p:-1:1 % ponechame jen uzly v L-shape oblasti
        if (xx1(k)<L1/2 && xx2(k)<L2/2) xx1(k)=[]; xx2(k)=[]; Nuzel = Nuzel-1; end;    
    end;
end;

for kk = 1:3 % hlavni cyklus s adaptivni siti %<<<<<<<<<<<<<<<<<<<<

elem = delaunay(xx1,xx2); % sada elementu (trojuhelniku) pro dane uzly
Nelem = size(elem,1); % pocet elementu
% KresliElem(elem,xx1,xx2,Nelem);

% uprava oblasti do L-shape - NENI TRIVIALNI - DISKUSE !
if (L_shape==1)
    p = Nelem;
    for k = Nelem:-1:1 % vymaze elementy mimo L-shape
        pom = elem(k,:);
        x1s = sum(xx1(pom))/3;
        x2s = sum(xx2(pom))/3;
        if (x1s<L1/2 && x2s<L2/2) elem(k,:)=[]; Nelem=Nelem-1; end;
    end;
end;

% oznaceni hranicnich uzlu:
if (L_shape==1)
    uzel_hra = Najdi_uzel_hra_L(xx1,xx2,Nuzel,L1,L2,NeumannOP); % najdeme hranicni uzly
else
    uzel_hra = Najdi_uzel_hra(xx1,xx2,Nuzel,L1,L2); % najdeme hranicni uzly
end;

% sestaveni matice A a prave strany B:
A = zeros(Nuzel);
B = zeros(Nuzel,1);
% jdu po prvcich:
for k = 1:Nelem
    cu = elem(k,:); % cisla uzlu
    A(cu,cu) = A(cu,cu)+Integruj_a(k,elem,xx1,xx2)+Integruj_p(k,elem,xx1,xx2);
    B(cu) = B(cu)+Integruj_f(k,elem,xx1,xx2);
end;
upom = zeros(Nuzel,1);
for k = Nuzel:-1:1 % uprava matice A a prave strany B
    if (uzel_hra(k)==1) 
        B = B-A(:,k)*FcegD(xx1(k),xx2(k));
        A(k,:)=[]; A(:,k)=[]; B(k)=[];
        upom(k) = FcegD(xx1(k),xx2(k)); % dosazeni Dirichlet. OP
    end;
end;
u = A\B; % reseni soustavy lin. rovnic

zde = 1;
for k = 1:Nuzel
    if (uzel_hra(k)==0)
        upom(k) = u(zde); % doplneni reseni ve vnitrnich uzlech
        zde = zde+1;
    end;
end;

% adaptivita site:
[elem,Nelem,xx1,xx2,Nuzel] = Zmena_site(elem,Nelem,upom,xx1,xx2,Nuzel,uzel_hra,Q);
Nuzel % vypise se, kolik ted mame uzlu
end;

KresliReseni(elem,xx1,xx2,Nelem,upom); % vykresleni posledniho reseni

% save('MaticeMKP1',"A","B");


%====================================================================

function z = KresliElem(elem,xx1,xx2,Nelem)
cla; hold on;
for k = 1:Nelem
    t = elem(k,:);
    plot(xx1([t,t(1)]),xx2([t,t(1)]));
end;

function z = Najdi_uzel_hra(xx1,xx2,Nuzel,L1,L2)
epsi = 1e-6;
uzel_hra = zeros(Nuzel,1);
for k = 1:Nuzel
    if (abs(xx1(k)-L1)<epsi || abs(xx1(k))<epsi || abs(xx2(k)-L2)<epsi || abs(xx2(k))<epsi)
        uzel_hra(k) = 1;
    end;
end;
z = uzel_hra;

function z = Najdi_uzel_hra_L(xx1,xx2,Nuzel,L1,L2,Neumann)
epsi = 1e-8;
if (Neumann>0) cc = 1; else cc = 0; end;
uzel_hra = zeros(Nuzel,1);
for k = 1:Nuzel
    if (abs(xx1(k)-L1)<epsi || abs(xx1(k)-L1*cc)<epsi || ...
            abs(xx2(k)-L2)<epsi || abs(xx2(k))<epsi)
        uzel_hra(k) = 1; continue;
    end;
    if ((abs(xx1(k)-L1/2)<epsi && xx2(k)<=L2/2) || ...
            (abs(xx2(k)-L2/2)<epsi && xx1(k)<=L1/2))
        uzel_hra(k) = 1; continue;
    end;
end;
z = uzel_hra;

function z = KresliReseni(elem,xx1,xx2,Nelem,u)
cla; hold on;
for k = 1:Nelem
    t = elem(k,:);
    surf([xx1(t(1)),xx1(t(2));xx1(t(3)),xx1(t(3))],...
        [xx2(t(1)),xx2(t(2));xx2(t(3)),xx2(t(3))],...
        [u(t(1)),u(t(2));u(t(3)),u(t(3))],'FaceAlpha',0.5);
end;
view(10,15);

function z = Integruj_a(k,elem,xx1,xx2)
t = elem(k,:);
x1c = sum(xx1(t)/3);
x2c = sum(xx2(t)/3);
pom  = [xx1(t),xx2(t),[1;1;1]];
d = [pom(1,1:2)-pom(3,1:2);pom(2,1:2)-pom(3,1:2)];
d = abs(det(d))/2;
pom = inv(pom);
pom = pom(1:2,:);
z = pom'*Fcea(x1c,x2c)*pom*d;

function z = Integruj_p(k,elem,xx1,xx2)
t = elem(k,:);
x1c = sum(xx1(t)/3);
x2c = sum(xx2(t)/3);
pom  = [xx1(t),xx2(t)];
d = [pom(1,:)-pom(3,:);pom(2,:)-pom(3,:)];
d = abs(det(d))/2;
z = Fcep(x1c,x2c)*[2,1,1;1,2,1;1,1,2]*d/12;

function z = Integruj_f(k,elem,xx1,xx2)
t = elem(k,:);
x1c = sum(xx1(t)/3);
x2c = sum(xx2(t)/3);
pom  = [xx1(t),xx2(t)];
d = [pom(1,:)-pom(3,:);pom(2,:)-pom(3,:)];
d = abs(det(d))/2;
z = Fcef(x1c,x2c)*d/3;

function z = Fcea(x1,x2)
z = eye(2)/5;%+x1^2;
function z = Fcep(x1,x2)
z = 0.1;
function z = Fcef(x1,x2)
z = 1;%-x1-x2;
function z = FcegD(x1,x2)
z = 1+x2/30+x1/30;
z = 0;
function z = FcegN(x1,x2)
z = 1+x1/10;

function [elem,Nelem,xx1,xx2,Nuzel] = Zmena_site(elem,Nelem,upom,xx1,xx2,Nuzel,uzel_hra,q)
max_deriv = 0;
for k = 1:Nelem
    p = elem(k,:);
    x1p = xx1(p);
    x2p = xx2(p);
    up = upom(p);
    a = [x1p,x2p,ones(3,1)];
    pom = inv(a)*up;
    mm = norm([pom(1),pom(2)]);
    if (mm>max_deriv) max_deriv=mm; end;
end;
for k = 1:Nelem
    p = elem(k,:);
    x1p = xx1(p);
    x2p = xx2(p);
    up = upom(p);
    a = [x1p,x2p,ones(3,1)];
    pom = inv(a)*up;
    mm = norm([pom(1),pom(2)]);
    if (mm>max_deriv*2/4) 
        x1 = [sum(x1p([1,2]))]/2;%,sum(x1p([1,3])),sum(x1p([2,3]))]/2;
        x2 = [sum(x2p([1,2]))]/2;%,sum(x2p([1,3])),sum(x2p([2,3]))]/2;
        x1s = sum(x1p)/3;
        x2s = sum(x2p)/3;
        % xx1 = [xx1;x1s];
        % xx2 = [xx2;x2s];
        xx1 = [xx1;x1s];
        xx2 = [xx2;x2s];
        Nuzel = Nuzel+1;
        if (sum(uzel_hra(p([1,2])))==2) 
            x1 = [sum(x1p([1,2]))]/2;
            x2 = [sum(x2p([1,2]))]/2;
            xx1 = [xx1;x1];
            xx2 = [xx2;x2];
            Nuzel = Nuzel+1;
        end;
        if (sum(uzel_hra(p([1,3])))==2) 
            x1 = [sum(x1p([1,3]))]/2;
            x2 = [sum(x2p([1,3]))]/2;
            xx1 = [xx1;x1];
            xx2 = [xx2;x2];
            Nuzel = Nuzel+1;
        end;
        if (sum(uzel_hra(p([2,3])))==2) 
            x1 = [sum(x1p([2,3]))]/2;
            x2 = [sum(x2p([2,3]))]/2;
            xx1 = [xx1;x1];
            xx2 = [xx2;x2];
            Nuzel = Nuzel+1;
        end;
        
    end;
end;




