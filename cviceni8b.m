% MKP pro 2d ulohu
% -div.(a grad u)+pu=f na (0,L1)x(0,L2), 
% na hranici Dirichlet u=gD (nebo Neumann a*grad(u)*n=gN)

% v tomto programu je pouze Dirichletova OP

function abc
clc;
L1 = 1; % delka oblasti ve smeru x1
L2 = 1; % delka oblasti ve smeru x2
N1 = 18; % pocet intervalu ve smeru x1
N2 = 17; % pocet intervalu ve smeru x2
x1 = linspace(0,L1,N1+1); % sit bodu ve smeru x1
x2 = linspace(0,L2,N2+1); % sit bodu ve smeru x2
[xx1,xx2] = meshgrid(x1,x2); % vsechny uzly
xx1 = xx1(:); % x1 souradnice srovnane do vektoru
xx2 = xx2(:); % x2 souradnice srovnane do vektoru
elem =delaunay(xx1,xx2); % sada elementu (trojuhelniku)
Nelem = size(elem,1); % pocet elementu
Nuzel = size(xx1,1); % pocet uzlu
% KresliElem(elem,xx1,xx2,Nelem);
uzel_hra = Najdi_uzel_hra(xx1,xx2,Nuzel,L1,L2); % najdeme hranicni uzly
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
for k = Nuzel:-1:1 % urava matice A a rave strany B
    if (uzel_hra(k)==1) 
        B = B-A(:,k)*FcegD(xx1(k),xx2(k));
        A(k,:)=[]; A(:,k)=[]; B(k)=[];
        upom(k) = FcegD(xx1(k),xx2(k)); % doplneni hranicnich hodnot
    end;
end;
u = A\B; % reseni
zde = 1;
for k = 1:Nuzel
    if (uzel_hra(k)==0)
        upom(k) = u(zde); % doplneni reseni ve vnitrnich uzlech
        zde = zde+1;
    end;
end;
KresliReseni(elem,xx1,xx2,Nelem,upom);

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
z = eye(2)/5+x1^2;
function z = Fcep(x1,x2)
z = 0.1;
function z = Fcef(x1,x2)
z = 1-x1-x2;
function z = FcegD(x1,x2)
z = 1+x2/30+x1/30;
function z = FcegN(x1,x2)
z = 1+x1/10;




