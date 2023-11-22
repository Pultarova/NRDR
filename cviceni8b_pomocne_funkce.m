% MKP pro 2d ulohu
% -div.(a grad u)+pu=f na (0,L1)x(0,L2), 
% na hranici Dirichlet u=gD nebo Neumann a*grad(u)*n=gN

function abc
L1 = 1; 
L2 = 1;
N1 = 12;
N2 = 12;



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
        [u(t(1)),u(t(2));u(t(3)),u(t(3))],'FaceColor','blue','FaceAlpha',0.5);
end;

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
z = 1;
function z = Fcep(x1,x2)
z = 0;
function z = Fcef(x1,x2)
z = 1;
function z = FcegD(x1,x2)
z = 1+x1/10-x2^2/7;
z = 0;
function z = FcegN(x1,x2)
z = 1+x1/10-x2^2/7;
z = 0;




