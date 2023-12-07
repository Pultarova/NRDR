% MKP pro u_t = (au_x)_x + f, u(0)=OP1, u(L)=OP2
function abc
global elem
global uzly
L = 2; % interval (0,L)
OP1 = 0; % Dirichletova okrajova podminka vlevo
OP2 = -1; % Dirichletova okrajova podminka vpravo
Nx = 29; % pocet elementu
dx = L/Nx; % delka prostoroveho kroku
T = 4.4; % celkovy cas
Nt = 20; % pocet casovych kroku
dt = T/Nt; % delka casoveho kroku
uzly = linspace(0,L,Nx+1); % seznam uzlovych hodnot
elem = zeros(Nx,2); % seznam poradovych cisel uzlu na jednotlivych elementech
for k = 1:Nx 
    elem(k,:) = [k,k+1];
end;
U0 = zeros(Nx+1,1); % pocatecni stav (v case t=0)
for k = 1:Nx+1
    U0(k) = FceU0((k-1)*dx);
end;

VseU = zeros(Nt+1,Nx+1); % pole pro VSECHNA casova reseni
VseU(1,:) = U0'; % pocatecni reseni

for kt = 1:Nt % cylus pro casove kroky

A = zeros(Nx+1); % matice A tuhosti
M = zeros(Nx+1); % matice M hmotnosti
B = zeros(Nx+1,1); % vektor se zdrojem
for k = 1:Nx % jdeme po elementech
    % cisla bazovych funkci (= cisla uzlu) na tomto elementu:
    kde = elem(k,:); 
    % pricteni integralu z tohoto elementu na nektere pozice matic A,M:
    A(kde,kde) = A(kde,kde) + Integruj_a(k,kt*dt); 
    M(kde,kde) = M(kde,kde) + Integruj_1(k,kt*dt); 
    % pricteni integralu z tohoto elementu na nektere pozice vektoru B:
    B(kde) = B(kde) + Integruj_f(k,kt*dt);
end;
AA = M + dt*A; % implicitni schema
B = M*U0+dt*B;
B = B - AA(:,1)*OP1; % uplatneni leve OP
B = B - AA(:,Nx+1)*OP2; % uplatneni prave OP
AA = AA(2:Nx,2:Nx); % vlastni matice soustavy (zmensena o krajni uzly)
B = B(2:Nx); % vlastni prava strana (zmensena o krajni uzly)
U1 = AA\B; % reseni
U1 = [OP1;U1;OP2]; % pridani obou OP
U0 = U1; % posunuti
VseU(kt+1,:) = U1'; % ulozeni tohoto casoveho reseni
end;
cla; hold on;

% plot(uzly,U1); % kresleni jen posledniho reseni
% plot(uzly,uzly*0,'k.'); % kresleni uzlu

t = linspace(0,T,Nt+1); % sada casovych bodu
[xx,tt] = meshgrid(uzly,t); % sada vsech uzlu typu [bod,cas]
surf(xx,tt,VseU,'FaceAlpha',0.5) % kresleni celeho prubehu reseni
view(60,10); % uhel pohledu na graf

% --------------------------------------------------
function z = Integruj_a(k,t) % numericka integrace integralu s "a"
global elem
global uzly
un = elem(k,:);
x = uzly(un);
h = x(2)-x(1);
xs = (x(1)+x(2))/2;
z = Fcea(xs,t)*[1,-1;-1,1]/h;

function z = Integruj_1(k,t) % numericka integrace integralu s "1"
global elem
global uzly
un = elem(k,:);
x = uzly(un);
h = x(2)-x(1);
xs = (x(1)+x(2))/2;
z = [1/3,1/6;1/6,1/3]*h;

function z = Integruj_f(k,t)  % numericka integrace integralu s "f"
global elem
global uzly
un = elem(k,:);
x = uzly(un);
h = x(2)-x(1);
xs = (x(1)+x(2))/2;
z = Fcef(xs,t)*[1;1]*h/2;

function z = Fcea(x,t) % vodivost
z = (x-0.5)^2+0.001;
function z = Fcef(x,t) % zdroj tepla
z = 0.5+t/4-x*t/2;
z = 0;
if (x>1.5) z=1; end;
function z = FceU0(x) % pocatecni rozlozeni teploty
z = sin(x*pi*3/4);


