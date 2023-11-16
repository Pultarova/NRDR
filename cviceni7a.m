% MKP pro -(au')'+pu=f, 
% u(0)= OP1, u(L)=OP2 Dirichlet
% nebo pro
% u(0) = OP1, u'(L)+OP2N1*u(L)=OP2N2 Newton
function abc
global elem
global uzly
L = 2; % interval (0,L)
OP1 = -1; 
OP2 = 0;
OP2N1 = 1;
OP2N2 = -1;
ZDA_OP2_D_N = 2; % D..1, N..2
N = 31; % pocet elementu
uzly = linspace(0,L,N+1); % seznam uzlovych hodnot
N2 = 20;
uzly2 = linspace(0.4,0.6,N2);
% pridane uzly (neosetreny duplikovane hodnoty!!)
uzly = (sort([uzly,uzly2]));
N = N + N2;
elem = zeros(N,2); 
for k = 1:N % seznam uzlu na jednotlivych elementech
    elem(k,:) = [k,k+1];
end;
A = zeros(N+1); % matice soustavy
B = zeros(N+1,1); % prava strana soustavy
for k = 1:N
    % cisla bazovych funkci (= cisla uzlu) na tomto elementu:
    kde = elem(k,:); 
    % pricteni integralu z tohoto elementu na nektere pozice matice A:
    A(kde,kde) = A(kde,kde) + Integruj_a(k) + Integruj_p(k); 
    % pricteni integralu z tohoto elementu na nektere pozice vektoru B:
    B(kde) = B(kde) + Integruj_f(k);
end;
B = B - A(:,1)*OP1; % uplatneni leve OP
if (ZDA_OP2_D_N==1) % Dirichletova OP vpravo
    B = B - A(:,N+1)*OP2; % uplatneni prave OP
    A = A(2:N,2:N); % vlastni matice soustavy 
    B = B(2:N); % vlastni prava strana
    U = A\B; % reseni
    cla; hold on;
    plot(uzly,[OP1;U;OP2]); % kresleni reseni
    plot(uzly,uzly*0,'k.'); % kresleni uzlu
else % Newtonova OP vpravo
    A(N+1,N+1) = A(N+1,N+1)+Fcea(L)*OP2N1;
    B(N+1) = B(N+1)+Fcea(L)*OP2N2;
    A = A(2:N+1,2:N+1); % vlastni matice soustavy 
    B = B(2:N+1); % vlastni prava strana
    U = A\B; % reseni
    cla; hold on;
    plot(uzly,[OP1;U]); % kresleni reseni
    plot(uzly,uzly*0,'k.'); % kresleni uzlu
end

% --------------------------------------------------
function z = Integruj_a(k) % numericka integrace integralu s "a"
global elem
global uzly
un = elem(k,:);
x = uzly(un);
h = x(2)-x(1);
xs = (x(1)+x(2))/2;
z = Fcea(xs)*[1,-1;-1,1]/h;

function z = Integruj_p(k) % numericka integrace integralu s "p"
global elem
global uzly
un = elem(k,:);
x = uzly(un);
h = x(2)-x(1);
xs = (x(1)+x(2))/2;
z = Fcep(xs)*[1/3,1/6;1/6,1/3]*h;

function z = Integruj_f(k)  % numericka integrace integralu s "f"
global elem
global uzly
un = elem(k,:);
x = uzly(un);
h = x(2)-x(1);
xs = (x(1)+x(2))/2;
z = Fcef(xs)*[1;1]*h/2;

function z = Fcea(x)
z = (x-0.5)^2+0.001;
function z = Fcep(x)
z = 0;
function z = Fcef(x)
z = 1;

