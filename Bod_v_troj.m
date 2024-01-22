
function z = Bod_v_troj(v1,v2,v3,b)
% Fce vrati 1, jestlize bod b lezi uvnitr trojuhelnika s vrcholy v1,v2,v3.
% Fce vrati 0, jestlize bod b nelezi uvnitr trojuhelnika s vrcholy v1,v2,v3.
% Bod b a vrcholy v1,v2,v3 jsou radkove vektory o 2 slozkach.
z = 1; % predpokladame, ze bod b je v trojuhelniku
if (det([v2-v1;v3-v1])*det([v2-v1;b-v1])<0) z = 0;
elseif (det([v3-v2;v1-v2])*det([v3-v2;b-v2])<0) z = 0;
elseif (det([v1-v3;v2-v3])*det([v1-v3;b-v3])<0) z = 0;
end;