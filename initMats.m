% w = wM*[x;u;alpha;beta]+wm

xpM = zeros(2,3,s.N);
xpm = zeros(2,1,s.N);


wM = zeros(2,3,s.N);
wm = zeros(2,1,s.N);

uK = zeros(1,2,s.N);
uk = zeros(1,1,s.N);


muM = zeros(2,3,s.N);
mum = zeros(2,1,s.N);
muZ = zeros(2,1,s.N);


lambdaM = zeros(length(s.g),3,s.N);
lambdam = zeros(length(s.g),1,s.N);
lambdaZ = zeros(length(s.g),1,s.N);

zetaM = zeros(1,3,s.N);
zetam = zeros(1,1,s.N);
zetaZ = zeros(1,1,s.N);

kappaK = cell(1,s.N);
kappak = cell(1,s.N);
kappaZ = cell(1,s.N);

rhoK = zeros(1,2,s.N);
rhok = zeros(1,1,s.N);
rhoZ = zeros(1,1,s.N);