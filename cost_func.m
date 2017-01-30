% Aysel Akgemci 
% Fall 2016 Term Project 
% cost function
function result = cost_func(x)
%% initilizations
% variables
lm = x(1);            % magnet lenght
Dia = x(2);           % bore diameter
PP = x(3);            % # of pole pairs
pole = 2*PP;          % # of poles

% constant values
cu_density = 8960;    % densities in kg/m^3
iron_density = 7870;
magnet_density = 7550;
g = 1.5*10^-3;          % air gap
Pout = 50000;           % output power in W
q = 22000;              % stator electric loading in A/m
ns = 60;                % speed in rpm
Ea = 220;               % stator phase voltage in V
kw = 1;                 % winding factor
kf = 0.4;               % fill factor
kstk = 0.95;            % stacking factor
m = 3;                  % phase number
Br = 1.2;               % remenant flux density in T
Bsat = 1.4;             % saturation flux density of the core in T
Bt = 1.8;               % tooth flux density
mu_rec = 1.044;         % recoil permeability
J = 3e6;                % conductor current density in A/m^2
kml = 1.05;             % magnet leakage flux
mu0 = 4*pi*10^(-7);     % vacuum permeability
N_ph = 1;               % turns per pole per phase
S = 0.5;                % tooth&slot/pitch ratio( tooth width = slot width)
mag_arc = 0.67;
% minumum values determined before applying penalty 
min_tw = 0.01;
min_L = 0.1;
min_Di = 0.2;
%% calculations of penalty and dependent variables
N_slot = N_ph*pole*m;          % turns per slot
freq = ns*pole/120;            % frequency
we = 2*pi*freq;                % speed in electrical rad/s
wm = we/PP;                    % speed in mechanical rad/s
Tem = Pout/wm;                 % electromechanical torque

% calculation of air gap flux density
Am1 = (2/3)*pi*(Dia-2*g-lm)/pole;            % magnet area per unit of axial length
Pm01 = mu0*mu_rec*Am1/lm;                    % magnet internal permeance per unit of axial length
Pm1 = kml*Pm01;                              % magnet total permeance per unit of axial length

Ag1 = (2/3)*pi*(Dia-g)/pole+2*g;             % Air gap area per unit of axial length
Rg1 = g/(mu0*Ag1);                           % Air gap reluctance per unit of axial length

Bg = (Am1/Ag1)*(1/(1+Pm1*Rg1))*Br;           % Flat top value of air gap flux demsity
Bavg = 8*Bg*sin(pi/3)/pi^2;                  % average value of air gap flux density

%axial length calculation
L_axial = Tem*4*sqrt(2)/(pi^2*Bavg*q*kw*Dia^2);    % axial length of generator

    if L_axial<=min_L
    L_pen=0.1-L_axial;        % determination of penalty for length;
    else
    L_pen=0;
    end

%pole flux
Ap = pi*(Dia-g)*L_axial/pole;
phim = Bavg*Ap;
%back core length
hrbc = phim/(2*Bsat*L_axial);
hsbc = phim/(2*Bsat*L_axial);
%tooth length and width
tw = pi*(Dia)*(1-S)/(N_slot);
    if tw<=min_tw
    tw_pen=min_tw-tw;          % determination of penalty for tooth width;
    else
    tw_pen=0;
    end

h12 = 1.5e-3;
h13 = h12;
h1 = 0.5*(N_slot*tw/pi-Dia-g-2*h12-2*h13+sqrt((Dia+2*h12+2*h13+g-N_slot*tw/pi).^2+4*Dia*q/(J*kf)));

%outer and inner diameters
Di = Dia-2*(lm+g+hrbc);

    if Di<=min_Di
    Di_pen=min_Di-Di;        % determination of penalty for diameter;
    else
    Di_pen=0;
    end

Do = Dia+2*(hsbc+h1+h12+h13);

N = Ea/(4.44*kw*freq*phim);           % number of turns
Aco = q*pi*Dia/(6*N*J);
Lend =(N_ph+1)*(tw/2);
LMC = 2*L_axial+2*pi*(Dia+h1+2*h12+2*h13+2*g)/pole+4*Lend;

B_peak_tooth=(pi/2)*2*Bavg;         % peak value of tooth flux density

    if B_peak_tooth>=Bt
    Bt_pen=B_peak_tooth-Bt;         % determination of penalty tooth magnet flux;
    else
    Bt_pen=0;
    end
    
%% Mass calculation    
%copper mass
M_cu = 3*cu_density*N*LMC*Aco;
%iron mass
M_iron = iron_density*L_axial*(0.25*pi*(Do^2-(Dia+2*h1+2*h12+2*h13)^2)+N_slot*tw*(h1+h12+h13)+0.25*pi*((Dia-2*g-2*lm)^2-Di^2));
%magent mass
M_magnet = magnet_density*L_axial*mag_arc*pi*0.25*((Dia-2*g)^2-(Dia-2*g-2*lm)^2);
%penalty coefficients
L_pc=10000;
tw_pc=100000;
Di_pc=10000;
Bt_pc=10000;

penalty = tw_pc*tw_pen+L_pc*L_pen+Di_pc*Di_pen+Bt_pc*Bt_pen;
M_total = M_cu+M_iron+M_magnet;

result = M_total + penalty;

end
