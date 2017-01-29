function M=RF_mass_NoStrc(x)

global j M_t

Lm=x(1);
D=x(2);
ps=x(3);
p=2*ps;

ro_cu=1.72e-8;

cu_density=8960;
iron_density=7870;
magnet_density=7550;

mag_arc=0.67;
min_tw=0.01;
min_L=0.1;
min_Di=.2;

%penalty coefficients
L_pc=10000;
tw_pc=100000;
Di_pc=10000;
Bt_pc=10000;


Pout=50e3;
q=22000;
n=60;
f=n*p/120;
we=2*pi*f;
wm=we/(p/2);
Tem=Pout/wm;
g=1.5e-3;
Bbc=1.4;
Bt=1.8;
Ea=220;
Ia=Pout/(3*Ea);
kw=1;
m=3;
%Lend=5e-3;

%magnet length
Br=1.2;                                             %magnet Remanent Flux Density
urec=1.044;                                         %magnet Recoil permeability
kml=1.05;                                           %magnet leakage flux
u0=4*pi*1e-7;

%core losses coefficients
kh=0.0383;
ke=5.8e-5;
nn=2;

so=.5;
Ns_p_ph=1;

N_slot=Ns_p_ph*p*m;
J=3e6;
Kcu=0.4;

%Calculation of air gap flux density
Am1=(2/3)*pi*(D-2*g-Lm)/p;                          %magnet area per unit of axial length
Pm01=u0*urec*Am1/Lm;                                %magnet internal permeance per unit of axial length
Pm1=kml*Pm01;                                       %magnet total permeance per unit of axial length

Ag1=(2/3)*pi*(D-g)/p+2*g;                           %Air gap area per unit of axial length
Rg1=g/(u0*Ag1);                                     %Air gap reluctance per unit of axial length

Bg=(Am1/Ag1)*(1/(1+Pm1*Rg1))*Br;                    %Flat top value of air gap flux demsity
Bavg=8*Bg*sin(pi/3)/pi^2;                           %average value of air gap flux density

%axial length calculation
Lstk=Tem*4*sqrt(2)/(pi^2*Bavg*q*kw*D^2);             %axial length of generator
kr=Lstk/D;

if Lstk<=min_L
    L_pen=0.1-Lstk;
else
    L_pen=0;
end

%pole flux
Ap=pi*(D-g)*Lstk/p;
phim=Bavg*Ap;


%back core length
hrbc=phim/(2*Bbc*Lstk);
hsbc=phim/(2*Bbc*Lstk);



%tooth length and width

tw=pi*(D)*(1-so)/(N_slot);
if tw<=min_tw
    tw_pen=min_tw-tw;
else
    tw_pen=0;
end

h12=1.5e-3;
h13=h12;
h1=.5*(N_slot*tw/pi-D-g-2*h12-2*h13+sqrt((D+2*h12+2*h13+g-N_slot*tw/pi).^2+4*D*q/(J*Kcu)));


%outer and inner diameters
Di=D-2*(Lm+g+hrbc);

if Di<=min_Di
    Di_pen=min_Di-Di;
else
    Di_pen=0;
end

Do=D+2*(hsbc+h1+h12+h13);

%Number of turns
N=Ea/(4.44*kw*f*phim);
N_p=(N/(p/2));                                  
Nph=(p/2)*round(N_p);

%Phase resistance
Aco=q*pi*D/(6*N*J);
A_slot=N_p*Aco/Kcu;

Lend=(Ns_p_ph+1)*(tw/2);
%Lend=10e-3;
LMC=2*Lstk+2*pi*(D+h1+2*h12+2*h13+2*g)/p+4*Lend;
Rph=ro_cu*LMC*N/Aco;

%phase inductance
Am=Am1*Lstk;
Ag=Ag1*Lstk;
Rg=g/(u0*Ap);
Rm=Lm/(u0*urec*Ap);
Lph=(p/2)*N_p^2/(2*(Rg+Rm));
xs=1.5*2*pi*f*Lph;


%Copper and core losses
Pcu=3*Rph*Ia^2;                                             %Copper losses

B_peak_tooth=(pi/2)*2*Bavg;                                 %peak value of tooth flux density

if B_peak_tooth>=Bt
    Bt_pen=B_peak_tooth-Bt;
else
    Bt_pen=0;
end


pc_tooth=kh*f*B_peak_tooth^nn+ke*f^2*B_peak_tooth^2;         %tooth core losses
pc_backcore=kh*f*Bbc^nn+ke*f^2*Bbc^2;                        %back core core losses

Mt=iron_density*Lstk*N_slot*tw*(h1+h12+h13);                %tooth mass
Mbc=iron_density*Lstk*(0.25*pi*(Do^2-(D+2*h1+2*h12+2*h13)^2));  %back core mass

Pcore=Mt*pc_tooth+Mbc*pc_backcore;                         %total core loasses

Ploss=Pcu+Pcore;                                           %total losses

%performance calculation
pf=Ea/sqrt(Ea^2+(xs*Ia)^2);                               %power factor
efficiency=Pout/(Pout+Ploss);                             %efficiency

%copper mass
M_cu=3*cu_density*N*LMC*Aco;

%iron mass
M_iron=iron_density*Lstk*(0.25*pi*(Do^2-(D+2*h1+2*h12+2*h13)^2)+N_slot*tw*(h1+h12+h13)+0.25*pi*((D-2*g-2*Lm)^2-Di^2));

%magent mass
M_magnet=magnet_density*Lstk*mag_arc*pi*0.25*((D-2*g)^2-(D-2*g-2*Lm)^2);


penalty_fun=tw_pc*tw_pen+L_pc*L_pen+Di_pc*Di_pen+Bt_pc*Bt_pen;

M_active=M_cu+M_iron+M_magnet;



%% Objective function
M=M_active+penalty_fun;

j=j+1;
M_t(j)=M;



end 




