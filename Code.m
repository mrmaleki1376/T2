clc
clear All
close All
%%bename khoda
format long
%n_reza=0:1:6500;
m_e=1250;
Mu=1.0;
w_D=0.55;
n_N=5500;
k_SI=0.1;
M_M1=250;
M_M2=250;
n_M1=2000;
n_M2=3000;
n_idle=750;
n_max=6500;
n_reza=0:1:6500;
w=n_reza*0.1047;
N_N=110;
n_N1=576; %5500rpm
i_gx =[14.5 , 8.137 , 5.328 , 4.008 , 3.246 , 2.792 ] ;
%225/45R17
D_W=.6343;
R_d=0.485*D_W;
A= 2.18;
C_d=0.31;
ro=1.2;
i_f=3.77;
m_t=1800;
g=9.8;
G_t=m_t*g;
p_t=241316.5052585;
eta_tr=0.93;
     
syms e_n e_M n_M M_M M_N
E=e_n*e_M;
E1=(n_N/n_M)*(M_M/M_N);
N_nmax = ((-0.6+k_SI)*(n_max/n_N)+(1.6-k_SI))*N_N;
%N_nmax_CI = ((-1.2+k_CI)*(n_max/n_N)+(2.2-k_CI))*N_N;

n=750:6500;
%n=v_t1*((i_gx(1)*60)/(R_d*2*pi))
v_t1=n*(R_d*2*pi)/(i_gx(1)*60);
v_t2=n*(R_d*2*pi)/(i_gx(2)*60);
v_t3=n*(R_d*2*pi)/(i_gx(3)*60);
v_t4=n*(R_d*2*pi)/(i_gx(4)*60);
v_t5=n*(R_d*2*pi)/(i_gx(5)*60);
v_t6=n*(R_d*2*pi)/(i_gx(6)*60);

N_M1=(M_M1*n_M1)/9550;
N_M2=(M_M2*n_M2)/9550;
M_N=9550*(N_N/n_N);
p=(M_N*(n_M2-n_N))/(n_N*(M_N-M_M2));
%if n ==(750:1:2000)
%   N=(N_M1/(n_M1-n_idle))*(n-n_idle);
%elseif n==(2001:1:3000);
  %  N=((N_M2-N_M1)/(n_M2-n_M1))*(n-n_M1)+ N_M1;
%elseif n== (3001:1:5500)
  %  N=(((M_N-M_M2)/(n_N-n_M2)^p)*(n-n_M2)^p+ N_M2)*(n/9550);
%elseif n== (5501:1:6500)
 %   N=((N_nmax-N_N)/(n_max-n_N)^2)*((n-n_N)^2)+ N_N;
%end
N = ((N_M1/(n_M1- n_idle )).*(n - n_idle )).*(n>=750 & n<= 2000) + (((N_M2-N_M1)/(n_M2-n_M1)).*(n-n_M1)+ N_M1).*( n>2000 & n<=3000 )+((((M_N-M_M2)/(n_N-n_M2).^p).*((n-n_M2).^p)+ M_M2).*(n/9550)).* ( n> 3000 & n<=5500) + (((N_nmax-N_N)/(n_max-n_N).^2).*((n-n_N).^2)+ N_N).* ( ( n>5500)) ;
M = 10000*((((N_M1/(n_M1- n_idle )).*(n - n_idle ))./n).*(n>=750 & n<= 2000) + ((((N_M2-N_M1)/(n_M2-n_M1)).*(n-n_M1)+ N_M1)./n).*( n>2000 & n<=3000 )+(((((M_N-M_M2)/(n_N-n_M2).^p).*((n-n_M2).^p)+ M_M2).*(n/9550))./n).* ( n> 3000 & n<=5500) + ((((N_nmax-N_N)/(n_max-n_N).^2).*((n-n_N).^2)+ N_N)./n).* ( ( n>5500)));
M_1= 10000*((((N_M1/(n_M1- n_idle )).*((v_t1*((i_gx(1)*60)/(R_d*2*pi))) - n_idle ))./(v_t1*((i_gx(1)*60)/(R_d*2*pi)))).*(n>=750 & n<= 2000) + ((((N_M2-N_M1)/(n_M2-n_M1)).*((v_t1*((i_gx(1)*60)/(R_d*2*pi)))-n_M1)+ N_M1)./(v_t1*((i_gx(1)*60)/(R_d*2*pi)))).*( n>2000 & n<=3000 )+(((((M_N-M_M2)/(n_N-n_M2).^p).*(((v_t1*((i_gx(1)*60)/(R_d*2*pi)))-n_M2).^p)+ M_M2).*((v_t1*((i_gx(1)*60)/(R_d*2*pi)))/9550))./(v_t1*((i_gx(1)*60)/(R_d*2*pi)))).* ( n> 3000 & n<=5500) + ((((N_nmax-N_N)/(n_max-n_N).^2).*(((v_t1*((i_gx(1)*60)/(R_d*2*pi)))-n_N).^2)+ N_N)./(v_t1*((i_gx(1)*60)/(R_d*2*pi)))).* ( ( n>5500)));
M_2= 10000*((((N_M1/(n_M1- n_idle )).*((v_t2*((i_gx(2)*60)/(R_d*2*pi))) - n_idle ))./(v_t2*((i_gx(2)*60)/(R_d*2*pi)))).*(n>=750 & n<= 2000) + ((((N_M2-N_M1)/(n_M2-n_M1)).*((v_t2*((i_gx(2)*60)/(R_d*2*pi)))-n_M1)+ N_M1)./(v_t2*((i_gx(2)*60)/(R_d*2*pi)))).*( n>2000 & n<=3000 )+(((((M_N-M_M2)/(n_N-n_M2).^p).*(((v_t2*((i_gx(2)*60)/(R_d*2*pi)))-n_M2).^p)+ M_M2).*((v_t2*((i_gx(2)*60)/(R_d*2*pi)))/9550))./(v_t2*((i_gx(2)*60)/(R_d*2*pi)))).* ( n> 3000 & n<=5500) + ((((N_nmax-N_N)/(n_max-n_N).^2).*(((v_t2*((i_gx(2)*60)/(R_d*2*pi)))-n_N).^2)+ N_N)./(v_t2*((i_gx(2)*60)/(R_d*2*pi)))).* ( ( n>5500)));
M_3= 10000*((((N_M1/(n_M1- n_idle )).*((v_t3*((i_gx(3)*60)/(R_d*2*pi))) - n_idle ))./(v_t3*((i_gx(3)*60)/(R_d*2*pi)))).*(n>=750 & n<= 2000) + ((((N_M2-N_M1)/(n_M2-n_M1)).*((v_t3*((i_gx(3)*60)/(R_d*2*pi)))-n_M1)+ N_M1)./(v_t3*((i_gx(3)*60)/(R_d*2*pi)))).*( n>2000 & n<=3000 )+(((((M_N-M_M2)/(n_N-n_M2).^p).*(((v_t3*((i_gx(3)*60)/(R_d*2*pi)))-n_M2).^p)+ M_M2).*((v_t3*((i_gx(3)*60)/(R_d*2*pi)))/9550))./(v_t3*((i_gx(3)*60)/(R_d*2*pi)))).* ( n> 3000 & n<=5500) + ((((N_nmax-N_N)/(n_max-n_N).^2).*(((v_t3*((i_gx(3)*60)/(R_d*2*pi)))-n_N).^2)+ N_N)./(v_t3*((i_gx(3)*60)/(R_d*2*pi)))).* ( ( n>5500)));
M_4= 10000*((((N_M1/(n_M1- n_idle )).*((v_t4*((i_gx(4)*60)/(R_d*2*pi))) - n_idle ))./(v_t4*((i_gx(4)*60)/(R_d*2*pi)))).*(n>=750 & n<= 2000) + ((((N_M2-N_M1)/(n_M2-n_M1)).*((v_t4*((i_gx(4)*60)/(R_d*2*pi)))-n_M1)+ N_M1)./(v_t4*((i_gx(4)*60)/(R_d*2*pi)))).*( n>2000 & n<=3000 )+(((((M_N-M_M2)/(n_N-n_M2).^p).*(((v_t4*((i_gx(4)*60)/(R_d*2*pi)))-n_M2).^p)+ M_M2).*((v_t4*((i_gx(4)*60)/(R_d*2*pi)))/9550))./(v_t4*((i_gx(4)*60)/(R_d*2*pi)))).* ( n> 3000 & n<=5500) + ((((N_nmax-N_N)/(n_max-n_N).^2).*(((v_t4*((i_gx(4)*60)/(R_d*2*pi)))-n_N).^2)+ N_N)./(v_t4*((i_gx(4)*60)/(R_d*2*pi)))).* ( ( n>5500)));
M_5= 10000*((((N_M1/(n_M1- n_idle )).*((v_t5*((i_gx(5)*60)/(R_d*2*pi))) - n_idle ))./(v_t5*((i_gx(5)*60)/(R_d*2*pi)))).*(n>=750 & n<= 2000) + ((((N_M2-N_M1)/(n_M2-n_M1)).*((v_t5*((i_gx(5)*60)/(R_d*2*pi)))-n_M1)+ N_M1)./(v_t5*((i_gx(5)*60)/(R_d*2*pi)))).*( n>2000 & n<=3000 )+(((((M_N-M_M2)/(n_N-n_M2).^p).*(((v_t5*((i_gx(5)*60)/(R_d*2*pi)))-n_M2).^p)+ M_M2).*((v_t5*((i_gx(5)*60)/(R_d*2*pi)))/9550))./(v_t5*((i_gx(5)*60)/(R_d*2*pi)))).* ( n> 3000 & n<=5500) + ((((N_nmax-N_N)/(n_max-n_N).^2).*(((v_t5*((i_gx(5)*60)/(R_d*2*pi)))-n_N).^2)+ N_N)./(v_t5*((i_gx(5)*60)/(R_d*2*pi)))).* ( ( n>5500)));
M_6= 10000*((((N_M1/(n_M1- n_idle )).*((v_t6*((i_gx(6)*60)/(R_d*2*pi))) - n_idle ))./(v_t6*((i_gx(6)*60)/(R_d*2*pi)))).*(n>=750 & n<= 2000) + ((((N_M2-N_M1)/(n_M2-n_M1)).*((v_t6*((i_gx(6)*60)/(R_d*2*pi)))-n_M1)+ N_M1)./(v_t6*((i_gx(6)*60)/(R_d*2*pi)))).*( n>2000 & n<=3000 )+(((((M_N-M_M2)/(n_N-n_M2).^p).*(((v_t6*((i_gx(6)*60)/(R_d*2*pi)))-n_M2).^p)+ M_M2).*((v_t6*((i_gx(6)*60)/(R_d*2*pi)))/9550))./(v_t6*((i_gx(6)*60)/(R_d*2*pi)))).* ( n> 3000 & n<=5500) + ((((N_nmax-N_N)/(n_max-n_N).^2).*(((v_t6*((i_gx(6)*60)/(R_d*2*pi)))-n_N).^2)+ N_N)./(v_t6*((i_gx(6)*60)/(R_d*2*pi)))).* ( ( n>5500)));
N_reza=(((N_N/n_N1)*w)+((N_N/(n_N1^2))*(w.^2))+(((-N_N)/((n_N1)^3))*(w.^3)));
M_reza=1000*(((N_N/n_N1))+((N_N/(n_N1^2))*(w))+(((-N_N)/((n_N1)^3))*(w.^2)));

figure(1)
   plot (n_reza , N_reza ,'--', 'linewidth' , 2.5)
   hold on;
   plot ( n_reza , M_reza,'--', 'linewidth' , 2.5)
   hold on;  
   plot(n , N, 'linewidth' , 2.5);
   hold on;
   plot(n , M , 'linewidth' , 2.5);
   hold on;
  
   %%Forces %%

for i=1:6
syms   teta


%%%%%%%%%%%%%%%%%%%%%%F_D
v=n_reza * 2 * pi * R_d * 3.6 /( i_gx(i)* 60);
F_D1= M_reza*i_gx(i)* eta_tr / R_d;
F_D2=m_t*g*Mu*w_D;
figure(9)


plot(v,F_D1,'linewidth' , 3) 
hold on;
end
for i=1:6
syms   teta


%%%%%%%%%%%%%%%%%%%%%%F_D
v=n * 2 * pi * R_d * 3.6 /( i_gx(i)* 60);
F_D1= M*i_gx(i)* eta_tr / R_d;
F_D2=m_t*g*Mu*w_D;
figure(2)


plot(v,F_D1,'linewidth' , 3) 
hold on;
end
%%%%%%%%%%%%%%%%%%%%%%%F_aero

V=0:1:300;

F_aero= 0.5 * ro * C_d * A * ((V/3.6).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%F_roll
for teta=0:0.1:0.6
alpha=atand(teta);
f_r=0.0008*(5.1+((5.5+90*G_t)/(p_t))+((1100+0.0388*G_t)/(p_t))*((V/3.6).^2));
F_roll=m_t*g*cosd(alpha)*f_r;
%N_roll=(m_t*g*cosd(alpha)*f_r).*v;



%%%%%%%%%%%%%%%%%%%%%F_slope

F_slope= m_t * g * sind(alpha);
%N_slope=(m_t * g * sind(alpha)).* v ;
%%%%%%%%%%%%%%%%%%%F_tot
F_tot=(F_aero+F_roll+F_slope);
figure(2)
       plot(V , F_tot,'--','linewidth' , 1.5);
       
hold on;
N_tot=F_tot .* V;
figure(9)
       plot(V , F_tot,'--','linewidth' , 1.5);
       
hold on;
end
%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%fuel consumption
%syms ge_1 ge_2
ge_min=225;
n0=2500;
M0=212.5;
n_alpha=2000;
n_P1=4500;
M_P2=50;
ge_P1=245;
ge_P2=305;
G_V=0.6;
ro_f=750;
M_idle=0.5;
M_M=250;
ge=226;

%syms nge Mge x
%a1=(33.3/M_M)*log10((9.55e06*G_V*ro_f*(((M_P2-M0)^2))/(M_idle*n_idle*(ge_P2-ge_min)*Mge^2))-1);
%ge_1=(1-0.5^(0.2*Mge))*((ge_P1-ge_min)/((n_P1-n0)^2))*((nge-(n_alpha+((n0-n_alpha)/M0)))^2)+ ge_min;
%ge_2=((ge_P2-ge_min)/((M_P2-M0)^2))*((Mge-M0)^2)*(1+0.5^(a1*(Mge-0.1*M_M)));
%x=solve('((((ge_P2-ge_min)/((M_P2-M0)^2))*((Mge-M0)^2)*(1+0.5^(a1*(Mge-0.1*M_M))))+((1-0.5^(0.2*Mge))*((ge_P1-ge_min)/((n_P1-n0)^2))*((nge-(n_alpha+((n0-n_alpha)/M0)))^2)+ ge_min))=ge');
%ge=ge_1+ge_2;
%ge=226;
%plot(nge ,Mge )




%%%%%%%%%%%%%%%%%%%%%%%

%%%%226
MUe226 = 2500;
SIG2262 = 450.9;
SIG2261 = 432;
N226=2090:1:2910; 
x2261 = (N226-MUe226)/SIG2261;
x2262 = (N226-MUe226)/SIG2262;
A2261=     2.887 *x2261.^4 + -1.615e-14*x2261.^3 +   -27.31 *x2261.^2 + 4.167e-14 *x2261 +   234.3;
A2262=    0.7678 *x2262.^4 +3.777e-14*x2262.^3 +   23.96 *x2262.^2 + -4.922e-14 *x2262 +    188.2  ;
%%%%230
MUe230 = 2683 ;
SIG230 =837.7;
N230=1500:1:3860; 
x230 = (N230-MUe230)/SIG230;
A230=12.4*x230.^4 + 7.204 *x230.^3 +     0.946*x230.^2 + 3.3*x230 +174;
%%%%240
MUe240 = 2922 ;
SIG240 =1173;
N240=1250:1:4410; 
x240 = (N240-MUe240)/SIG240;
A240=17.67 *x240.^4 + 14.41*x240.^3 +     -1.606 *x240.^2 +  8.868 *x240 +150.3;
%%%%260
MUe260 = 3225 ;
SIG260 =1474;
N260=1050:1:5300; 
x260 = (N260-MUe260)/SIG260;
A260=3.203*x260.^4 + 7.655*x260.^3 +     17.49 *x260.^2 +15.79 *x260 +114.3;
%%%%290
MUe290 = 3500 ;
SIG290 =1658;
N290=1000:1:6030; 
x290 = (N290-MUe290)/SIG290;
A290=8.708*x290.^4 +  10.79 *x290.^3 +     1.554  *x290.^2 +8.254  *x290 +76.59;
%%%%340
MUe340 = 3500 ;
SIG340 =1658;
N340=1000:1:6000; 
x340 = (N340-MUe340)/SIG340;
A340=3.42*x340.^4 + 5.152*x340.^3 +     3.163  *x340.^2 +6.045  *x340 +39.09;
%%%%1000
MUe1000 = 3500 ;
SIG1000 =1658;
N1000=1000:1:6000; 
x1000 = (N1000-MUe1000)/SIG1000;
A1000= 0.04231*x1000.^4 + 0.01843*x1000.^3    -0.1949 *x1000.^2 +0.938  *x1000 +18.74;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3)
plot(n , M,'R','linewidth' , 2.5)
hold on;
plot(N226,A2261,'k -.','linewidth' , 1.5 )
hold on;
plot(N1000,A1000,' -.','linewidth' , 1.5 )   
hold on;
plot(N340,A340,' -.','linewidth' , 1.5 )   
hold on;
plot(N290,A290,' -.' ,'linewidth' , 1.5)   
hold on;
plot(N260,A260,' -.','linewidth' , 1.5 )   
hold on;
plot(N230,A230,' -.','linewidth' , 1.5 )
hold on;
plot(N240,A240,' -.' ,'linewidth' , 1.5)
hold on;
plot(N226,A2262,'k -.','linewidth' , 1.5)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%First gear ratio calculation
%syms l_f l_r l h_c
%For front wheel drive
%tand(alpha_max)=(((l_r*Mu)-(l_f*f_r))/(l+(Mu+f_r)*h_c));
%For rear wheel drive
%tand(alpha_max)=(((l_f*Mu)-(l_r*f_r))/(l+(Mu+f_r)*h_c));
%M_M*i_t1*eta_tr/R_d>m_g*(g*sind(alpha)+g*cosd(alpha)*f_r+a_start)
%C=i_t1>3.7*(m_g*R_d/M_M);




%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_start=[1500,2000,2500];
c=[0.5,0.75];
for i=1:1:3
    for j=1:1:2
i_t1=12.2:0.1:17.8;
m_tu=1400;
m_teu=1500;
M_Mc=250;
s_tu=3;
s_teu=0.2;
alpha_u=0;
alpha_eu=5;
t_clutch_max=0.7*((5.5*1e10*((M_Mc).^(2/3))*((i_t1/R_d*n_start(i)).^2))/((2.5*m_tu*((1.5+(g*sind(alpha_u)))/1.4)*s_tu*c(j))+(m_teu*((1.5+(g*sind(alpha_eu)))/1.4)*s_teu*(1-c(j)))))/1e16;
figure(4)
    plot(i_t1 , t_clutch_max,'linewidth' , 2.5 );
hold on;
    end
    hold on;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n_vmax=[5000 , 5500 ,6500]
N_vmax=(((((M_N-M_M2)/(n_N-n_M2).^p).*((n_vmax-n_M2).^p)+ M_M2).*(n_vmax/9550)).*(n_vmax<n_N)+(((N_nmax-N_N)/(n_max-n_N).^2).*((n_vmax-n_N).^2)+ N_N).*(n_vmax>=n_N));
v_max=(45.7.*(N_vmax.^(1/3)))-(71.6*C_d*A)-0.019*(m_e+100)+59.7;
f_r_vmax=0.0008*(5.1+((5.5+90*G_t)/(p_t))+((1100+0.0388*G_t)/(p_t))*((v_max/3.6).^2));
i_tvmax=(R_d*n_vmax)/2.653*v_max;
v_vmax=0.86*50000*(n * 2 * pi * R_d * 3.6 /( i_tvmax * 60));
N_w=0.95*(((0.5*ro*C_d*A.*(v_vmax/3.6).^2)+(m_t*g*f_r_vmax)).*(v_vmax/3.6))/1000;
F_aero_vmax= 0.5 * ro * C_d * A * ((v_vmax/3.6).^2);
F_roll_vmax=m_t*g*f_r_vmax;
N_tot_vmax=2.9*(F_aero_vmax+F_roll_vmax).*v_vmax/10000;
N_v_max =0.9*(((N_M1/(n_M1- n_idle )).*(n - n_idle )).*(n>=750 & n<= 2000) + (((N_M2-N_M1)/(n_M2-n_M1)).*(n-n_M1)+ N_M1).*( n>2000 & n<=3000 )+((((M_N-M_M2)/(n_N-n_M2).^p).*((n-n_M2).^p)+ M_M2).*(n/9550)).* ( n> 3000 & n<=5500) + (((N_nmax-N_N)/(n_max-n_N).^2).*((n-n_N).^2)+ N_N).* ( ( n>5500)));
figure(5)
plot ( v_vmax , N_v_max, 'linewidth' , 2.5)
hold on;

end
hold on;

figure(5)
plot ( v_vmax , N_w , '-.', 'linewidth' , 1.5)

hold on;
axis([50 250 0 120])
grid on
title('Three possibilities of engine speed n_v_m_a_x at maximum vehicle speed v_m_a_x')
xlabel('v[km/h]')
ylabel('N_w[kW]')
legend( {'eco' ,'normal','dynamic' ,'N_r_e_s_i_s_t'} ,'location'  , 'bestoutside' )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Selection of intermediate gears (and overdrive)
k=5;
b=3.5;
n_vmax_2=6000;
v_max_2=204;
v_N1=44;
a=1046;
%v_N=(R_d*n_vmax_2)/2.653*i_t(i);
%a=((((n_N/n_vmax_2)*v_max_2)/((n_N^(-((b-1^(k-1))-(b^(k-1)))/(b^(k-2))))*(v_N1^(((b-1)^(k-1))/(b^(k-1))))))^((b^(k-2))/(((b-1)^(k-1))-(b^(k-1)))));
n_12=a*(v_N1^(1/b));
i_t2=(R_d*n_12)/(2.653*v_N1);
v_N2=(R_d*n_N)/(2.653*i_t2);
n_23=a*(v_N2^(1/b));
i_t3=(R_d*n_23)/(2.653*v_N2);
v_N3=(R_d*n_N)/(2.653*i_t3);
n_34=a*(v_N3^(1/b));
i_t4=(R_d*n_34)/(2.653*v_N3);
v_N4=(R_d*n_N)/(2.653*i_t4);
n_45=a*(v_N4^(1/b));
i_t5=(R_d*n_45)/(2.653*v_N4);
v_N5=(R_d*n_N)/(2.653*i_t5);
n_56=a*(v_N5^(1/b));
i_t6=(R_d*n_56)/(2.653*v_N5);
v_N6=(R_d*n_N)/(2.653*i_t6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_n_N=25:240;
y_vmax=0:6000;
x_n_N2=25:212;
y_v1=3084:5500;
y_v2=3639:5500;
y_v3=4094:5500;
y_v4=4454:5500;
y_v5=4731:5500;

MUe1 = 22;
SIG1 = 31.11;
V1=0:44; 
x1 = (V1-MUe1)/SIG1;
N11=   3889*x1 +     2750;

MUe2 = 39.2;
SIG2 = 55.44;
V2=0:78.4; 
x2 = (V2-MUe2)/SIG2;
N22=   3889*x2 +     2750;

MUe3 = 59.3 ;
SIG3 = 83.86;
V3=0:118.6; 
x3 = (V3-MUe3)/SIG3;
N33=   3889*x3 +     2750;

MUe4 = 79.65 ;
SIG4 = 112.6;
V4=0:159.3; 
x4 = (V4-MUe4)/SIG4;
N44=   3889*x4 +     2750;

MUe5 = 98.35 ;
SIG5 =139.1;
V5=0:212; 
x5 = (V5-MUe5)/SIG5;
N55=   3889*x5 +   2750 ;

MUe6 = 114.3 ;
SIG6 =161.7;
V6=0:228.7; 
x6 = (V6-MUe6)/SIG6;
N66=   3889*x6 +     2750;
V7=0:1:250; 
N77=   a*V7.^(1/b);

figure(6)
plot(V1,N11,'linewidth' , 2.5)
hold on;
plot(V2,N22,'-.','linewidth' , 2.5)   
hold on;
 plot(V3,N33 ,'-.','linewidth' , 2.5)   
hold on;
 plot(V4,N44,'-.','linewidth' , 2.5 )   
hold on;
 plot(V5,N55,'linewidth' , 2.5) ;  
hold on;
 plot(V6,N66,'-.','linewidth' , 2.5);
hold on;
 plot(V7,N77,'linewidth' , 1.5) ;  
hold on;
 plot(x_n_N,ones(size(x_n_N))*5500 ,'k','linewidth' , 0.9);
hold on;
 plot(x_n_N2,ones(size(x_n_N2))*6000 ,'k -.','linewidth' , 0.9);
hold on;
 plot(ones(size(y_vmax))*212,y_vmax ,'k -.','linewidth' , 0.9);
 hold on;
 plot(ones(size(y_v1))*44,y_v1 ,'k -.','linewidth' , 0.8);
  hold on;
 plot(ones(size(y_v2))*78.4,y_v2 ,'k -.','linewidth' , 0.8);
  hold on;
 plot(ones(size(y_v3))*118.6,y_v3 ,'k -.','linewidth' , 0.8);
  hold on;
 plot(ones(size(y_v4))*159.3,y_v4 ,'k -.','linewidth' , 0.8);
  hold on;
 plot(ones(size(y_v5))*196.7,y_v5 ,'k -.','linewidth' , 0.8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms Vt
i_t11=14.5;
i_t22=8.137;
i_t33=5.328;
i_t44=4.008;
i_t55=3.246;
i_t66=2.792;
V11=linspace(0,44,5751);
V22=linspace(44,78.4,5751);
V33=linspace(78.4,118.6,5751);
V44=linspace(118.6,159.3,5751);
V55=linspace(159.3,196.7,5751);
V66=linspace(196.7,228.7,5751);
M11=500;
%M22=(ones(1,5751)*M11);
%f11=(((m_t)*(1+0.04+0.0025*(i_t11)^2)));
%f22=( M.*(i_t11*eta_tr )./ R_d);
%f33=(ones(1,5751)*((m_t*g*0.02)-(0.5 * ro * C_d * A * ((Vt/3.6).^2))));
%fun1=@(Vt)(f11./(f22-f33));
%?(1800*(1+0.04+0.0025*(i^2)))*dV/(((M*i*0.93)/0.307)-((1800*9.8*0.02)+(0.5*1.2*0.31*2.18*((V/3.6)^2))))
t1=0.005*((-59689.3*atand(v_t1/((11276.2-(96.8234*i_t11)*M_1).^(0.5))))./(((11276.2-(96.8234*i_t11)*M_1).^(0.5))));
t2=0.005*((-59689.3*atand(v_t2/((11276.2-(96.8234*i_t22)*M_2).^(0.5))))./(((11276.2-(96.8234*i_t22)*M_2).^(0.5))));
t3=0.005*((-59689.3*atand(v_t3/((11276.2-(96.8234*i_t33)*M_3).^(0.5))))./(((11276.2-(96.8234*i_t33)*M_3).^(0.5))));
t4=0.005*((-59689.3*atand(v_t4/((11276.2-(96.8234*i_t44)*M_4).^(0.5))))./(((11276.2-(96.8234*i_t44)*M_4).^(0.5))));
t5=0.005*((-59689.3*atand(v_t5/((11276.2-(96.8234*i_t55)*M_5).^(0.5))))./(((11276.2-(96.8234*i_t55)*M_5).^(0.5))));
t6=0.005*((-59689.3*atand(v_t6/((11276.2-(96.8234*i_t66)*M_6).^(0.5))))./(((11276.2-(96.8234*i_t66)*M_6).^(0.5))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y_100=linspace(0,100,100);
X_100=0:0.5:10.5;
F1=((eta_tr*i_gx(1))/R_d*1000000)*M_1;
F2=((eta_tr*i_gx(2))/R_d*1000000)*M_2;
F3=((eta_tr*i_gx(3))/R_d*1000000)*M_3;
F4=((eta_tr*i_gx(4))/R_d*1000000)*M_4;
F5=((eta_tr*i_gx(5))/R_d*1000000)*M_5;
F6=((eta_tr*i_gx(6))/R_d*1000000)*M_6;
FR=(eta_tr*N_nmax*1000)/(75.000301744176809);
DV1(1:5751)=zeros;
Dt1(1:5751)=zeros;
t_1(1:5751)=zeros;
DV2(1:5751)=zeros;
Dt2(1:5751)=zeros;
t_2(1:5751)=zeros;
DV3(1:5751)=zeros;
Dt3(1:5751)=zeros;
t_3(1:5751)=zeros;
DV4(1:5751)=zeros;
Dt4(1:5751)=zeros;
t_5(1:5751)=zeros;
DV5(1:5751)=zeros;
Dt5(1:5751)=zeros;
t_5(1:5751)=zeros;
DV6(1:5751)=zeros;
Dt6(1:5751)=zeros;
t_6(1:5751)=zeros;

for i=2:5751
DV1(i)=v_t1(i)-v_t1(i-1);
Dt1(i)=((m_t*DV1(i))/(F1(i)-FR));
t_1(i)=0.8*Dt1(i)+t_1(i-1);
end
t_2(2899)=t_1(5751);
for i=2900:5751
DV2(i)=v_t2(i)-v_t2(i-1);
Dt2(i)=((m_t*DV2(i))/(F2(i)-FR));
t_2(i)=1.35*Dt2(i)+t_2(i-1);
end
t_3(3507)=t_2(5751);
for i=3508:5751
DV3(i)=v_t3(i)-v_t3(i-1);
Dt3(i)=((m_t*DV3(i))/(F3(i)-FR));
t_3(i)=1.44927* Dt3(i)+t_3(i-1);
end
t_4(4140)=t_3(5751);
for i=4141:5751
DV4(i)=v_t4(i)-v_t4(i-1);
Dt4(i)=((m_t*DV4(i))/(F4(i)-FR));
t_4(i)=1.57121*Dt4(i)+t_4(i-1);
end
t_5(4515)=t_4(5751);
for i=4516:5751
DV5(i)=v_t5(i)-v_t5(i-1);
Dt5(i)=((m_t*DV5(i))/(F5(i)-FR));
t_5(i)=1.7857*Dt5(i)+t_5(i-1);
end
t_6(4841)=t_5(5751);
for i=4842:5751
DV6(i)=v_t6(i)-v_t6(i-1);
Dt6(i)=((m_t*DV6(i))/(F6(i)-FR));
t_6(i)=1.8*Dt6(i)+t_6(i-1);
end

figure(7)
plot(t_1*1e06*0.630769230769231,v_t1*3.6* 0.755692307692308,'linewidth' , 2.5);
hold on
plot(t_2(2899:5751)*1e06*0.630769230769231,v_t2(2899:5751)*3.6* 0.755692307692308 ,'linewidth' , 2.5);
hold on
plot(t_3(3507:5751)*1e06*0.630769230769231,v_t3(3507:5751)*3.6* 0.755692307692308,'linewidth' , 2.5 );
hold on
plot(t_4(4140:5751)*1e06*0.630769230769231,v_t4(4140:5751)*3.6* 0.755692307692308 ,'linewidth' , 2.5);
hold on
plot(t_5(4515:5751)*1e06*0.630769230769231,v_t5(4515:5751)*3.6* 0.755692307692308 ,'linewidth' , 2.5);
hold on
plot(t_6(4841:5751)*1e06*0.630769230769231,v_t6(4841:5751)*3.6* 0.755692307692308 ,'linewidth' , 2.5);
hold on
plot(X_100,ones(size(X_100))*100 ,'k -.','linewidth' , 0.9);
hold on;
plot(ones(size(y_100))*10.5,y_100 ,'k -.','linewidth' , 0.8);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_n_NR=10:350;
y_vmaxR=0:6000;
x_n_N2R=10:350;
y_v1R=3000:5500;
y_v2R=3000:5500;
y_v3R=3000:5500;
y_v4R=3000:5500;
y_v5R=3000:5500;

MUe1R =  8.3 ;
SIG1R =  11.74 ;
V1R=0:0.1:16.6; 
x1R = (V1R-MUe1R)/SIG1R;
N11R=      3889*x1R +     2750;

MUe2R = 15.19;
SIG2R =  21.48;
V2R=16.6:0.1:30.38; 
x2R = (V2R-MUe2R)/SIG2R;
N22R=   3889*x2R +     2750;

MUe3R = 27.86 ;
SIG3R= 39.41;
V3R=30.38:0.1:55.728; 
x3R = (V3R-MUe3R)/SIG3R;
N33R=   3889*x3R +     2750;

MUe4R =  51.05 ;
SIG4R = 72.19;
V4R=55.728:0.1:102.09; 
x4R = (V4R-MUe4R)/SIG4R;
N44R=   3889*x4R +     2750;

MUe5R = 93.6 ;
SIG5R =132.4;
V5R=102.09:0.1:187.2; 
x5R = (V5R-MUe5R)/SIG5R;
N55R=   3889*x5R +   2750 ;

MUe6R = 172.1 ;
SIG6R =243.3;
V6R=187.2:0.1:344.1; 
x6R = (V6R-MUe6R)/SIG6R;
N66R=   3889*x6R +     2750;
%V7R=0:1:250; 
%N77R=   a*V7.^(1/b);

figure(8)
plot(V1R,N11R,'linewidth' , 2.5)
hold on;
plot(V2R,N22R,'-.','linewidth' , 2.5)   
hold on;
 plot(V3R,N33R ,'-.','linewidth' , 2.5)   
hold on;
 plot(V4R,N44R,'-.','linewidth' , 2.5 )   
hold on;
 plot(V5R,N55R,'linewidth' , 2.5) ;  
hold on;
 plot(V6R,N66R,'-.','linewidth' , 2.5);
hold on;
% plot(V7R,N77R,'linewidth' , 1.5) ;  
hold on;
 plot(x_n_NR,ones(size(x_n_NR))*5500 ,'k','linewidth' , 0.9);
hold on;
 plot(x_n_N2R,ones(size(x_n_N2R))*3000 ,'k ','linewidth' , 0.9);
hold on;
 %plot(ones(size(y_vmaxR))*212,y_vmaxR ,'k -.','linewidth' , 0.9);
 hold on;
 plot(ones(size(y_v1R))*4.6*3.6,y_v1R ,'k -.','linewidth' , 0.6);
  hold on;
 plot(ones(size(y_v2R))*8.44*3.6 ,y_v2R ,'k -.','linewidth' , 0.6);
  hold on;
 plot(ones(size(y_v3R))*3.6*15.48,y_v3R ,'k -.','linewidth' , 0.6);
  hold on;
 plot(ones(size(y_v4R))*3.6*28.36,y_v4R ,'k -.','linewidth' , 0.6);
  hold on;
 plot(ones(size(y_v5R))*52*3.6,y_v5R ,'k -.','linewidth' , 0.6);

 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)

title('\it Engine power and torque performance curves','fontsize',18,'color',[1 1 0]);
set(gca,'color',[0.75 0.75 0.75],'Xlim',[1000 7000],'Xcolor',[1 1 0],'Ylim',[0 350],'Ycolor',[1 1 0],'fontsize',11,'linewidth',2);
set(gcf,'color',[0 0 0])
xlabel('\it\bf n_[_r_p_m_]','fontsize', 11, 'color',[1 1 0]);
ylabel('\it\bf M_[_N_m_] & N_[_k_W_]','fontsize', 11, 'color',[1 1 0]);
legend( {'Power(B) ' , 'Torque(B)', 'Power(Rp)' , 'Torque(Rp)' })
grid minor  

%%%%%%%%%%%

figure(2)

title('\it Tractive effort of a combustion engine and 6-speed transmission','fontsize',18,'color',[1 1 0]);
set(gca,'color',[0.75 0.75 0.75],'Xlim',[0 250],'Xcolor',[1 1 0],'Ylim',[0 12000],'Ycolor',[1 1 0],'fontsize',11,'linewidth',2);
set(gcf,'color',[0 0 0])
xlabel('\it\bf V_[_k_m_/_h_]','fontsize', 11, 'color',[1 1 0]);
ylabel('\it\bf F_[_N_]','fontsize', 11, 'color',[1 1 0]);
legend( {' 1st gear' , ' 2nd gear' , ' 3rd gear' , ' 4th gear' ,' 5th gear' ,'6th gear' , ' alpha 0% ' , ' alpha 10% ' , ' alpha 20% ' , ' alpha 30% ' , ' alpha 40 % ' , ' alpha 50% ' ,' alpha 60%'  }  )
grid minor  
   
%%%%%%%%%

figure(3)

text(2400,210,'226')
text(3250,210,'230')
text(3650,190,'240')
text(4200,160,'260')
text(4600,110,'290')
text(5100,75,'340')
text(5500,30,'1000')
title('\it Brake specific fuel consumption map','fontsize',18,'color',[1 1 0]);
set(gca,'color',[0.75 0.75 0.75],'Xlim',[1000 6000],'Xcolor',[1 1 0],'Ylim',[0 300],'Ycolor',[1 1 0],'fontsize',11,'linewidth',2);
set(gcf,'color',[0 0 0])
xlabel('\it\bf n_[_r_p_m_]','fontsize', 11, 'color',[1 1 0]);
ylabel('\it\bf M_[_N_m_]','fontsize', 11, 'color',[1 1 0]);
legend( {'Torque[Nm]' ,'g_e[g/kWh]'}  )
grid minor

%%%%%%%%%%%%%%%%%%%%%

figure(4)
title('\it Clutch durability depending on the first gear ratio and traffic conditions','fontsize',18,'color',[1 1 0]);
set(gca,'color',[0.75 0.75 0.75],'Xlim',[11 22.5],'Xcolor',[1 1 0],'Ylim',[0 600],'Ycolor',[1 1 0],'fontsize',11,'linewidth',2);
set(gcf,'color',[0 0 0])
xlabel('\it\bf i_t_1','fontsize', 11, 'color',[1 1 0]);
ylabel('\it\bf t(clutch-max) [x1000km]','fontsize', 11, 'color',[1 1 0]);
grid minor
text(18,550,'c=50% , n_s_t_a_r_t=1500 rpm')
text(18,380,'c=75% , n_s_t_a_r_t=1500 rpm')
text(18,350,'c=50% , n_s_t_a_r_t=2000 rpm')
text(18,240,'c=75% , n_s_t_a_r_t=2000 rpm')
text(18,195,'c=50% , n_s_t_a_r_t=2500 rpm')
text(18,130,'c=75% , n_s_t_a_r_t=2500 rpm')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(5)

title('\it Three possibilities of engine speed n_v_m_a_x at maximum vehicle speed v_m_a_x','fontsize',18,'color',[1 1 0]);
set(gca,'color',[0.75 0.75 0.75],'Xlim',[50 250],'Xcolor',[1 1 0],'Ylim',[0 120],'Ycolor',[1 1 0],'fontsize',11,'linewidth',2);
set(gcf,'color',[0 0 0])
xlabel('\it\bf V_[_k_m_/_h_]','fontsize', 11, 'color',[1 1 0]);
ylabel('\it\bf N_[_k_W_]','fontsize', 11, 'color',[1 1 0]);
grid minor
legend( {'eco' ,'normal','dynamic' ,'N_r_e_s_i_s_t'} ,'location'  , 'southeast' )

%%%%%%%%%%%%%%%%%%%%

figure(6)

title('\it Method of calculating intermediate gears and overdrive','fontsize',18,'color',[1 1 0]);
set(gca,'color',[0.75 0.75 0.75],'Xlim',[0 250],'Xcolor',[1 1 0],'Ylim',[0 6500],'Ycolor',[1 1 0],'fontsize',11,'linewidth',2);
set(gcf,'color',[0 0 0])
xlabel('\it\bf V_[_k_m_/_h_]','fontsize', 11, 'color',[1 1 0]);
ylabel('\it\bf n_[_r_p_m_]','fontsize', 11, 'color',[1 1 0]);
grid minor
legend( {'1' ,'2','3' ,'4','5','6','n_x_y'} ,'location'  , 'southeast' )

%%%%%%%%%%%%%%%%%%%%%%

figure(7)

title('\it Vehicle speed vs time during acceleration.','fontsize',18,'color',[1 1 0]);
set(gca,'color',[0.75 0.75 0.75],'Xlim',[0.5 45],'Xcolor',[1 1 0],'Ylim',[0 225],'Ycolor',[1 1 0],'fontsize',11,'linewidth',2);
set(gcf,'color',[0 0 0])
xlabel('\it\bf t_[_s_]','fontsize', 11, 'color',[1 1 0]);
ylabel('\it\bf V_[_k_m_/_h_]','fontsize', 11, 'color',[1 1 0]);
grid minor
legend( {'1' ,'2','3' ,'4','5','6'} ,'location'  , 'southeast' )

%%%%%%%%%%%%%%%%%%%%%%%

figure(8)

title('\it Method of calculating intermediate gears and overdrive','fontsize',18,'color',[1 1 0]);
set(gca,'color',[0.75 0.75 0.75],'Xlim',[0 360],'Xcolor',[1 1 0],'Ylim',[0 6500],'Ycolor',[1 1 0],'fontsize',11,'linewidth',2);
set(gcf,'color',[0 0 0])
xlabel('\it\bf V_[_k_m_/_h_]','fontsize', 11, 'color',[1 1 0]);
ylabel('\it\bf n_[_r_p_m_]','fontsize', 11, 'color',[1 1 0]);
grid minor
legend( {'1' ,'2','3' ,'4','5','6'} ,'location'  , 'southeast' )

%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(9)

title('\it Tractive effort of a combustion engine and 6-speed transmission','fontsize',18,'color',[1 1 0]);
set(gca,'color',[0.75 0.75 0.75],'Xlim',[0 250],'Xcolor',[1 1 0],'Ylim',[0 12000],'Ycolor',[1 1 0],'fontsize',11,'linewidth',2);
set(gcf,'color',[0 0 0])
xlabel('\it\bf V_[_k_m_/_h_]','fontsize', 11, 'color',[1 1 0]);
ylabel('\it\bf F_[_N_]','fontsize', 11, 'color',[1 1 0]);
legend( {' 1st gear' , ' 2nd gear' , ' 3rd gear' , ' 4th gear' ,' 5th gear' ,'6th gear' , ' alpha 0% ' , ' alpha 10% ' , ' alpha 20% ' , ' alpha 30% ' , ' alpha 40 % ' , ' alpha 50% ' ,' alpha 60%'  }  )
grid minor 