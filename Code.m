clc
clear All
close All
%% sec.1 Initial Values

format long
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
n_mreza=0:1:6500;
w=n_mreza*0.1047;
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
     
%% Sec.2 Engine Power & Torque
syms e_n e_M n_M M_M M_N
E=e_n*e_M;
E1=(n_N/n_M)*(M_M/M_N);
N_nmax = ((-0.6+k_SI)*(n_max/n_N)+(1.6-k_SI))*N_N;

n=750:6500;
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

N = ((N_M1/(n_M1- n_idle )).*(n - n_idle )).*(n>=750 & n<= 2000) + (((N_M2-N_M1)/(n_M2-n_M1)).*(n-n_M1)+ N_M1).*( n>2000 & n<=3000 )+((((M_N-M_M2)/(n_N-n_M2).^p).*((n-n_M2).^p)+ M_M2).*(n/9550)).* ( n> 3000 & n<=5500) + (((N_nmax-N_N)/(n_max-n_N).^2).*((n-n_N).^2)+ N_N).* ( ( n>5500)) ;
M = 10000*((((N_M1/(n_M1- n_idle )).*(n - n_idle ))./n).*(n>=750 & n<= 2000) + ((((N_M2-N_M1)/(n_M2-n_M1)).*(n-n_M1)+ N_M1)./n).*( n>2000 & n<=3000 )+(((((M_N-M_M2)/(n_N-n_M2).^p).*((n-n_M2).^p)+ M_M2).*(n/9550))./n).* ( n> 3000 & n<=5500) + ((((N_nmax-N_N)/(n_max-n_N).^2).*((n-n_N).^2)+ N_N)./n).* ( ( n>5500)));
M_1= 10000*((((N_M1/(n_M1- n_idle )).*((v_t1*((i_gx(1)*60)/(R_d*2*pi))) - n_idle ))./(v_t1*((i_gx(1)*60)/(R_d*2*pi)))).*(n>=750 & n<= 2000) + ((((N_M2-N_M1)/(n_M2-n_M1)).*((v_t1*((i_gx(1)*60)/(R_d*2*pi)))-n_M1)+ N_M1)./(v_t1*((i_gx(1)*60)/(R_d*2*pi)))).*( n>2000 & n<=3000 )+(((((M_N-M_M2)/(n_N-n_M2).^p).*(((v_t1*((i_gx(1)*60)/(R_d*2*pi)))-n_M2).^p)+ M_M2).*((v_t1*((i_gx(1)*60)/(R_d*2*pi)))/9550))./(v_t1*((i_gx(1)*60)/(R_d*2*pi)))).* ( n> 3000 & n<=5500) + ((((N_nmax-N_N)/(n_max-n_N).^2).*(((v_t1*((i_gx(1)*60)/(R_d*2*pi)))-n_N).^2)+ N_N)./(v_t1*((i_gx(1)*60)/(R_d*2*pi)))).* ( ( n>5500)));
M_2= 10000*((((N_M1/(n_M1- n_idle )).*((v_t2*((i_gx(2)*60)/(R_d*2*pi))) - n_idle ))./(v_t2*((i_gx(2)*60)/(R_d*2*pi)))).*(n>=750 & n<= 2000) + ((((N_M2-N_M1)/(n_M2-n_M1)).*((v_t2*((i_gx(2)*60)/(R_d*2*pi)))-n_M1)+ N_M1)./(v_t2*((i_gx(2)*60)/(R_d*2*pi)))).*( n>2000 & n<=3000 )+(((((M_N-M_M2)/(n_N-n_M2).^p).*(((v_t2*((i_gx(2)*60)/(R_d*2*pi)))-n_M2).^p)+ M_M2).*((v_t2*((i_gx(2)*60)/(R_d*2*pi)))/9550))./(v_t2*((i_gx(2)*60)/(R_d*2*pi)))).* ( n> 3000 & n<=5500) + ((((N_nmax-N_N)/(n_max-n_N).^2).*(((v_t2*((i_gx(2)*60)/(R_d*2*pi)))-n_N).^2)+ N_N)./(v_t2*((i_gx(2)*60)/(R_d*2*pi)))).* ( ( n>5500)));
M_3= 10000*((((N_M1/(n_M1- n_idle )).*((v_t3*((i_gx(3)*60)/(R_d*2*pi))) - n_idle ))./(v_t3*((i_gx(3)*60)/(R_d*2*pi)))).*(n>=750 & n<= 2000) + ((((N_M2-N_M1)/(n_M2-n_M1)).*((v_t3*((i_gx(3)*60)/(R_d*2*pi)))-n_M1)+ N_M1)./(v_t3*((i_gx(3)*60)/(R_d*2*pi)))).*( n>2000 & n<=3000 )+(((((M_N-M_M2)/(n_N-n_M2).^p).*(((v_t3*((i_gx(3)*60)/(R_d*2*pi)))-n_M2).^p)+ M_M2).*((v_t3*((i_gx(3)*60)/(R_d*2*pi)))/9550))./(v_t3*((i_gx(3)*60)/(R_d*2*pi)))).* ( n> 3000 & n<=5500) + ((((N_nmax-N_N)/(n_max-n_N).^2).*(((v_t3*((i_gx(3)*60)/(R_d*2*pi)))-n_N).^2)+ N_N)./(v_t3*((i_gx(3)*60)/(R_d*2*pi)))).* ( ( n>5500)));
M_4= 10000*((((N_M1/(n_M1- n_idle )).*((v_t4*((i_gx(4)*60)/(R_d*2*pi))) - n_idle ))./(v_t4*((i_gx(4)*60)/(R_d*2*pi)))).*(n>=750 & n<= 2000) + ((((N_M2-N_M1)/(n_M2-n_M1)).*((v_t4*((i_gx(4)*60)/(R_d*2*pi)))-n_M1)+ N_M1)./(v_t4*((i_gx(4)*60)/(R_d*2*pi)))).*( n>2000 & n<=3000 )+(((((M_N-M_M2)/(n_N-n_M2).^p).*(((v_t4*((i_gx(4)*60)/(R_d*2*pi)))-n_M2).^p)+ M_M2).*((v_t4*((i_gx(4)*60)/(R_d*2*pi)))/9550))./(v_t4*((i_gx(4)*60)/(R_d*2*pi)))).* ( n> 3000 & n<=5500) + ((((N_nmax-N_N)/(n_max-n_N).^2).*(((v_t4*((i_gx(4)*60)/(R_d*2*pi)))-n_N).^2)+ N_N)./(v_t4*((i_gx(4)*60)/(R_d*2*pi)))).* ( ( n>5500)));
M_5= 10000*((((N_M1/(n_M1- n_idle )).*((v_t5*((i_gx(5)*60)/(R_d*2*pi))) - n_idle ))./(v_t5*((i_gx(5)*60)/(R_d*2*pi)))).*(n>=750 & n<= 2000) + ((((N_M2-N_M1)/(n_M2-n_M1)).*((v_t5*((i_gx(5)*60)/(R_d*2*pi)))-n_M1)+ N_M1)./(v_t5*((i_gx(5)*60)/(R_d*2*pi)))).*( n>2000 & n<=3000 )+(((((M_N-M_M2)/(n_N-n_M2).^p).*(((v_t5*((i_gx(5)*60)/(R_d*2*pi)))-n_M2).^p)+ M_M2).*((v_t5*((i_gx(5)*60)/(R_d*2*pi)))/9550))./(v_t5*((i_gx(5)*60)/(R_d*2*pi)))).* ( n> 3000 & n<=5500) + ((((N_nmax-N_N)/(n_max-n_N).^2).*(((v_t5*((i_gx(5)*60)/(R_d*2*pi)))-n_N).^2)+ N_N)./(v_t5*((i_gx(5)*60)/(R_d*2*pi)))).* ( ( n>5500)));
M_6= 10000*((((N_M1/(n_M1- n_idle )).*((v_t6*((i_gx(6)*60)/(R_d*2*pi))) - n_idle ))./(v_t6*((i_gx(6)*60)/(R_d*2*pi)))).*(n>=750 & n<= 2000) + ((((N_M2-N_M1)/(n_M2-n_M1)).*((v_t6*((i_gx(6)*60)/(R_d*2*pi)))-n_M1)+ N_M1)./(v_t6*((i_gx(6)*60)/(R_d*2*pi)))).*( n>2000 & n<=3000 )+(((((M_N-M_M2)/(n_N-n_M2).^p).*(((v_t6*((i_gx(6)*60)/(R_d*2*pi)))-n_M2).^p)+ M_M2).*((v_t6*((i_gx(6)*60)/(R_d*2*pi)))/9550))./(v_t6*((i_gx(6)*60)/(R_d*2*pi)))).* ( n> 3000 & n<=5500) + ((((N_nmax-N_N)/(n_max-n_N).^2).*(((v_t6*((i_gx(6)*60)/(R_d*2*pi)))-n_N).^2)+ N_N)./(v_t6*((i_gx(6)*60)/(R_d*2*pi)))).* ( ( n>5500)));
N_mreza=(((N_N/n_N1)*w)+((N_N/(n_N1^2))*(w.^2))+(((-N_N)/((n_N1)^3))*(w.^3)));
M_mreza=1000*(((N_N/n_N1))+((N_N/(n_N1^2))*(w))+(((-N_N)/((n_N1)^3))*(w.^2)));

figure(1)
   plot (n_mreza , N_mreza ,'--', 'linewidth' , 2.5)
   hold on;
   plot ( n_mreza , M_mreza,'--', 'linewidth' , 2.5)
   hold on;  
   plot(n , N, 'linewidth' , 2.5);
   hold on;
   plot(n , M , 'linewidth' , 2.5);
   hold on;
  
%% Sec.3 Tractive Forces
for i=1:6
syms   teta


%%%%%%%%%%%%%%%%%%%%%%F_D
v=n_mreza * 2 * pi * R_d * 3.6 /( i_gx(i)* 60);
F_D1= M_mreza*i_gx(i)* eta_tr / R_d;
F_D2=m_t*g*Mu*w_D;
figure(6)


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



%%%%%%%%%%%%%%%%%%%%%F_slope

F_slope= m_t * g * sind(alpha);
%%%%%%%%%%%%%%%%%%%F_tot
F_tot=(F_aero+F_roll+F_slope);
figure(2)
       plot(V , F_tot,'--','linewidth' , 1.5);
       
hold on;
N_tot=F_tot .* V;
figure(6)
       plot(V , F_tot,'--','linewidth' , 1.5);
       
hold on;
end
%%%%%%%%%%%%%%%%%

%% Sec.4 Specific Fuel Consumption Map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%% Sec.5 Clutch Durability

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

%% Sec.6 Power

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sec.7  Figures

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

title('\it Three possibilities of engine Power n_v_m_a_x at maximum vehicle speed v_m_a_x','fontsize',18,'color',[1 1 0]);
set(gca,'color',[0.75 0.75 0.75],'Xlim',[50 250],'Xcolor',[1 1 0],'Ylim',[0 120],'Ycolor',[1 1 0],'fontsize',11,'linewidth',2);
set(gcf,'color',[0 0 0])
xlabel('\it\bf V_[_k_m_/_h_]','fontsize', 11, 'color',[1 1 0]);
ylabel('\it\bf N_[_k_W_]','fontsize', 11, 'color',[1 1 0]);
grid minor
legend( {'eco' ,'normal','dynamic' ,'N_r_e_s_i_s_t'} ,'location'  , 'southeast' )

%%%%%%%%%%%%%%%%%%%%
figure(6)

title('\it Tractive effort of a combustion engine and 6-speed transmission','fontsize',18,'color',[1 1 0]);
set(gca,'color',[0.75 0.75 0.75],'Xlim',[0 250],'Xcolor',[1 1 0],'Ylim',[0 12000],'Ycolor',[1 1 0],'fontsize',11,'linewidth',2);
set(gcf,'color',[0 0 0])
xlabel('\it\bf V_[_k_m_/_h_]','fontsize', 11, 'color',[1 1 0]);
ylabel('\it\bf F_[_N_]','fontsize', 11, 'color',[1 1 0]);
legend( {' 1st gear' , ' 2nd gear' , ' 3rd gear' , ' 4th gear' ,' 5th gear' ,'6th gear' , ' alpha 0% ' , ' alpha 10% ' , ' alpha 20% ' , ' alpha 30% ' , ' alpha 40 % ' , ' alpha 50% ' ,' alpha 60%'  }  )
grid minor 