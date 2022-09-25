alpha=0.0000013;%alpha*rho<>mu*eta*gamma
betta=0.4;%multiplication factor of antigens;
gamma=0.0001;%coefficient associated with the probability of neutralization of the antigen by antibodies
etta=0.2;%the number of antibodies required to neutralize a single antigen;
sigma=0.00006;%some constant, its own for each disease
tau=15;%the time during which the formulation of cascade of plasma cells;
rho=1;%the rate of antibody production by a single plasma cell
t0=0;%initial time
tn=20;%time
dt=1/24;%taking time as days
n=(tn-t0)/dt; 
T_cells=7;%life Time Plasma Cells
miu_c=1/T_cells;
T_antib = 3.5;%antibodies Decay Time
miu_f=1/T_antib;
Recovery= 20;%organ Recovery Period 20 for organism is needed 3 0 days for recoverying
miu_m=1/Recovery;
C_st=1000;
F_st=rho*C_st/miu_f;
F=zeros(1,n);%antibody concentration;
V=zeros(1,n);%concentration of pathogenic propagating antigens;
C=zeros(1,n);%plasma cell concentration;
m=zeros(1,n);%relative characteristic of target organ damage, which is [0,1] defined;
C(1)=C_st;
V(1)=1000;
F(1)=F_st;
m(1)=0;
 
for i=2:n
    C(i)=C(i-1)+dt*(sin(m(i-1)*pi/2)*alpha*V(i-1)*(dt*i-tau)*F(i-1)*(dt*i-tau)-miu_c*(C(i-1)-C_st));
    V(i)=V(i-1)+dt*((betta-gamma*F(i-1))*V(i-1));
    F(i)=F(i-1)+dt*(rho*C(i-1) - (miu_f+etta*gamma*V(i-1))*F(i-1));
    m(i)=m(i-1)+dt*(sigma*V(i-1)-miu_m*m(i-1));
    iter=i; 
    if ((m(i) >= 1 || m(i) <= 0 || V(i) <= 0 || C(i) <= 0 || F(i) <= 0))
        break
   end
end 
time=dt*(1:n);
fig=figure(); 
plot(time, C, 'b-.', 'LineWidth',2)
hold on
plot(time, V, 'r.', 'LineWidth', 2)
hold on
plot(time, F, 'g--', 'LineWidth', 2)
set(gca, 'FontSize', 14)
set(fig, 'color', 'white')
grid on
xlabel ('time[days]')
ylabel ('development scale')
legend('C(plasma-cell)', 'V(antigen)', 'F(antibody)');
fig=figure();
plot(time, m, 'b--', 'LineWidth', 2) 
set(gca, 'FontSize', 14)
set(fig, 'color', 'white')
xlabel('time[days]')
ylabel('part')
legend('m(degree of damage to organ)');
