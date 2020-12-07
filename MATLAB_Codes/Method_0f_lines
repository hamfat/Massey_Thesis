%We use the method of lines to solve the PDE, this code generates different spatiotemporal patterns depending on parameter values

%model parameters
global d glbar gkbar gcabar vlbar vkbar vcabar v1 v2 v3 v4 psi a b 
d=.0001; %diffusion coefficient
v1=-0.245;
v2=0.3125;
v3=-0.1405;
v4=0.1812;
vkbar=-1.125;
vcabar=1.0;
vlbar=-0.875;
gkbar=1.0;
glbar=0.25;
gcabar=0.4997;
psi=0.1665;
A_0=.3; 

% Find steady state v_0 and  N_0 using fsolve
VN= fsolve(@(X) dimless(X),[-0.7654;0.001]);
V_0 = VN(1);
N_0 = VN(2);

%simulation time
tf=500;
Tf =1000;
tspan=linspace(0,tf,Tf);
%M is the number of discretisation
M =1000;
sigma=0.1;
a=-3.0;
b=3.0 ;
x=linspace(a,b,M);
dx=(b-a)/(M-1);

%%%%%%%%%%%%%%% initial conditions%%%%%%%%%%%%%%%%%%%
v_initial=V_0*ones(size(x))+A_0*exp(-((x)/sigma).^2).*ones(size(x)); 
n_initial=N_0*ones(size(x));

initial=[v_initial, n_initial];

%%%%%%%% to obtain the solutions for the state variables V and N %%%%%%%%
[tsol,xsol]=ode15s(@(t,x) myworkinternal(t,x,M),tspan,initial);

V=xsol(:,1:M)';
N=xsol(:,M+1:end)';
%%%%%%%%%%%%%%% Figures %%%%%%%%%%%%%%%%%%%%%%%%

[T,X] = meshgrid(tsol,x);

figure(1)
 %subplot(1,2,1);
mesh(T,X,V,'FaceLighting','gouraud','LineWidth',0.5)
xlabel('\itT\rm')
 ylabel('\itX\rm')
zlabel('\itV\rm')
set ( gca, 'xdir', 'reverse' )
colorbar
subplot(1,2,2);
mesh(T,X,N,'FaceLighting','gouraud','LineWidth',0.5)
xlabel('\itT\rm')
ylabel('\itX\rm')
zlabel('\itN\rm')
set ( gca, 'xdir', 'reverse' )

f=figure(3);
imagesc(x,tsol,V');
%title(sprintf('$v_1$=%4.4f', v1),'interpreter', 'latex')
title(['v_1=',num2str(v1)])
ylabel('T','Fontsize', 20, 'interpreter', 'latex')
xlabel('X','Fontsize', 20, 'interpreter', 'latex')
colorbar
f.CurrentAxes.YDir='normal';
colormap jet

figure(4)
plot(x,V(:,550))
ylabel('V','Fontsize', 20, 'LineWidth', 1.5, 'interpreter', 'latex')
xlabel('X','Fontsize', 20, 'LineWidth', 1.5, 'interpreter', 'latex')

%%%%%%%%%% Discretization in space %%%%%%%
function xdot = myworkinternal(~,x,M)
global d glbar gkbar gcabar vlbar vkbar v1 v2 v3 v4  psi a b vcabar
V = x(1:M);
N = x(M+1:end);
dx=(b-a)/(M-1);
mu=d/dx^2;
dVdt = zeros(M,1);
dNdt = zeros(M,1);
for i=1:M
    minf=0.5*(1+tanh((V(i)-v1)/v2));
    ninf=0.5*(1+tanh((V(i)-v3)/v4));
    lambda_v=cosh((V(i)-v3)/(2*v4));
    if i==1
			dVdt(i)= 2*mu*(V(i+1)-V(i))-(glbar*(V(i)-vlbar)+gkbar*N(i)*(V(i)-vkbar)+...
			...gcabar*minf*(V(i)-vcabar));
			dNdt(i)= psi*lambda_v*(ninf-N(i));
    end
    if i==M
			dVdt(i)= 2*mu*(V(i-1)-V(i))-(glbar*(V(i)-vlbar)+gkbar*N(i)*(V(i)-vkbar)+...
			...gcabar*minf*(V(i)-vcabar));

			dNdt(i)= psi*lambda_v*(ninf-N(i));
    end
    if i>1 && i<M
			dVdt(i)= mu*(V(i+1)-2*V(i)+V(i-1))-(glbar*(V(i)-vlbar)+gkbar*N(i)*(V(i)-vkbar)+...
			...gcabar*minf*(V(i)-vcabar));
			dNdt(i)= psi*lambda_v*(ninf-N(i));
    end
end
xdot = [dVdt;dNdt];
end
%%%%%%%%%%% functions used to obtain the steady state V0 and N0 above %%%%%
function ica = I_Ca(V)
global  gcabar vcabar
    ica=gcabar * m_inf(V)*(V - vcabar);
end
function ik = I_K(V,N)
global gkbar vkbar 
    ik=gkbar  * N * (V - vkbar);
end
function ileak = I_Leak(V)
global glbar vlbar 
    ileak=glbar* (V - vlbar);
end
function VN = dimless(X)
global  psi 
    V=X(1);
    N=X(2);
    VN(1) = - I_Leak(V) - I_K(V,N) - I_Ca(V);
    VN(2) = psi*lamda(V)*(n_inf(V) - N) ;
end
function out = m_inf(V)
global  v1 v2 
	out = 0.5*(1+tanh((V-v1)/v2));
end
function out =n_inf(V)
global  v3 v4
	out = 0.5*(1+tanh((V-v3)/v4));
end
function out =lamda(V)
global  v3 v4 
	out = cosh((V-v3)/(2*v4));
end
