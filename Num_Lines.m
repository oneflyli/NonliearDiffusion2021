% Filename: Num_Lines.m
% Author: Yifei Li
% Queensland University of Technology, Brisbane, Australia, Nov 2021
% Reference:  Y. Li, P.R. Buenzli, M.J. Simpson (2021) 
% Interpreting how nonlinear diffusion affects the fate of bistable populations using a discrete modelling framework
% The script contains:
%   - one call to the function Columndensity_Num to generate 
%     Figures 6(b), the averaged column density at time t in the continuum model.
%   - one call to the function Totaldensity_Num to generate 
%     Figures 6(c), the total population density in the continuum model.

D=1/4; %diffusivity constant
ini=0.2; %initial width is ini*L, where L=100
p=0.006; %P
A=0.4;%Allee threshold

%Obtain column density 
figure
Columndensity_Num(ini,D,p,A)

%obtain total density
figure
Totaldensity_Num(ini,D,p,A);


function Columndensity_Num(ini,D,p,A)
t1=200;
[x,u1]=Export_columndensity1D_Num(ini,t1,D,p,A);
plot(x,u1,'b')
end

function Totaldensity_Num(ini,D,p,A)
t=10000;
total=totaldensity1D(ini,t,D,p,A);
plot(total(:,1),total(:,2))
end

function [x,column]=Export_columndensity1D_Num(ini,T,D,p,A)   
    L = 100;
    dx = 0.5;
    N=L/dx;
    u_initial = 1;
    u0 = zeros(N+1,1);
    len=ini*N/2;
    right=round(N/2+len);
    left=N-right;
    for i=left+1:right
        u0(i)=u_initial;
    end
    u0=reshape(u0,[],1);
    tspan = 0:T/2:T;
    [t,u] = ode45(@(t,u) LineApproach_reaction_diffusion_1D(t,u,N+1,D,dx,p,A), tspan, u0);
    u = reshape(u, [], N+1,1);
    column=u(end,:);
    x=(0:length(column)-1)*dx;
end

function total=totaldensity1D(ini,T,D,p,A)
    tstep=100;
    L = 100;
    dx = 0.5;
    N=L/dx;
    u0 = zeros(N+1,1);
    u_initial=1;
    len=ini*N/2;
    right=round(N/2+len);
    left=N-right;
    for i=left+1:right
        u0(i)=u_initial;
    end
    u0=reshape(u0,[],1);
    tspan = 0:round(T/tstep):T;
    [t,u] = ode45(@(t,u) LineApproach_reaction_diffusion_1D(t,u,N+1,D,dx,p,A), tspan, u0);
    u = reshape(u, [], N+1,1);
    total=zeros(1,length(tspan));
    for i=1:length(tspan)
        total(i)=sum(sum(u(i,:)))/(N+1);
    end
    a=0:length(total)-1;
    total=[(a/10)',total'];
end
