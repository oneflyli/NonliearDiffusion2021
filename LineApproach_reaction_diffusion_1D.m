% Filename: LineApproach_reaction_diffusion_1D.m
% Author: Yifei Li
% Queensland University of Technology, Brisbane, Australia, Jan 2021
% Reference:  Y. Li, P.R. Buenzli, M.J. Simpson (2021) 
% The script contains the function 'LineApproach_reaction_diffusion_1D', which applies
% the method of lines for solving a 1-dimensional RDE

function dudt = LineApproach_reaction_diffusion_1D(t,u,N,D,dx,p,A)
u=reshape(u,N,1);
dudt1=zeros(N,1);
for i=2:N-1
        dudt1(i)=D*((Func_Du(u(i+1))+Func_Du(u(i)))*(u(i+1)-u(i))-(Func_Du(u(i))+Func_Du(u(i-1)))*(u(i)-u(i-1)))/(2*dx^2)+p*Func_Ru(u(i),A);
end

i=1;
dudt1(i)=D*((Func_Du(u(i+1))+Func_Du(u(i)))*(u(i+1)-u(i))-(Func_Du(u(i))+Func_Du(u(N)))*(u(i)-u(N)))/(2*dx^2)+p*Func_Ru(u(i),A);

i=N;
dudt1(i)=D*((Func_Du(u(1))+Func_Du(u(i)))*(u(1)-u(i))-(Func_Du(u(i))+Func_Du(u(i-1)))*(u(i)-u(i-1)))/(2*dx^2)+p*Func_Ru(u(i),A);

dudt=dudt1;
end
%Diffusion term
function a=Func_Du(C)
    a=1; %relate to G(C)=1-C
    %a=-C^2/2+C+1; %relate to G(C)=(1-C)(1+C/2)
    %a=C^2/2-C+1; %relate to G(C)=(1-C)(1-C/2)
    %a=C;
    %a=C^2;
    %a=C^3;
end
%Bistable source term
function a=Func_Ru(C,A)
    a=2.5*C*(1-C)*(C-A);
end