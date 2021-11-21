% Filename: Dis.m
% Author: Yifei Li
% Queensland University of Technology, Brisbane, Australia, Nov 2021
% Reference:  Y. Li, P.R. Buenzli, M.J. Simpson (2021) 
% Interpreting how nonlinear diffusion affects the fate of bistable populations using a discrete modelling framework
% The script contains:
%   - one call to the function Columndensity_Dis to generate the averaged
%     data correspondign to Figures 6(b), the averaged column density at time t.
%   - one call to the function Totaldensity_Dis to generate the averaged
%     data corresponding to Figures 6(c), the total population density.

realisations=10; %times of realisations
r=4; % the concentric ring of the spatial template has diameter r
neighbour=3*r*(r+1); %number of neighbouring sites for the spatial template where the concentric ring has diameter r
ini=0.2; %initial density, the initial width is ini*L, where L=100;
p=0.006; %Intrinsic growth rate

%Produce column density

figure
Columndensity_Dis(realisations,ini,r,neighbour,p)


%Produce total population density

figure
Totaldensity_Dis(realisations,ini,r,neighbour,p)
%}

toc
function Columndensity(realisations,ini,r,neighbour,p)
t=200; %time for generating the column density
[x,columndensity]=Obtain_columndensity(realisations,t,ini,r,neighbour,p);
plot(x,columndensity)
end

function Totaldensity(realisations,ini,r,neighbour,p)
MaxT=10000; %time for generating the total density
tstep=MaxT/100; % time step for recording data
totaldensity=Obtain_totaldensity(realisations,MaxT,tstep,ini,r,neighbour,p);
plot(totaldensity)
end

function totaldensity_final=Obtain_totaldensity(realisations,MaxT,tstep,ini,r,neighbour,p)
    totaldensity=Produce_totaldensity(MaxT,tstep,ini,r,neighbour,p);
    total=zeros(realisations,length(totaldensity));
    total(1,:)=totaldensity;
    parfor i=2:realisations
        total(i,:)=Produce_totaldensity(MaxT,tstep,ini,r,neighbour,p);
    end
    totaldensity_final=sum(total)./realisations;
end
function totaldensity=Produce_totaldensity(MaxT,tstep,ini,r,neighbour,p)
  [N,M,A]=initial(ini);
  [count,indexA]=Initialize_agentindex(A,N,M);
  totalnumber=(N*M)/2;
  totaldensity=ones(1,MaxT/tstep+1);
  coor=1;
  totaldensity(1)=count;
  for i=1:MaxT
      if count>totalnumber-1
          totaldensity(coor+1:end)=totalnumber;
          break;
      elseif count<1
          totaldensity(coor+1:end)=0;
          break;
      end
      [count,indexA,A]=iteration(A,N,M,indexA,count,r,neighbour,p);
      if i>coor*tstep-1
          coor=coor+1;
          totaldensity(coor)=count;
      end
  end
  totaldensity=totaldensity./totalnumber;
end
function [x,columndensity]=Obtain_columndensity(realisations,t,ini,r,neighbour,p)
    totaldensity=Produce_singlecolumndensity(t,ini,r,neighbour,p);
    total=zeros(realisations,length(totaldensity));
    total(1,:)=totaldensity;
    parfor i=2:realisations
        total(i,:)=Produce_singlecolumndensity(t,ini,r,neighbour,p);
    end  
    columndensity=sum(total)./realisations;
    x=(1:length(columndensity))-0.5;
end
function columndensity=Produce_singlecolumndensity(t,ini,r,neighbour,p)
  [N,M,A]=initial(ini);
  [count,indexA]=Initialize_agentindex(A,N,M);
  totalnumber=N*M/2-1;
  %iteration of movement events
  for i=1:t
      if count>totalnumber||count<10
          break
      end
      [count,indexA,A]=iteration(A,N,M,indexA,count,r,neighbour,p);
  end
  columndensity=Producedensity(A,N,M);
end
function density=Producedensity(A,N,M)
    len=M/2;
    density=zeros(1,len);
    i=1;
    while i<len+1
        density(i)=(sum(sum(A(:,i*2-1:i*2))))/N;
        i=i+1;
    end
end
function [count,indexA,A]=iteration(A,N,M,indexA,count,r,neighbour,p)
    j=1;
    randomchoose=randi(count,count,1);
    while j<count+1
            row=indexA(randomchoose(j),1);
            col=indexA(randomchoose(j),2);
            [newrow,newcol]=judge_M(A,row,col,N,M);
            if newrow>0
                A(row,col)=0;
                A(newrow,newcol)=1;
                indexA(randomchoose(j),1)=newrow;
                indexA(randomchoose(j),2)=newcol;
            end
        j=j+1;
    end
    j=1;
    count_fix=count;
    while j<count_fix+1
        if rand(1)<p
            randomchooseindex=randi(count);
            row=indexA(randomchooseindex,1);
            col=indexA(randomchooseindex,2);
            [newrow,newcol,occu]=judge_P(A,row,col,N,M,r,neighbour);
            Pii=Func_P(occu);
            if Pii>0
                if rand(1)<Pii
                    A(newrow,newcol)=1;
                    count=count+1;
                    indexA(count,1)=newrow;
                    indexA(count,2)=newcol;
                end
            elseif Pii<0
                if rand(1)<-Pii
                    A(row,col)=0;
                    indexA(randomchooseindex,:)=[];
                    count=count-1;
                end
            end
        end
        j=j+1;
    end
end
function [count,indexA]=Initialize_agentindex(A,N,M)
    count=0;
    indexA=zeros(1,2);
    for i=1:N
        for j=1:M
            if A(i,j)>0
                count=count+1;
                indexA(count,1)=i;
                indexA(count,2)=j;
            end
        end
    end
end
function [N,M,A]=initial(ini)
    N=116;%this has to be an even number so that the periodic boundary conditions works for the hexagonal lattice
    M=100*2;%this has to be an even number as well
    N_up=0;
    N_down=N;
    len=ini*M/4;
    M_right=round(M/4+len)*2;
    M_left=M-M_right;
    A=zeros(N,M);
    for i=N_up+1:N_down
        for j=M_left+1:M_right
            if mod(i+j,2)>0
                A(i,j)=1;
            end
        end
    end
end
function [newrow,newcol]=judge_M(A,row,col,N,M)
    r=1;
	neighbour=6;
    coordinate=zeros(neighbour,2);
    count=0;
    for i=row-r:row-1
        for j=col+row-2*r-i:2:col-row+2*r+i
            [newi,newj]=isoverBC(i,j,N,M);
            if A(newi,newj)<1
                count=count+1;
                coordinate(count,1)=newi;
                coordinate(count,2)=newj;
            end
        end
    end
    for j=col-2*r:2:col+2*r
        i=row;
        if j~=col
            [newi,newj]=isoverBC(i,j,N,M);
            if A(newi,newj)<1
               count=count+1;
               coordinate(count,1)=newi;
               coordinate(count,2)=newj;
            end   
        end
    end
    for i=row+1:row+r
        for j=col-row-2*r+i:2:col+row+2*r-i
            [newi,newj]=isoverBC(i,j,N,M);
            if A(newi,newj)<1
                count=count+1;
                coordinate(count,1)=newi;
                coordinate(count,2)=newj;
            end
        end
    end
    newrow=0;
    newcol=0;
    if count>0
        occu=(neighbour-count)/(neighbour);
        move=Func_M(occu);
        if move>rand(1)
            selectedsite=randi(count);
            newrow=coordinate(selectedsite,1);
            newcol=coordinate(selectedsite,2);
        end
    end
end
%movement crowding function
function a=Func_M(C)
    a=1-C; %leading to D(C)=1
    %a=(1-C)*(1+C/2); %leading to D(C)=1+C(1-C/2)
    %a=(1-C)*(1-C/2); %leading to D(C)=1-C(1-C/2)
end
function [newrow,newcol,occu]=judge_P(A,row,col,N,M,r,neighbour)
    coordinate=zeros(neighbour,2);
    count=0;
    for i=row-r:row-1
        for j=col+row-2*r-i:2:col-row+2*r+i
            [newi,newj]=isoverBC(i,j,N,M);
            if A(newi,newj)<1
                count=count+1;
                coordinate(count,1)=newi;
                coordinate(count,2)=newj;
            end
        end
    end
    for j=col-2*r:2:col+2*r
        i=row;
        if j~=col
          [newi,newj]=isoverBC(i,j,N,M);
          if A(newi,newj)<1
              count=count+1;
              coordinate(count,1)=newi;
              coordinate(count,2)=newj;
          end
        end
    end
    for i=row+1:row+r
        for j=col-row-2*r+i:2:col+row+2*r-i
            [newi,newj]=isoverBC(i,j,N,M);
            if A(newi,newj)<1
                count=count+1;
                coordinate(count,1)=newi;
                coordinate(count,2)=newj;
            end
        end
    end
    occu=(neighbour-count)/(neighbour);
    if count>0
        selectedsite=randi(count);
        newrow=coordinate(selectedsite,1);
        newcol=coordinate(selectedsite,2);
    else
        newrow=0;
        newcol=0;
    end
end
%growth crowding function
function a=Func_P(C)
    a=2.5*(1-C)*(C-0.4);
end
function [newrow,newcol]=isoverBC(newrow,newcol,N,M)
    if newrow<1
        newrow=N+newrow;
    end
    if newrow>N
        newrow=newrow-N;
    end
    if newcol<1
        newcol=M+newcol;
    end
    if newcol>M
        newcol=newcol-M;
    end
end