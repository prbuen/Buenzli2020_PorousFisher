% Nondimensionalised Porous-Fisher model applied to tissue engineering pore infilling.
%
% Authors: Pascal R Buenzli, Matthew J Simpson, 2019-2020

function dudt = odefn(t,u,N,alpha2,dx)
u=reshape(u,N,N);
dudt1=zeros(N,N);

for i=2:N-1
    for j=2:N-1
        dudt1(i,j)=alpha2*((u(i+1,j)+u(i,j))*(u(i+1,j)-u(i,j))-(u(i,j)+u(i-1,j))*(u(i,j)-u(i-1,j)))/(2.*dx^2) ...
            +alpha2*((u(i,j+1)+u(i,j))*(u(i,j+1)-u(i,j))-(u(i,j)+u(i,j-1))*(u(i,j)-u(i,j-1)))/(2.*dx^2) ...
            +u(i,j)*(1.0-u(i,j));
    end
end

i=1;
for j=2:N-1
    dudt1(i,j)=alpha2*((u(i,j)+u(i+1,j))*(u(i+1,j)-u(i,j))-(u(i,j)+u(i+1,j))*(u(i+1,j)-u(i,j)))/(2.*dx^2) ...
        +alpha2*((u(i,j+1)+u(i,j))*(u(i,j+1)-u(i,j))-(u(i,j)+u(i,j-1))*(u(i,j)-u(i,j-1)))/(2.*dx^2) ...
        +u(i,j)*(1.0-u(i,j));
end

i=N;
for j=2:N-1
    dudt1(i,j)=alpha2*((u(i,j)+u(i-1,j))*(u(i,j)-u(i-1,j))-(u(i,j)+u(i-1,j))*(u(i,j)-u(i-1,j)))/(2.*dx^2) ...
        +alpha2*((u(i,j+1)+u(i,j))*(u(i,j+1)-u(i,j))-(u(i,j)+u(i,j-1))*(u(i,j)-u(i,j-1)))/(2.*dx^2) ...
        +u(i,j)*(1.0-u(i,j));
end

j=1;
for i=2:N-1
    dudt1(i,j)=alpha2*((u(i+1,j)+u(i,j))*(u(i+1,j)-u(i,j))-(u(i,j)+u(i-1,j))*(u(i,j)-u(i-1,j)))/(2.*dx^2) ...
        +alpha2*((u(i,j)+u(i,j+1))*(u(i,j+1)-u(i,j))-(u(i,j)+u(i,j+1))*(u(i,j+1)-u(i,j)))/(2.*dx^2) ...
        +u(i,j)*(1.0-u(i,j));
end

j=N;
for i=2:N-1
    dudt1(i,j)=alpha2*((u(i+1,j)+u(i,j))*(u(i+1,j)-u(i,j))-(u(i,j)+u(i-1,j))*(u(i,j)-u(i-1,j)))/(2.*dx^2) ...
        +alpha2*((u(i,j)+u(i,j-1))*(u(i,j)-u(i,j-1))-(u(i,j)+u(i,j-1))*(u(i,j)-u(i,j-1)))/(2.*dx^2) ...
        +u(i,j)*(1.0-u(i,j));
end

i=1;
j=1;
dudt1(i,j)=alpha2*((u(i,j)+u(i+1,j))*(u(i+1,j)-u(i,j))-(u(i,j)+u(i+1,j))*(u(i+1,j)-u(i,j)))/(2.*dx^2) ...
    +alpha2*((u(i,j)+u(i,j+1))*(u(i,j+1)-u(i,j))-(u(i,j)+u(i,j+1))*(u(i,j+1)-u(i,j)))/(2.*dx^2) ...
    +u(i,j)*(1.0-u(i,j)); 
i=N;
j=1;
dudt1(i,j)=alpha2*((u(i,j)+u(i-1,j))*(u(i,j)-u(i-1,j))-(u(i,j)+u(i-1,j))*(u(i,j)-u(i-1,j)))/(2.*dx^2) ...
    +alpha2*((u(i,j)+u(i,j+1))*(u(i,j+1)-u(i,j))-(u(i,j)+u(i,j+1))*(u(i,j+1)-u(i,j)))/(2.*dx^2) ...
    +u(i,j)*(1.0-u(i,j));
i=1;
j=N;
dudt1(i,j)=alpha2*((u(i,j)+u(i+1,j))*(u(i+1,j)-u(i,j))-(u(i,j)+u(i+1,j))*(u(i+1,j)-u(i,j)))/(2.*dx^2) ...
    +alpha2*((u(i,j)+u(i,j-1))*(u(i,j)-u(i,j-1))-(u(i,j)+u(i,j-1))*(u(i,j)-u(i,j-1)))/(2.*dx^2) ...
    +u(i,j)*(1.0-u(i,j));
i=N;
j=N;
dudt1(i,j)=alpha2*((u(i,j)+u(i-1,j))*(u(i,j)-u(i-1,j))-(u(i,j)+u(i-1,j))*(u(i,j)-u(i-1,j)))/(2.*dx^2) ...
    +alpha2*((u(i,j)+u(i,j-1))*(u(i,j)-u(i,j-1))-(u(i,j)+u(i,j-1))*(u(i,j)-u(i,j-1)))/(2.*dx^2) ...
    +u(i,j)*(1.0-u(i,j)); 
dudt=reshape(dudt1,[],1);
