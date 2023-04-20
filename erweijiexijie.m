 %%%%%%%% 参数定义 %%%%%%%%%%
h0=50; %埋深
r0=1; %异常体半径
p1=4; %围岩电阻率
p2=0.04; %异常体电阻率
a=[0]; %角度范围
% a=[0,60,90]
% uu=0.001;
x=-100:1:100;
% M=2*p1/(p1+2*p2)*r0^2*uu;
figure(1)
M=118
Ixdl=30
 for i=1:length(a)
    for j=1:length(x)
        U(i,j)=M*(x(j)*cosd(a(i))-h0*sind(a(i)))/((h0^2+x(j)^2)^(3/2));
        UU(j)=Ixdl*x(j)/(4*pi*p2*(x(j)^2+h0^2)^(3/2))+(p1-p2)/(p1+p2)*Ixdl*x(j)/(4*pi*p2*(x(j)^2+h0^2)^(3/2));
        err(j)=(UU(j)-U(i,j))/U(i,j);
    end
 end
%U_3D=u2(:,101);
%U_3D=flipud(U_3D);
%plot(-x,U,'-r',-x,U_3D,'-b');
plot(-x,U,'-r',-x,UU,'--b');

figure(2)

plot(-x,err,'-y');
U2=U';
U3=UU';
err_2=err'*10;

