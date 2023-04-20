%%%%% ��������ż���ӵ�λ %%%%%  
% a0           %�����絼��
% b0           %����������
a1=81;           %��ˮ������
b1=4;           %��ˮ�絼��
% a2=           %���������
b2=0.1;           %����絼��
h=0.5;         %��ˮ���
d=0.0;         %�������
x0=0.29;         %Դ��λ�� X
y0=1.2;         %Դ��λ�� Y
z0=0.45;        %Դ��λ�� Z
Ix=-40;
Iy=0;
Iz=-4;         %��ż����
y=0:0.1:2.5;
x=0:0.02:0.5;

z=-0.28;  %���߷�Χ
%% 
for i=1:length(x)     %��
    for j=1:length(y) %��
        for ii=1:length(z)
            n=(b1-b2)/(b1+b2);
            r1m(1,i,j,ii)=((x(i)-x0)^2+(y(j)-y0)^2+(z(ii)+2*(1-1)*h-z0)^2)^(1/2);
            r2m(1,i,j,ii)=((x(i)-x0)^2+(y(j)-y0)^2+(z(ii)+2*(1-1)*h+z0)^2)^(1/2);
            r1k(1,i,j,ii)=((x(i)-x0)^2+(y(j)-y0)^2+(z(ii)-2*(1)*h+z0)^2)^(1/2);
            r2k(1,i,j,ii)=((x(i)-x0)^2+(y(j)-y0)^2+(z(ii)-2*(1)*h-z0)^2)^(1/2);
            sum_x_k(1,i,j,ii)=((n^1*Ix*(x(i)-x0)/(4*pi*b1*r1k(1,i,j,ii)^3))+(n^1*Ix*(x(i)-x0)/(4*pi*b1*r2k(1,i,j,ii)^3)));
            sum_x_m(1,i,j,ii)=((n^0*Ix*(x(i)-x0)/(4*pi*b1*r1m(1,i,j,ii)^3))+(n^0*Ix*(x(i)-x0)/(4*pi*b1*r2m(1,i,j,ii)^3)));
            sum_y_k(1,i,j,ii)=((n^1*Iy*(y(j)-y0)/(4*pi*b1*r1k(1,i,j,ii)^3))+(n^1*Iy*(y(j)-y0)/(4*pi*b1*r2k(1,i,j,ii)^3)));
            sum_y_m(1,i,j,ii)=((n^0*Iy*(y(j)-y0)/(4*pi*b1*r1m(1,i,j,ii)^3))+(n^0*Iy*(y(j)-y0)/(4*pi*b1*r2m(1,i,j,ii)^3)));
            sum_z_k(1,i,j,ii)=-((n^1*Iz*(z(ii)-2*1*h+z0)/(4*pi*b1*r1k(1,i,j,ii)^3))+(n^1*Iz*(z(ii)-2*1*h-z0)/(4*pi*b1*r2k(1,i,j,ii)^3)));
            sum_z_m(1,i,j,ii)=((n^0*Iz*(z(ii)+2*0*h-z0)/(4*pi*b1*r1m(1,i,j,ii)^3))-(n^0*Iz*(z(ii)+2*0*h+z0)/(4*pi*b1*r2m(1,i,j,ii)^3)));
            ux(1,i,j,ii)=sum_x_k(1,i,j,ii)+sum_x_m(1,i,j,ii);
            uy(1,i,j,ii)=sum_y_k(1,i,j,ii)+sum_y_m(1,i,j,ii);
            uz(1,i,j,ii)=sum_z_k(1,i,j,ii)+sum_z_m(1,i,j,ii);
              for k=2:10
                  for m=2:11
                      r1m(m,i,j,ii)=((x(i)-x0)^2+(y(j)-y0)^2+(z(ii)+2*(m-1)*h-z0)^2)^(1/2);
                      r2m(m,i,j,ii)=((x(i)-x0)^2+(y(j)-y0)^2+(z(ii)+2*(m-1)*h+z0)^2)^(1/2);
                      r1k(k,i,j,ii)=((x(i)-x0)^2+(y(j)-y0)^2+(z(ii)-2*(k)*h+z0)^2)^(1/2);
                      r2k(k,i,j,ii)=((x(i)-x0)^2+(y(j)-y0)^2+(z(ii)-2*(k)*h-z0)^2)^(1/2);
                      
                      ux_k(k,i,j,ii)=((n^k*Ix*(x(i)-x0)/(4*pi*b1*r1k(k,i,j,ii)^3))+(n^k*Ix*(x(i)-x0)/(4*pi*b1*r2k(k,i,j,ii)^3)));
                      sum_x_k(k,i,j,ii)=ux_k(k,i,j,ii)+sum_x_k(k-1,i,j,ii);
                      uy_k(k,i,j,ii)=((n^k*Iy*(y(j)-y0)/(4*pi*b1*r1k(k,i,j,ii)^3))+(n^k*Iy*(y(j)-y0)/(4*pi*b1*r2k(k,i,j,ii)^3)));
                      sum_y_k(k,i,j,ii)=uy_k(k,i,j,ii)+sum_y_k(k-1,i,j,ii);
                      uz_k(k,i,j,ii)=-((n^k*Iz*(z(ii)-2*k*h+z0)/(4*pi*b1*r1k(k,i,j,ii)^3))+(n^k*Iz*(z(ii)-2*k*h-z0)/(4*pi*b1*r2k(k,i,j,ii)^3)));
                      sum_z_k(k,i,j,ii)=uz_k(k,i,j,ii)+sum_z_k(k-1,i,j,ii);
                      
                      ux_m(m,i,j,ii)=((n^(m-1)*Ix*(x(i)-x0)/(4*pi*b1*r1m(m,i,j,ii)^3))+(n^(m-1)*Ix*(x(i)-x0)/(4*pi*b1*r2m(m,i,j,ii)^3)));
                      sum_x_m(m,i,j,ii)=ux_m(m,i,j,ii)+sum_x_m(m-1,i,j,ii);
                      uy_m(m,i,j,ii)=((n^(m-1)*Iy*(y(j)-y0)/(4*pi*b1*r1m(m,i,j,ii)^3))+(n^(m-1)*Iy*(y(j)-y0)/(4*pi*b1*r2m(m,i,j,ii)^3)));
                      sum_y_m(m,i,j,ii)=uy_m(m,i,j,ii)+sum_y_m(m-1,i,j,ii);
                      uz_m(m,i,j,ii)=((n^(m-1)*Iz*(z(ii)+2*m*h-z0)/(4*pi*b1*r1m(m,i,j,ii)^3))-(n^(m-1)*Iz*(z(ii)+2*m*h+z0)/(4*pi*b1*r2m(m,i,j,ii)^3)));
                      sum_z_m(m,i,j,ii)=uz_m(m,i,j,ii)+sum_z_m(m-1,i,j,ii);
                      
                      end_k=max(k);
                      end_m=max(m);
                      end_sum_x_k(i,j,ii)=sum_x_k(end_k,i,j,ii);
                      end_sum_x_m(i,j,ii)=sum_x_m(end_m,i,j,ii);
                      end_sum_y_k(i,j,ii)=sum_y_k(end_k,i,j,ii);
                      end_sum_y_m(i,j,ii)=sum_y_m(end_m,i,j,ii);
                      end_sum_z_k(i,j,ii)=sum_z_k(end_k,i,j,ii);
                      end_sum_z_m(i,j,ii)=sum_z_m(end_m,i,j,ii);
                  end
              end
             ux=end_sum_x_k+end_sum_x_m;
             uy=end_sum_y_k+end_sum_y_m;
             uz=end_sum_z_k+end_sum_z_m;
             u=ux+uy+uz;

             
        end
    end
end
%% ���ݼ��ܡ���ͼ�͵���
%XX=result(:,2);YY=result(:,3);UU=result(:,1);

% [XX1,YY1]=meshgrid(-48:2:48,-12:0.5:12)
% tri=delaunay(XX1,YY1); 
% trisurf(tri,XX1,YY1,u);
% shading interp

[XX1,YY1]=meshgrid(0:0.02:0.5,0:0.1:2.5);

[XX2,YY2]=meshgrid(0:0.01:0.5,0:0.05:2.5);
u2=interp2(XX1,YY1,u,XX2,YY2,'spline');
tri=delaunay(XX2,YY2); 
figure(5)
trisurf(tri,XX2,YY2,u2);
shading interp
c=[];
result=[];
