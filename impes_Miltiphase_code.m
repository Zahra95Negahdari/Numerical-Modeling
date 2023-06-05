clear
clc
close all
format short
p=zeros(7,7,1);
p(2:6,2:6,1)=5000;
Sw=zeros(7,7,1);
So=Sw;
So(2:6,2:6,1)=1;
dt=10;
n=1;
water_cut(1)=0;
oil_vis=5;
water_vis=1;
Bo=1.2;
Bw=1;
V=1000*1000*100;
phi=.2;
C_phi=4*10^-6;
c_oil(1)=0;
c_water(1)=0;
c_qinj(1)=400;
qinj(1)=400;
tic
while water_cut(n)<=.5
    n=n+1;
    qinj(n)=400;
    s=0;
    m=zeros(1,36);
    for i=2:6
        for j=2:6
            s=s+1;
            m(s,:)=impes1(i,j,s,dt,Bo,Bw,V,phi,C_phi,p(i,j,n-1),So(i-1,j,n-1),Sw(i-1,j,n-1),So(i,j-1,n-1),Sw(i,j-1,n-1),So(i,j,n-1),Sw(i,j,n-1),oil_vis,water_vis);
        end
    end
    A=m(:,6:30);
    B=m(:,36);
    p(2:6,2:6,n)=(reshape(A\B,[5,5]))';
    s=0;
    for i=2:6
        for j=2:6
            s=s+1;
            Sw(i,j,n)=Sw(i,j,n-1)+(5.615*dt*Bw)/(V*phi*(1+C_phi*(p(i,j,n)-5000)))*impes2(i,j,s,dt,Bw,V,phi,C_phi,p(i-1,j,n),p(i,j-1,n),p(i,j+1,n),p(i+1,j,n),p(i,j,n),p(i,j,n-1),Sw(i-1,j,n-1),Sw(i,j-1,n-1),Sw(i,j,n-1),water_vis);
            So(i,j,n)=1-Sw(i,j,n);
        end
    end
    Jo=(2*pi*1.127*.090*So(i,j,n)*100)/(oil_vis*1.2*log(200/.5));
    Jw=(2*pi*1.127*.090*Sw(i,j,n)*100)/(water_vis*1*log(200/.5));
    qosc(n)=(Jo*(p(i,j,n)-3000));
    qwsc(n)=(Jw*(p(i,j,n)-3000));
    water_cut(n)=qwsc(n)/(qosc(n)+qwsc(n));
    c_oil(n)=c_oil(n-1)+qosc(n);
    c_water(n)=c_water(n-1)+qwsc(n);
    c_qinj(n)=c_qinj(n-1)+400;
end
toc
tic
figure(1)
plot(qosc)
xlabel('Time steps')
ylabel('Oil production rate')
grid on
%==================
figure(2)
plot(qwsc)
xlabel('Time steps')
ylabel('Water production rate')
grid on
%=================
figure(3)
plot(qinj)
xlabel('Time steps')
ylabel('Water injection rate')
grid on
%================
figure(4)
plot(c_oil)
xlabel('Time steps')
ylabel('Cumulative oil production')
grid on
%===================
figure(5)
plot(c_water)
xlabel('Time steps')
ylabel('Cumulative water production')
grid on
%===================
figure(6)
plot(c_qinj)
xlabel('Time steps')
ylabel('Cumulative water injection')
grid on
%===================
figure(7)
plot(water_cut)
xlabel('Time steps')
ylabel('Water cut')
grid on
%===================
figure(8)
surf(p(2:6,2:6,1))
ylabel('X Axis')
xlabel('Y Axis')
zlabel('Pressure distribution')
%===================
figure(9)
surf(p(2:6,2:6,2))
ylabel('X Axis')
xlabel('Y Axis')
zlabel('Pressure distribution')
%===================
figure(10)
surf(p(2:6,2:6,100))
ylabel('X Axis')
xlabel('Y Axis')
zlabel('Pressure distribution')
%===================
figure(11)
surf(p(2:6,2:6,1000))
ylabel('X Axis')
xlabel('Y Axis')
zlabel('Pressure distribution')
%===================
figure(12)
surf(Sw(2:6,2:6,1))
ylabel('X Axis')
xlabel('Y Axis')
zlabel('Water saturation distribution')
%===================
figure(13)
surf(Sw(2:6,2:6,2))
ylabel('X Axis')
xlabel('Y Axis')
zlabel('Water saturation distribution')
%===================
figure(14)
surf(Sw(2:6,2:6,100))
ylabel('X Axis')
xlabel('Y Axis')
zlabel('Water saturation distribution')
%===================
figure(15)
surf(Sw(2:6,2:6,1000))
ylabel('X Axis')
xlabel('Y Axis')
zlabel('Water saturation distribution')
%====================
figure(16)
surf(Sw(2:6,2:6,n))
ylabel('X Axis')
xlabel('Y Axis')
zlabel('Water saturation distribution')
%====================
figure(17)
surf(So(2:6,2:6,1))
ylabel('X Axis')
xlabel('Y Axis')
zlabel('Oil saturation distribution')
%===================
figure(18)
surf(So(2:6,2:6,2))
ylabel('X Axis')
xlabel('Y Axis')
zlabel('Oil saturation distribution')
%===================
figure(19)
surf(So(2:6,2:6,100))
ylabel('X Axis')
xlabel('Y Axis')
zlabel('Oil saturation distribution')
%===================
figure(20)
surf(So(2:6,2:6,n))
ylabel('X Axis')
xlabel('Y Axis')
zlabel('Oil saturation distribution')
%=================
figure(21)
surf(Sw(2:6,2:6,1))
ylabel('X Axis')
xlabel('Y Axis')
zlabel('Krw')
%===================
figure(22)
surf(Sw(2:6,2:6,2))
ylabel('X Axis')
xlabel('Y Axis')
zlabel('Krw')
%===================
figure(23)
surf(Sw(2:6,2:6,100))
ylabel('X Axis')
xlabel('Y Axis')
zlabel('Krw')
%===================
figure(24)
surf(Sw(2:6,2:6,1000))
ylabel('X Axis')
xlabel('Y Axis')
zlabel('krw')

toc
