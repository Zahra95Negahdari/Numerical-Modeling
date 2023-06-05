function [ m ] = impes1(i,j,s,dt,Bo,Bw,V,phi,C_phi,p,Son,Swn,Sow,Sww,Soc,Swc,oil_vis,water_vis)
K=zeros(7,7);
K(2:6,2:6)=[.100 .050 .070 .060 .070;.080 .100 .080 .050 .090;.085 .070 .100 .070 .085;.050 .060 .070 .100 .080;.080 .070 .060 .070 .090];

Kn=(2*K(i-1,j)*K(i,j))/(K(i-1,j)+K(i,j));
Kw=(2*K(i,j-1)*K(i,j))/(K(i,j-1)+K(i,j));
Ke=(2*K(i,j+1)*K(i,j))/(K(i,j+1)+K(i,j));
Ks=(2*K(i+1,j)*K(i,j))/(K(i+1,j)+K(i,j));

Ton=(1.127*10^2)/(oil_vis*1.2)*Son*Kn;
Twn=(1.127*10^2)/(water_vis*1)*Swn*Kn;

Tow=(1.127*10^2)/(oil_vis*1.2)*Sow*Kw;
Tww=(1.127*10^2)/(water_vis*1)*Sww*Kw;

Toe=(1.127*10^2)/(oil_vis*1.2)*Soc*Ke;
Twe=(1.127*10^2)/(water_vis*1)*Swc*Ke;

Tos=(1.127*10^2)/(oil_vis*1.2)*Soc*Ks;
Tws=(1.127*10^2)/(water_vis*1)*Swc*Ks;
%===================
d_phi=phi*C_phi;
Cop=V/(5.615*dt)*(d_phi/Bo)*(1-Swc);
Cwp=V/(5.615*dt)*(d_phi/Bw)*Swc;
%===================
qosc=zeros(1,25);
qwsc=zeros(1,25);
Jo=(2*pi*1.127*.090*Soc*100)/(oil_vis*1.2*log(200/.5));
Jw=(2*pi*1.127*.090*Swc*100)/(water_vis*1*log(200/.5));
qosc(25)=-Jo*(p-3000);
qwsc(1)=400;
qwsc(25)=-Jw*(p-3000);
%===================
m(1,s)=Bo*Ton+Bw*Twn;
m(1,s+4)=Bo*Tow+Bw*Tww;
m(1,s+5)=-(Bo*Cop+Bw*Cwp)-(Bo*Ton+Bw*Twn+Bo*Tow+Bw*Tww+Bo*Toe+Bw*Twe+Bo*Tos+Bw*Tws);
m(1,s+6)=Bo*Toe+Bw*Twe;
m(1,s+10)=Bo*Tos+Bw*Tws;
m(1,36)=-(Bo*Cop+Bw*Cwp)*p-Bo*qosc(s)-Bw*qwsc(s);
end
