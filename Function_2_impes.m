function [m] = impes2(i,j,s,dt,Bw,V,phi,C_phi,Pn,Pw,Pe,Ps,Pc,p,Swn,Sww,Swc,water_vis)
K=zeros(7,7);
K(2:6,2:6)=[.100 .050 .070 .060 .070;.080 .100 .080 .050 .090;.085 .070 .100 .070 .085;.050 .060 .070 .100 .080;.080 .070 .060 .070 .090];

Kn=(2*K(i-1,j)*K(i,j))/(K(i-1,j)+K(i,j));
Kw=(2*K(i,j-1)*K(i,j))/(K(i,j-1)+K(i,j));
Ke=(2*K(i,j+1)*K(i,j))/(K(i,j+1)+K(i,j));
Ks=(2*K(i+1,j)*K(i,j))/(K(i+1,j)+K(i,j));

Twn=(1.127*1e+2)/(water_vis*1)*Swn*Kn;

Tww=(1.127*1e+2)/(water_vis*1)*Sww*Kw;

Twe=(1.127*1e+2)/(water_vis*1)*Swc*Ke;

Tws=(1.127*1e+2)/(water_vis*1)*Swc*Ks;
%===================
qwsc=zeros(1,25);
qwsc(1)=400;
Jw=(2*pi*1.127*.090*Swc*100)/(water_vis*1*log(200/.5));
qwsc(25)=-Jw*(Pc-3000);
%===================
d_phi=phi*C_phi;
Cwp=V/(5.615*dt)*(d_phi/Bw)*Swc;
m=(Twn*(Pn-Pc)+Tww*(Pw-Pc)+Twe*(Pe-Pc)+Tws*(Ps-Pc))+qwsc(s)-Cwp*(Pc-p);
end

