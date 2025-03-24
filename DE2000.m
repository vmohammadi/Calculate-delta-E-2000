function DeltaE2000=DE2000(Lab_r,Lab_Pr,K)

% This function (DE2000) calculates color difference between an Lab_r (from Actual reflectance) and Lab_Pr (from XYZ predicted) 
% Input data can be single values or multiple values arranged in columns

% set the values of parametric weighting factors KL,KC,KH
% Construct a metric of the CIE Delta E 2000 recommendation, with weighting parameters kl, kc and kh as provided for in the recommendation. When not provided, these parameters default to 1.

if nargin>2
   if length(K)>2
      kL=K(1);kC=K(2);kH=K(3);
   end
else
   kL=1;kC=1;kH=1;
end

%___________________________________________________________________

L=Lab_r(:,1);a=Lab_r(:,2);b=Lab_r(:,3);
C=(a.^2+b.^2).^0.5;h=hue_angle(a,b);

Ls=Lab_Pr(:,1);as=Lab_Pr(:,2);bs=Lab_Pr(:,3);
Cs=(as.^2+bs.^2).^0.5;hs=hue_angle(as,bs);

% find G and recompute a', C' and h'
Cm=(C+Cs)/2;
G=0.5*(1-(Cm.^7./(Cm.^7+25^7)).^0.5);
a=(1+G).*a;as=(1+G).*as;
C=(a.^2+b.^2).^0.5;h=hue_angle(a,b);
Cs=(as.^2+bs.^2).^0.5;hs=hue_angle(as,bs);

% find the mean chroma and hue for each reference/sample pair
Cm=(C+Cs)/2;
hm=(h+hs)/2;
j=find(abs(h-hs)>180);hm(j)=hm(j)-180;

rad=pi/180;
Dh=angle_diff(h,hs);
DL=abs(L-Ls);
DC=abs(C-Cs);
DH=2*((C.*Cs).^0.5).*sin(rad*(Dh)/2);

% calculate T
T=1-0.17*cos(rad*(hm-30))+0.24*cos(rad*2*hm)+0.32*cos(rad*(3*hm+6))-0.2*cos(rad*(4*hm-63));

% calculate weighting factors SL, SC, SH
SL=1+(0.015.*((L+Ls)./2-50).^2)./(20+((L+Ls)./2-50).^2).^.5;
SC=1+0.045.*Cm;
SH=1+0.015.*Cm.*T;

Dt=30*exp(-(((hm-275)/25).^2));
RC=2.*((Cm.^7)./(Cm.^7+25.^7)).^.5;
RT=-sin(2*rad*Dt).*RC;

DeltaE2000=((DL./(SL.*kL)).^2+(DC./(SC.*kC)).^2+(DH./(SH.*kH)).^2+RT.*(DC./(SC.*kC)).*(DH./(SH.*kH))).^0.5;

