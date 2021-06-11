clc;clear

%position

% o is the origin
r2=6*10^-2; %m

theta2 =    45   ;

%oa= r2*exp(1i*theta2);
xa=r2*cosd(theta2);
ya=r2* sind(theta2);
r3=24.5*10^-2; %m
%ob= oa + ab
%xb+iyb= xa+i*ya + r3*exp(i*theta3)


yb=   0   ;


xb=xa+ sqrt((r3^2)-((yb-ya)^2));    %m
cosdtheta3=(xb-xa)/r3;
sindtheta3=(yb-ya)/r3;
if cosdtheta3 > 0;  %rob3 2wl
    if sindtheta3 >0;
        theta3=acosd(cosdtheta3);
    end
 
elseif sindtheta3 > 0    %rob3 tany
    if cosdtheta3 < 0
        theta3 =180- asind(sindtheta3);
    end
elseif sindtheta3 <0    %rob3 talt
    if cosdtheta3 < 0
        theta3=180+asind(-1*sindtheta3);
    end
end
if cosdtheta3 > 0 ;   %rob3 rab3
    if sindtheta3 < 0 ;  
        theta3= 360 - (acosd(cosdtheta3));
    end
end
if sindtheta3 == 0 ;
    theta3=asind(sindtheta3);
end
if cosdtheta3 == 0 ;
    theta3=acosd(cosdtheta3);
end 


%velocity



%omega2= 50 rbm
omega2=50*2*pi / 60 ; %rad/s
%va= vax+i*vay= i*omega2*r2*exp(i*theta2)
vax= -1*omega2*r2*sind(theta2);
vay= omega2*r2*cosd(theta2);
%vb*exp(i*alpha)= i*omega2*r2*exp(i*theta2)+i*omega3*r3*exp(i*theta3)  
%vb =(vax+ivay)*exp(-1*i*alpha)+ i*omega3*r3*exp(i*(theta3-alpha))
alpha=0;
%img
%0= -1*vax*sind(alpha)+vay*cod(alpha)+ omega3*r3*cosd((theta3-alpha))
omega3= (vax*sind(alpha)-(vay*cosd(alpha)))/ r3*cosd((theta3-alpha)) ;   %rad/s
%real
vb= vax*cosd(alpha)+ vay*sind(alpha)- omega3*r3*sind((theta3-alpha)); %m/s




% acceleration

%aa=aax+iaay=(-1*omega2^2*r2+i*alpha2*r2)*exp(i*theta2)


alpha2=   0   ;  % angular acceleration for cranck


aax= (-1*omega2^2*r2*cosd(theta2))- (alpha2*r2*sind(theta2));
aay=-1*omega2^2*r2*sind(theta2)+ (alpha2*r2*cosd(theta2));
%Ab = aa + Aba
%Ab*exp(i*alpha) = aax+iaay +(-1*omega^2*r3+i*alpha2*r3)*exp(i*theta3)
%img
%0=-1*aax*sind(alpha)+aay*cosd(alpha)-omega3^2*r3*sind(theta3-alpha)+alpha3*r3*cosd(theta3-alpha)
alpha3=(aax*sind(alpha)-aay*cosd(alpha)+omega3^2*r3*sind(theta3-alpha))/(r3*cosd(theta3-alpha));  %rad/s^2
%real
ab=(aax*cosd(alpha)+aay*sind(alpha))-(omega3^2*r3*cosd(theta3-alpha))-(alpha3*r3*sind(theta3-alpha)); %m/s^2


%force


m2=0.032 ; %kg
m3=0.094; %kg
m4=0.6; %kg slider
g2=r2/2; %m
g3=r3/2; %m
m12=2.452; %n.m
i2=(1/3)*m2*(r2/2)^2; %kg.m^2
i3=(1/3)*m3*(r3/2)^2; %kg.m^2
f4=m4*ab;

%ag2=ag2x+ iag2y = (-omega2^2*g2 +i alpha2*g2)*exp(i*theta2);
ag2x= -1*omega2^2*g2*cosd(theta2)- alpha2*g2*sind(theta2);
ag2y= -1*omega2^2*g2*sind(theta2)- alpha2*g2*cosd(theta2);
%ag3=aa+ag3a
%ag3x+i ag3y= aax+iaay + (-omega3^2*g3+i alpha3*g3)*exp(i*theta3);
ag3x=aax + (-1*omega3^2*g3*cosd(theta3))- alpha3*g3*sind(theta3);
ag3y= aay + (-1*omega3^2*g3*sind(theta3))+ (alpha3*g3*cosd(theta3));
x2=-1*m2*ag2x;
y2=-1*m2*ag2y;
t2=-1*i2*alpha2;
x3=-1*m3*ag3x;
y3=-1*m3*ag3y;
t3=-1*i3*alpha3;

%link 3 & slider
%sigma moment b = 0 c.c.w
phi=360-theta3;
%t3-y3*g3*cosd(phi)-x3*g3*sind(phi)- y23*r3*cosd(phi)-x23*r3*sind(phi);
p3=t3-y3*g3*cosd(phi)-x3*g3*sind(phi); %el7d elmotlk 
p1=-1*r3*cosd(phi); % mo3amel y23 fe mo3adla 1 
p2=-1*r3*sind(phi); % mo3amel x23 fe mo3adla 1
%(1) p1 * y23 + p2 * x23 + p3 = 0


%link 2
%sigma moment o = 0 c.c.w
%t2+m12+y2*g2*cosd(theta2)-x2*g2*sind(theta2)+ y32*r2*cosd(theta2)-x32*r2*sind(theta2)=0
k3= t2+m12+y2*g2*cosd(theta2)-x2*g2*sind(theta2); %el7d elmotlk fe mo3adla (2)
k1=r2*cosd(theta2);  % mo3amel y23 fe mo3adla 2
k2=-1*r2*sind(theta2);  % mo3amel x23 fe mo3adla 2
%(2) k1 * y32 + k2 * x32 + k3 = 0 ; 
% y32= -1* y23 &&&&&&&& x32= -1* x23
%(2) -1*k1 * y23 + -1*k2 * x23 + k3 = 0 ;


n=[p1,p2;-k1,-k2];   %matrix mo3amelat elmgahil
e=[-p3;-k3];        %matrix el7dod elmotlka
N=n^-1 * e;           %matrix elnwateg
y23=N(1,1);
x23=N(2,1);

%sigma fx=0 right 
%x12+x2+x32=0
x32=-x23;
y32=-y23;
x12=-x2-x32;

%sigma fy=0
%y12+y2+y32=0
y12=-y2-y32;

%link3 & slider
%sigma fx=0 right
%x23+x3+(f4+p4)*cosd(alpha)-f14*cosd(90-alpha)=0
% x23+x3+ f4*cosd(alpha) + p4*cosd(alpha) - f14* cosd(90-alpha)
c3= x23+x3+ f4*cosd(alpha); %el7d elmotlk fe mo3dla 3 
c1= cosd(alpha);  %mo3amel p4 fe mo3adla 3
c2=-1* cosd(90-alpha);     %mo3amel f14 fe mo3adla 3
%(3) c1 * p4 + c2 * f14 +c3 = 0

%sigma fy=o up
%y23 + y3 + (f4+p4)*sind(alpha) +f14* sind(90-alpha);
%y23 + y3 + f4*sind(alpha)+ p4*sind(alpha) +f14* sind(90-alpha);
l3= y23 + y3 + f4*sind(alpha); %el7d elmotlk fe mo3adla 4
l1= sind(alpha);     %mo3amel p4 fe mo3adla 4
l2 = sind(90-alpha);     %mo3amel f14 fe mo3adla 4
%(4) l1 * p4 + l2 * f14 + l3 = 0


u=[c1,c2;l1,l2];  %matrix mo3amelat elmgahil
q=[-c3;-l3];       %matrix el7dod elmotlka
answer=u^-1*q;     %matrix elnwateg
p4=answer(1,1);
f14=answer(2,1);


%for the slider
%sigma fx=0 right
%x34+f4*cosd(alpha)+p4*cosd(alpha)-f14*cosd(90-alpha)=0;
x34= -1*f4*cosd(alpha)- p4*cosd(alpha) + f14*cosd(90-alpha);
x43= - x34;

%sigma fy=0 up 
%y34+f4*sind(alpha)+ p4*sind(alpha) + f14*sind(90-alpha)=0;
y34= -f4*sind(alpha) - p4*sind(alpha)-f14*sind(90-alpha);
y43= -y34;














