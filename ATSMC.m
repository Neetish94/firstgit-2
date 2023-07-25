AAA=520


function [dx]= ATSMC(t,x)

t
m10=0.4;m20=1.2;L1=1;L2=0.8;J1=5;J2=5;x0=300;x1=400;x2=200;
d1=0.2*sin(3*t)+0.02*sin(26*pi*t);
d2=0.1*sin(2*t)+0.01*sin(26*pi*t);
if t<3
    m1=0.44;m2=1.32;
else if t<5
    m1=0.4;m2=1.2;
    else if t<6
            m1=0.36;m2=1;
        else if t<8
                m1=0.44;m2=1.32;
            else
                m1=0.4;m2=1.2;
            end
        end
    end
end
   M11=(m1+m2)*L1^2+m2*L2^2+2*m2*L1*L2*cos(x(2))+J1;
   M12=m2*L2^2+m2*L1*L2*cos(x(2));M22=m2*L2^2+J2;
   M110=(m10+m20)*L1^2+m20*L2^2+2*m2*L1*L2*cos(x(2))+J1;
   M120=m20*L2^2+m20*L1*L2*cos(x(2));M220=m20*L2^2+J2;
   
   C1=-m2*L1*L2*sin(x(2))*x(3)^2-2*m2*L1*L2*sin(x(2))*x(3)*x(4);
   C2=m2*L1*L2*sin(x(2))*x(4);
   
   G1=(m1+m2)*L1*cos(x(2))+m2*L2*cos(x(1)+x(2));
   G2=m2*L2*cos(x(1)+x(2));
   
   C10=-m20*L1*L2*sin(x(2))*x(3)^2-2*m20*L1*L2*sin(x(2))*x(3)*x(4);
   C20=m20*L1*L2*sin(x(2))*x(4);
   
   G10=(m10+m20)*L1*cos(x(2))+m20*L2*cos(x(1)+x(2));
   G20=m20*L2*cos(x(1)+x(2));
   
   M=[M11 M12;M12 M22];M0=[M110 M120;M120 M220];
   
   G=[G1;G2];C=[C1 C2];G0=[G10 G20];C0=[C10;C20];D=[d1;d2];
   
   
   
%COMPUTE DESIRED TRAJECTORY qd(t), qdot(t), qddot(t) 
qd1=1.25-(7/5)*exp(-t)+(7/20)*exp(-4*t);
qd2=1.4-(7/5)*exp(-t)+7/20*exp(-4*t); 
qd=[qd1;qd2];
qdot1=(7/5)*exp(-t)-(28/20)*exp(-4*t);
qdot2=(7/5)*exp(-t)-(28/20)*exp(-4*t); 
qddot=[qdot1;qdot2];
qddot1=-(7/5)*exp(-t)+(28/5)*exp(-4*t);
qddot2=-(7/5)*exp(-t)+(28/5)*exp(-4*t); 
qdddot=[qddot1;qddot2];

c=diag([2,2]); a=5;b=7;del=0.0005;
e1=[x(1);x(2)]-qd;
e2=[x(3);x(4)]-qddot;

S=e2+c*e1.^(a/b);

ueq=M0*(qdddot-(a/b)*c*(e1.^(a/b)-1))+C0+G0;
if norm(S'*inv(M0))>=del
    du=-(S'*inv(M0))'/norm(S'*inv(M0))^2*norm(S)*norm(inv(M0))*(x(5)+x(6)*norm([x(1);x(2)])+x(7)*norm([x(3);x(4)]));
else
     du=-(S'*inv(M0))'/del^2*norm(S)*norm(inv(M0))*(x(5)+x(6)*norm([x(1);x(2)])+x(7)*norm([x(3);x(4)]));
end


Det = (M11*M22)-(M12*M12); 
 
%Inverse Matrix 
Im11=M22/Det; Im12=-M12/Det; Im21=-M12/Det; Im22=M11/Det; 
INVM=[Im11 Im12;Im21 Im22];


U=ueq+du;
x11dot=x(3);
x12dot=x(4);
x21dot=Im11*(U(:,1)+d1-G1-C1)+Im12*(U(:,2)+d2-G2-C2);
x22dot=Im12*(U(:,1)+d1-G1-C1)+Im22*(U(:,2)+d2-G2-C2);
b0dot=x0*norm(S)*norm(inv(M0));
b1dot=x1*norm(S)*norm(inv(M0))*norm([x(1);x(2)]);
b2dot=x2*norm(S)*norm(inv(M0))*norm([x(3);x(4)]);
dx=[x11dot;x12dot;x21dot;x22dot;b0dot;b1dot;b2dot];


fid=fopen('ATSMCcon1.dat','a');
fprintf(fid,'%f %f \n',ueq,t);
fclose(fid);
% fid=fopen('ATSMCcon2.dat','a');
% fprintf(fid,'%f %f \n',U(:,2),t);
% fclose(fid);

    


