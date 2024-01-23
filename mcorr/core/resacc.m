function [OVERSHIFT]=resacc(junk,LAG,OXC)

CC=OXC;
PORT=CC(LAG-2:LAG+2);
P = polyfit(LAG-2:LAG+2,PORT',2);
X=LAG-2:0.001:LAG+2;
Y=((X.^2)*P(1))+(X*P(2))+P(3);

[i,j]=max(Y);
OVALMAX=i;
OLAGMAX=j;

[m,n]=min(Y);
OVALMIN=m;
OLAGMIN=n;

if abs(OVALMAX) > abs(OVALMIN)
    OLAG=OLAGMAX;
else
    OLAG=OLAGMIN;
end
OVERSHIFT=X(OLAG);
