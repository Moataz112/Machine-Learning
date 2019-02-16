clc
close all

ds = tabularTextDatastore('house_prices_data_training_data.csv','TreatAsMissing','NA',.....
    'MissingValue',0,'ReadSize',25000);

T = read(ds);
s=size(T); %size of available data
Alpha=0.01; %typical values (0.01,0.03,0.1,0.3)

m=length(T{:,1});
U0=T{:,2}; %cost
U=T{:,4:19}; %features

U1=T{:,20:21};
U2=U.^2;
U3=sqrt(U);

X1=[ones(m,1) U U1 U.^2 U.^3];%1st hypothesis
X2=[ones(m,1) U U1 U3]; %2nd hypothesis 
X3=[ones(m,1) U U1]; %3rd hypothesis
X4=[ones(m,1) U U1 U.^3]; %4th hypothesis


n1=length(X1(1,:)); %length of features vector
n2=length(X2(1,:));
n3=length(X3(1,:));
n4=length(X4(1,:));

Y=T{:,3}/mean(T{:,3});
Theta1=zeros(n1,1);
Theta2=zeros(n2,1);
Theta3=zeros(n3,1);
Theta4=zeros(n4,1);

k=1;
j=1;
i=1;
r=1;

E1(k)=(1/(2*m))*sum((X1*Theta1-Y).^2);
E2(k)=(1/(2*m))*sum((X2*Theta2-Y).^2);
E3(k)=(1/(2*m))*sum((X3*Theta3-Y).^2);
E4(k)=(1/(2*m))*sum((X4*Theta4-Y).^2);

for w=2:n1
    if max(abs(X1(:,w)))~=0
    X1(:,w)=(X1(:,w)-mean((X1(:,w))))./std(X1(:,w));
    end
end

R=1;
while R==1

Alpha=Alpha*1;
Theta1=Theta1-(Alpha/m)*X1'*(X1*Theta1-Y);
k=k+1;
E1(k)=(1/(2*m))*sum((X1*Theta1-Y).^2); %cost function

if E1(k-1)-E1(k)<0 %error doesnot decrease
    break
end 
q=(E1(k-1)-E1(k))./E1(k-1); %percentage error
if q <.000001;
    R=0;
end
end

for w=2:n2
    if max(abs(X2(:,w)))~=0
    X2(:,w)=(X2(:,w)-mean((X2(:,w))))./std(X2(:,w));
    end
end

R=1;
while R==1

Alpha=Alpha*1;
Theta2=Theta2-(Alpha/m)*X2'*(X2*Theta2-Y);
j=j+1;
E2(j)=(1/(2*m))*sum((X2*Theta2-Y).^2); %cost function

if E2(j-1)-E2(j)<0 %error doesnot decrease
    break
end 
q=(E2(j-1)-E2(j))./E2(j-1); %percentage error
if q <.000001;
    R=0;
end
end

for w=2:n3
    if max(abs(X3(:,w)))~=0
    X3(:,w)=(X3(:,w)-mean((X3(:,w))))./std(X3(:,w));
    end
end

R=1;
while R==1

Alpha=Alpha*1;
Theta3=Theta3-(Alpha/m)*X3'*(X3*Theta3-Y);
i=i+1;
E3(i)=(1/(2*m))*sum((X3*Theta3-Y).^2); %cost function

if E3(i-1)-E3(i)<0 %error doesnot decrease
    break
end 
q=(E3(i-1)-E3(i))./E3(i-1); %percentage error
if q <.000001;
    R=0;
end
end

for w=2:n4
    if max(abs(X4(:,w)))~=0
    X4(:,w)=(X4(:,w)-mean((X4(:,w))))./std(X4(:,w));
    end
end

R=1;
while R==1

Alpha=Alpha*1;
Theta4=Theta4-(Alpha/m)*X4'*(X4*Theta4-Y);
r=r+1;
E4(r)=(1/(2*m))*sum((X4*Theta4-Y).^2); %cost function

if E4(r-1)-E4(r)<0 %error doesnot decrease
    break
end 
q=(E4(r-1)-E4(r))./E4(r-1); %percentage error
if q <.000001;
    R=0;
end
end

figure(1)
plot(1:k,E1)
title('1st hypothesis')

figure(2)
plot(1:j,E2)
title('2nd hypothesis')

figure(3)
plot(1:i,E3)
title('3rd hypothesis')

figure(4)
plot(1:r,E4)
title('4th hypothesis')
