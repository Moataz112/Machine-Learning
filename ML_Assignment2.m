clc
close all

ds = tabularTextDatastore('heart_DD.csv','TreatAsMissing','NA',.....
    'MissingValue',0,'ReadSize',25000);

T = read(ds);
s=size(T); %%size of available data
Alpha=0.03; %typical values (0.01,0.03,0.1,0.3)


m=length(T{:,1});
U0=T{:,14};
U=T{:,1:13};
U2=U.^2;
U3=sqrt(U);

%Note:
%I had no time to copy paste and add the hypothesis so i commented them as
%shown below to easily switch and get different graphs for different X
%uncommenting each X and commenting the previous one will be correct

X=[ones(m,1) U U.^2 U.^3]; %1st hypothesis
%X=[ones(m,1) U U3]; %2nd hypothesis 
%X=[ones(m,1) U U.^2]; %3rd hypothesis
%X=[ones(m,1) U]; %4th hypothesis

n=length(X(1,:));

for w=2:n
    if max(abs(X(:,w)))~=0
    X(:,w)=(X(:,w)-mean((X(:,w))))./std(X(:,w));
    end
end

Y=T{:,14};
Theta=zeros(n,1);
k=1;
G=(1./(1+exp(-X*Theta)));
E(k)=(1/(m))*sum((-Y.'*log(1./(1+exp(-X*Theta))))-((1-Y).'*log(1-(G)))); %cost function

R=1;
while R==1

Alpha=Alpha*1;
Theta=Theta-(Alpha/m)*X'*((1./(1+exp(-X*Theta)))-Y);
k=k+1;
E(k)=(1/(m))*sum((-Y.'*log(1./(1+exp(-X*Theta))))-((1-Y).'*log(1-(1./(1+exp(-X*Theta))))));

if E(k-1)-E(k)<0
    break
end 
q=(E(k-1)-E(k))./E(k-1);
if q <.000001;
    R=0;
end
end

figure(1)
plot(1:k,E)
title('hypothesis')