clc
clear
close all
data2=xlsread('gene_ecoli.xlsx');
data1=xlsread('DIP.xlsx');
data3=xlsread('ECOLI(7).xlsx');
for i=1:8
   data(:,i)=mapminmax(data2(:,i)',0,1)';
end

a=zeros(2727,2727);
for row = 1:length(data1(:,1))
   a(data1(row,1),data1(row,2)) = 1;    
   a(data1(row,2),data1(row,1)) = 1;
end
degree=full(sum(a,2));
for i=1:2727
    E_mean(i,1)=mean(data(i,:),2);
    E_sigma(i,1)=std(data(i,:));
    E_var(i,1)=var(data(i,:));
    T(i,1)=E_var(i,1)/(E_var(i,1)+1);
    G(i,1)=E_mean(i,1)+2*E_sigma(i,1)*T(i,1);                                                                                                        
end

A=zeros(2727,8);
for i=1:2727
   for j=1:8
      if data(i,j)>G(i,1)
         A(i,j)=1; 
      end
   end
end
P=zeros(2727,2727);
J=zeros(2727,2727);
Ecc=zeros(2727,2727);
for i=1:2727
   for j=1:2727
      if a(i,j)==1 && i~=j
         P(i,j)=corr(data(i,1:8)',data(j,1:8)','type','pearson');
         D=[A(i,:);A(j,:)];
         J(i,j)=1-pdist(D,'Jaccard');
          
        temp=0;
           for k=1:2727
              if a(i,k)==1 && a(j,k)==1 && i~=k~=j
                  temp=temp+1;
              end
           end
           z(i,j)=temp;
           m(i,j)=min(degree(i,1)-1,degree(j,1)-1);
           Ecc(i,j)=z(i,j)./m(i,j);
           Pec(i,j)=Ecc(i,j)*P(i,j);
           PeC_J(i,j)=Pec(i,j)*J(i,j); 
      end
   end
end
J(isnan(J)==1)=0;
P(isnan(P)==1)=0;
Ecc(isnan(Ecc)==1)=0;
Pec(isnan(Pec)==1)=0;
PeC_J(isnan( PeC_J)==1)=0;
ecc=full(sum(Ecc,2));
pec=full(sum(Pec,2));
pec_j=full(sum( PeC_J,2));
j=full(sum(J,2));
p=full(sum(P,2));
for i=4:9
   all_j(:,i-3)=data3(:,i).*j;
   all_p(:,i-3)=data3(:,i).*p; 

end
a0=[];
for i=linspace(0,1,11); 
    WDC_J=((i*Ecc)+((1-i)*P)).*J;
    WDC_J(isnan(WDC_J)==1)=0;
    wdc_j=full(sum(WDC_J,2));
    E=[a0,wdc_j]; 
    a0=E; 
end
b0=[];
for i=linspace(0,1,11); 
    wdc=((i*Ecc)+((1-i)*P));
    wdc(isnan(wdc)==1)=0;
    wdc=full(sum(wdc,2));
    F=[b0,wdc]; 
    b0=F; 
end 

