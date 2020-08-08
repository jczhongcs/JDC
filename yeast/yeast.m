clc
clear
close all
%prepare
data=xlsread('Yeast_topology.xlsx');%topology
data1=xlsread('Yeast_PPI_network.xlsx');%PPI network
data3=xlsread('geneexpress.xlsx');%gene expression
for i=1:36
   data2(:,i)=data3(:,i);
end

a=zeros(5093,5093);
for row = 1:length(data1(:,1))
   a(data1(row,1),data1(row,2)) = 1;    
   a(data1(row,2),data1(row,1)) = 1;
end
degree=full(sum(a,2));
degree_gene=full(sum(data2,2));
for i=1:5093
    E_mean(i,1)=mean(data2(i,1:36),2);
    E_sigma(i,1)=std(data2(i,1:36));
    E_var(i,1)=var(data2(i,1:36));
    T(i,1)=E_var(i,1)/(E_var(i,1)+1);
    BD(i,1)=E_sigma(i,1)/E_mean(i,1);
    G(i,1)=E_mean(i,1)+2*E_sigma(i,1)*T(i,1);
end
A=zeros(5093,36);
for i=1:5093
   for j=1:36
      if data2(i,j)>G(i,1)
         A(i,j)=1; 
      end
   end
end
%calculate Pearson Jaccard,Ecc,pec,pec_j
 Ecc=zeros(5093,5093);
for i=1:5093
    for j=1:5093
       if a(i,j)==1 && i~=j 
           P(i,j)=corr(data3(i,1:36)',data3(j,1:36)','type','pearson');%pearson
           D=[A(i,:);A(j,:)];
           J(i,j)=1-pdist(D,'Jaccard');%Jaccard
            temp=0;
           for k=1:5093
              if a(i,k)==1 && a(j,k)==1 && i~=k~=j
                  temp=temp+1;
              end
           end
           z(i,j)=temp;
           m(i,j)=min(degree(i,1)-1,degree(j,1)-1);
           Ecc(i,j)=z(i,j)./m(i,j);
           PeC(i,j)=Ecc(i,j)*P(i,j);%PeC
           PeC_J(i,j)=PeC(i,j)*J(i,j);%PeC_J      
      end
    end
end
Ecc(isnan(Ecc)==1)=0;
PeC(isnan(PeC)==1)=0;
P(isnan(P)==1)=0;
J(isnan(J)==1)=0;
PeC_J(isnan(PeC_J)==1)=0;
pec_j=full(sum(PeC_J,2));
pec=full(sum(PeC,2));
ecc=full(sum(Ecc,2));
p=full(sum(P,2));
j=full(sum(J,2));
for i=1:6
   all_p(:,i)=data(:,i).*p; 
   all_j(:,i)=data(:,i).*j;
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
    WDC=((i*Ecc)+((1-i)*P));
    WDC(isnan(WDC)==1)=0;
    wdc=full(sum(WDC,2));
    F=[b0,wdc]; 
    b0=F; 
end