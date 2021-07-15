%%Bifurcation code
clc;
clear all

len1=250; 
len2=250;
trails=1;

cpy_ei=linspace(0.1,0.9,len1);
ci1_ei=linspace(0.2,0.9,len2);


State=zeros(len1,len2);
FD=zeros(len1,len2);

for i=1:len1
    i
    for j=1:len2
        [p1,p2,p3,p4,p5,p6]=state(cpy_ei(i),ci1_ei(j));
        State(i,j)=p5;
        FD(i,j)=p6;
    end
end

save datatFig2DE
%%
figure(1)
imagesc(ci1_ei,cpy_ei,State),
set(gca,'ydir','normal'),colorbar,colormap(jet);
axis on;
xlabel('Ci1_ei','FontSize',20);
ylabel('Cpy_ei','FontSize',20);

figure(2)
imagesc(ci1_ei,cpy_ei,FD),
set(gca,'ydir','normal'),colorbar,colormap(jet);
axis on;
xlabel('Ci1_ei','FontSize',20);
ylabel('Cpy_ei','FontSize',20);


