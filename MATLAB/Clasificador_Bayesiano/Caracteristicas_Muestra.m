% 
%           Universidad de Guanajuato
%   Division de Ingenierias Campus Irapuato Salamanca
%   
%         Maestria en Ingenieria Eléctrica
%   
%                Detección de Patrones
%         
%            Mendoza Pinto Gustavo David 
%         
%             Selector de Caracteristicas - Practica 01
%           

clc;
clear;
close all;
%% Lectura de Imagenes 
if ispc
    d='\';
else
    d='/';
end

pht1 = [pwd d , 'Imagenes', d];
fl1=dir([pht1  '*.jpg']);

N=[1000 2000 3000 4000 5000 6000 7000 8000 9000 10000;% numero de puntos a extraer
    5 7 9 11 13 15 17 19 21 23];              %WW = [5 7 9 11 13 15 17 19 21 23];            % Anchos de Ventanas
Nm = 50;

nfl=length(fl1);
img_rgb = cell(nfl,1);
img_hsl = cell(nfl,1);
vecM = cell(nfl,1);
vecM1= cell(nfl,1);
% col = [1 0 0; 0 1 0; 0 0 1; 1 1 0];
%% Cliclo de Generacion de puntos
for i1 = 1:nfl
     %Cclass = ceil(i1/4);
    img_rgb{i1} = imread([pht1 fl1(i1).name]);
    [lv,lu,ch] = size(img_rgb{i1});
    img_hsl{i1} = rgb2hsl(img_rgb{i1});
%     figure
%     imshow(img_rgb{i1})
%     figure
%    imshow(img_hsl{i1})
    for i2=1:length(N)  % inicia generacion de puntos y ventanas
        mxu = (lu - 1-N(2,i2));
        mxv = (lv - 1-N(2,i2));
        mn = (1+N(2,i2));
        Window = zeros(N(2,i2));
        for i3=1:Nm
            % generar N puntos aleatorios
            pts = zeros(N(1,i2),2);
            pts(:,1) = round((lu - 1)*rand(N(1,i2),1)+1);
            pts(:,2) = round((lv - 1)*rand(N(1,i2),1)+1);
            OWu = floor((mxu-mn)*rand()+mn);
            OWv = floor((mxv-mn)*rand()+mn);
            
            Window =  img_hsl{i1}(OWv-floor(N(2,i2)/2):OWv+floor(N(2,i2)/2),OWu-floor(N(2,i2)/2):OWu+floor(N(2,i2)/2),3);
            
            vecM{i1}(1,Nm*(i2-1)+i3)=mean(mean(img_hsl{i1}(pts(:,2),pts(:,1),1)));
            vecM{i1}(2,Nm*(i2-1)+i3)=mean(mean(img_hsl{i1}(pts(:,2),pts(:,1),2)));
            
            vecM1{i1}(1,Nm*(i2-1)+i3)=mean(var(Window));
            vecM1{i1}(2,Nm*(i2-1)+i3)=mean(kurtosis(Window));
            vecM1{i1}(3,Nm*(i2-1)+i3)=mean(skewness(Window));
            
        end
        
%         figure(2)
%         hold on;
%         plot(vecM{i1}(2,:), vecM{i1}(1,:), 'Marker','*', 'Color', col(Cclass,:), 'LineStyle', 'none');
%         xlabel('H');
%         ylabel('S');
%         grid on;
    end
    fprintf('Porcentaje de generacion de puntos: %2.2f porciento \n',i1/nfl*100);
end

%% Binary tree
vp =[];
vpl=[];
for i4=1:nfl
    vp = [vp,vecM{i4}];
    vpl = [vpl,vecM1{i4}];
end

 iG= [];
 iGL=[];
 vp = vp';
 vpl = vpl';
 NZ = length(vp);
 NZL = length(vpl);
 th= (NZ)*0.9;
 thl=(NZL)*0.9;
 
 Z = linkage(vp, 'single', 'euclidean');
 Zl = linkage(vpl, 'single', 'euclidean');
%  %%
% figure
% dendrogram(Z)
% figure
% dendrogram(Zl)
% puntos en el arbol que no son grupos compuestos
id1 = Z(:,1) <= NZ;
id2 = Z(:,2) <= NZ;
id1l = Z(:,1) <= NZL;
id2l = Z(:,2) <= NZL;
id = id1 + id2; % suma los puntos que no son grupos
idl = id1l + id2l; % suma los puntos que no son grupos
for i5=2:NZ-1 % acumulando numero de puntos por nivel
    id(i5)= id(i5)+id(i5-1);
    idl(i5)= idl(i5)+idl(i5-1);
end
id3 = find(id >= th);
id3l = find(idl >= th);
% nivel de Z con th puntos
Zniv = Z(1:id3(1),:);
Zlniv = Zl(1:id3l(1),:);
% recuperar los puntos validando no grupos compuestos
id1 = Zniv(:,1)<=NZ;
id2 = Zniv(:,2)<=NZ;
id1l = Zlniv(:,1)<=NZL;
id2l = Zlniv(:,2)<=NZL;
%% Relacion de caracteristicas
Vres = [vp(Zniv(id1,1),1) vp(Zniv(id1,1),2); vp(Zniv(id2,1),1) vp(Zniv(id2,1),2)];
VresL = [vpl(Zlniv(id1l,1),1) vpl(Zlniv(id1l,1),2); vpl(Zlniv(id2l,1),1) vpl(Zlniv(id2l,1),2)];

Ind_Vec = [Zniv(id1,1); Zniv(id2,1)];
Ind_Vecl = [Zlniv(id1,1); Zlniv(id2,1)];

GHS=zeros(4,10);
GL = zeros(4,10);
%%
pxC =length(vp)/4;
pxIm= pxC/(10);
pxM = pxIm/(length(N));
%%
for i6 = 1:nfl
    class = ceil(i6/10);
    
    for i7 = 1:length(N)
        LimInf = round(pxIm*(i6-1)+pxM*(i7-1)+1)
        LimSp = round(pxIm*(i6-1)+pxM*(i7))
        for i8=1:length(Ind_Vec(:,1))
            if(Ind_Vec(i8,1) >= LimInf && Ind_Vec(i8,1)<=LimSp)
                GHS(class,i7)=GHS(class,i7)+1;
            end
            if(Ind_Vecl(i8,1)>= LimInf && Ind_Vecl(i8,1)<=LimSp)
                GL(class,i7)=GL(class,i7)+1;
            end
        end
    end
end

%%
TP=sum(GHS*ones(10,1));
TPxC = sum(GHS');
TPxM = sum(GHS);
pTPxC = sum(GHS')/TP;
pTPxM = sum(GHS)/TP;

TPL=sum(GL*ones(10,1));
TPLxC = sum(GL');
TPLxM = sum(GL);
pTPLxC = sum(GHS')/TPL;
pTPLxM = sum(GHS)/TPL;

pxCxM=zeros(4,10);
pLxCxM = zeros(4,10);

for i9 = 1:4
    pxCxM(i9,:)=GHS(i9,:)/TPxC(i9);
    pLxCxM(i9,:)=GL(i9,:)/TPLxC(i9);
    
end

%%
fileID = fopen('DatosR1.txt','w');
fprintf(fileID,"Datos Relacionados a Muestreos \n");
fprintf(fileID,"Porcentaje de datos concervados por clase para H y S  \n");
fprintf(fileID,' %2.4f %2.4f %2.4f %2.4f \n',pTPxC*100);
fprintf(fileID,"Porcentaje de datos concervados por muestra para H y S  \n");
fprintf(fileID,' %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f\n',pTPxM*100);
fprintf(fileID,"Porcentaje de datos concervados por clase y por muestra para H y S  \n");
fprintf(fileID,'           %d      %d      %d      %d      %d      %d      %d      %d      %d      %d \n',N(1,:));
for i10=1:4
    fprintf(fileID,'clase %d:  ',i10);
    
    fprintf(fileID,' %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f\n',pxCxM(i10,:)*100);
end
fprintf(fileID,"Porcentaje de datos concervados por clase para Ventanas \n");
fprintf(fileID,' %2.4f %2.4f %2.4f %2.4f \n',pTPLxC*100);
fprintf(fileID,"Porcentaje de datos concervados por muestra para H y S  \n");
fprintf(fileID,' %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f\n',pTPxM*100);
fprintf(fileID,"Porcentaje de datos concervados por clase y por muestra para H y S  \n");
fprintf(fileID,'             %d       %d      %d     %d    %d     %d     %d      %d     %d      %d \n',N(2,:));
for i11=1:4
    fprintf(fileID,'clase %d:  ',i11);
    fprintf(fileID,' %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f\n',pLxCxM(i10,:)*100);
end

fclose(fileID);








