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
%             Bayes - Practica 01
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

N=[3000;% numero de puntos a extraer
   13]; % Anchos de Ventanas
Nm = [20 30 40 50 60 70 80 90 100 110];

ClassN = 4;
ImgxClass=10;
ImgN = ClassN*ImgxClass; % numero de imagenes procesadas 

nfl=length(fl1);
img_rgb = cell(nfl,1);
img_hsl = cell(nfl,1);
vecM = cell(nfl,length(Nm));
vecM1= cell(nfl,length(Nm));
% col = [1 0 0; 0 1 0; 0 0 1; 1 1 0];
%% Cliclo de Generacion de puntos
for i1 = 1:nfl
    img_rgb{i1} = imread([pht1 fl1(i1).name]);
    [lv,lu,ch] = size(img_rgb{i1});
    img_hsl{i1} = rgb2hsl(img_rgb{i1});
    
    mxu = (lu - 1-N(2,1));
    mxv = (lv - 1-N(2,1));
    mn = (1+N(2,1));
    Window = zeros(N(2,1));
    
    for i2=1:length(Nm)
        for i3=1:Nm(i2)
            pts = zeros(N(1,1),2);
            pts(:,1) = round((lu - 1)*rand(N(1,1),1)+1);
            pts(:,2) = round((lv - 1)*rand(N(1,1),1)+1);
            OWu = floor((mxu-mn)*rand()+mn);
            OWv = floor((mxv-mn)*rand()+mn);
            
            Window =  img_hsl{i1}(OWv-floor(N(2,1)/2):OWv+floor(N(2,1)/2),OWu-floor(N(2,1)/2):OWu+floor(N(2,1)/2),3);
            
            vecM{i1,i2}(1,i3)=mean(mean(img_hsl{i1}(pts(:,2),pts(:,1),1)));
            vecM{i1,i2}(2,i3)=mean(mean(img_hsl{i1}(pts(:,2),pts(:,1),2)));
            
            vecM1{i1,i2}(1,i3)=mean(var(Window));
            vecM1{i1,i2}(2,i3)=mean(kurtosis(Window));
            vecM1{i1,i2}(3,i3)=mean(skewness(Window));      
        end
    end
    fprintf('Generacion de muestras: %2.2f porciento \n',i1/nfl*100);
    
end
fprintf('Generacion de muestras completado ;) \n');
%%
VecCarxC = cell(ClassN,length(Nm));
Vmean = cell(ClassN,length(Nm));
Mcor = cell(ClassN,length(Nm));
HCaxC = zeros(ClassN,length(Nm));
Mi = cell(ClassN,1);

VecCarxCL = cell(ClassN,length(Nm));
VmeanL = cell(ClassN,length(Nm));
McorL = cell(ClassN,length(Nm));
HCaxCL = zeros(ClassN,length(Nm));
MiL = cell(ClassN,1);
for i4=1:ClassN
    Class = ceil(i4/ImgxClass);
    
    for i5=1:length(Nm)
        VP = [];
        VPL = [];
        for i6=1:ImgxClass
            VP = [VP vecM{ImgxClass*(i4-1)+i6,i5}];
            VPL = [VPL vecM1{ImgxClass*(i4-1)+i6,i5}];
        end
        VecCarxC{i4,i5}=VP';
        Vmean{i4,i5} = mean(VP');
        Mcor{i4,i5} = cov(VP');
        HCaxC(i4,i5) = (1/(2*pi*det(Mcor{i4,i5})));
        
        VecCarxCL{i4,i5}=VPL';
        VmeanL{i4,i5} = mean(VPL');
        McorL{i4,i5} = cov(VPL');
        HCaxCL(i4,i5) = (1/(pow2(2*pi,1.5)*det(McorL{i4,i5})));
    end
    [M,Id] =max(HCaxC(i4,:));
    Mi{i4} = [M,Id];
    [ML,IdL] =max(HCaxCL(i4,:));
    MiL{i4} = [ML,IdL];
    
end
fprintf('Entrenamiento Completado ... ;) \n');
%%
tA = 0;
PAi=zeros(ClassN,1);
PAiL=zeros(ClassN,1);

for i7 = 1:ClassN
    tA = tA + Nm(Mi{i7}(2));
    PAi(i7) = Nm(Mi{i7}(2));
end
PAi = PAi/tA;

save('P.mat','PAi')

for i8 = 1:ClassN
    FnM=sprintf('MCor_C%d.mat',i8);
    Matrix =Mcor{i8, Mi{i8}(2)};
    FnMean=sprintf('VMean%d.mat',i8);
    Vector =Vmean{i8, Mi{i8}(2)};
    save(FnM,'Matrix');
    save(FnMean,'Vector');
    
    FnML = sprintf('MCorL_C%d.mat',i8);
    MatrixL = McorL{i8, MiL{i8}(2)};
    FnMeanL = sprintf('VMeanL%d.mat',i8);
    VectorL = VmeanL{i8, MiL{i8}(2)};
    save(FnML,'MatrixL');
    save(FnMeanL,'VectorL');
    
end
fprintf('Datos Guardados ... ;) Buenas Noches Gustavo \n');
