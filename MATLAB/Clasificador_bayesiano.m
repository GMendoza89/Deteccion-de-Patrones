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

%clc;
%clear;
close all;
%% Lectura de datos
ClassN = 4;

MC1 = load('MCor_C1.mat');
MC2 = load('MCor_C2.mat');
MC3 = load('MCor_C3.mat');
MC4 = load('MCor_C4.mat');

VM1 = load('VMean1.mat');
VM2 = load('VMean2.mat');
VM3 = load('VMean3.mat');
VM4 = load('VMean4.mat');

MC1L = load('MCorL_C1.mat');
MC2L = load('MCorL_C2.mat');
MC3L = load('MCorL_C3.mat');
MC4L = load('MCorL_C4.mat');

VM1L = load('VMeanL1.mat');
VM2L = load('VMeanL2.mat');
VM3L = load('VMeanL3.mat');
VM4L = load('VMeanL4.mat');

PA = load('PA.mat');

%% Lectura de Imagenes 
if ispc
    d='\';
else
    d='/';
end

pht1 = [pwd d , 'Validacion', d];
fl1=dir([pht1  '*.jpg']);


N=[3000;% numero de puntos a extraer
   13];

Nm = 50;
nfl=length(fl1);
img_rgb = cell(nfl,1);
img_hsl = cell(nfl,1);
vecM = zeros(nfl,2);
vecM1= zeros(nfl,3);
PB=zeros(nfl,5);
PT=zeros(nfl,4);
vecML = zeros(nfl,2);
vecM1L= zeros(nfl,3);
PBL=zeros(nfl,5);
PTL=zeros(nfl,4);
for i1 = 1:nfl
    img_rgb{i1} = imread([pht1 fl1(i1).name]);
    [lv,lu,ch] = size(img_rgb{i1});
    img_hsl{i1} = rgb2hsl(img_rgb{i1});
    
    mxu = (lu - 1-N(2,1));
    mxv = (lv - 1-N(2,1));
    mn = (1+N(2,1));
   
    pts = zeros(N(1,1),2);
    pts(:,1) = round((lu - 1)*rand(N(1,1),1)+1);
    pts(:,2) = round((lv - 1)*rand(N(1,1),1)+1);
    OWu = floor((mxu-mn)*rand()+mn);
    OWv = floor((mxv-mn)*rand()+mn);
            
    Window =  img_hsl{i1}(OWv-floor(N(2,1)/2):OWv+floor(N(2,1)/2),OWu-floor(N(2,1)/2):OWu+floor(N(2,1)/2),3);

    vecM(i1,1)=mean(mean(img_hsl{i1}(pts(:,2),pts(:,1),1)));
    vecM(i1,2)=mean(mean(img_hsl{i1}(pts(:,2),pts(:,1),2)));

    vecM1L(i1,1)=mean(var(Window));
    vecM1L(i1,2)=mean(kurtosis(Window));
    vecM1L(i1,3)=mean(skewness(Window));
    PB(i1,1)=mvnpdf(vecM(i1,:),VM1.Vector,MC1.Matrix)*PA.PAi(1);
    PB(i1,2)=mvnpdf(vecM(i1,:),VM2.Vector,MC2.Matrix)*PA.PAi(2);
    PB(i1,3)=mvnpdf(vecM(i1,:),VM3.Vector,MC3.Matrix)*PA.PAi(3);
    PB(i1,4)=mvnpdf(vecM(i1,:),VM4.Vector,MC4.Matrix)*PA.PAi(4);
    PB(i1,5)=PB(i1,1)+PB(i1,2)+PB(i1,3)+PB(i1,4);
    
    PBL(i1,1)=mvnpdf(vecM1L(i1,:),VM1L.VectorL,MC1L.MatrixL);
    PBL(i1,2)=mvnpdf(vecM1L(i1,:),VM2L.VectorL,MC1L.MatrixL);
    PBL(i1,3)=mvnpdf(vecM1L(i1,:),VM3L.VectorL,MC1L.MatrixL);
    PBL(i1,4)=mvnpdf(vecM1L(i1,:),VM4L.VectorL,MC1L.MatrixL);
    PBL(i1,5)=PBL(i1,1)+PBL(i1,2)+PBL(i1,3)+PBL(i1,4);
    
    PBL(i1,:)=PBL(i1,:)/PBL(i1,5);
    
    [M,Id] = max(PB(i1,1:4));
    [ML,IdL] = max(PBL(i1,1:4));
    
    fprintf('Las propiedades de color de la figura, %s', fl1(i1).name);
    fprintf(' pertenece a la clase %d',Id);
    fprintf( ' Textura clase %d  \n',IdL);
    
end




