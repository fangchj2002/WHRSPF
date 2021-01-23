%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo of "Active contour driven by weighted hybrid region-based signed 
%  pressure force"(WHRSPF)
% Jiangxiong Fang
% East China University of Technology&&Nanchang University, Nanchang, China
% 8th, May, 2019
% Email: fangchj2002@163.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear all;close all;

%addpath 'D:\resFiles\Compared ACMs\image_Database\results'
addpath 'image'

index =3;

switch index
    case 1
        Img =imread('8.bmp');
        stat = 30;
        sigma = 1.5;
        rad = 10;
    case 2
        Img = imread('2.pgm');
        stat = 30;
        sigma = 1.5;
        rad = 5;
    case 3
        Img = imread('3.pgm');
        stat = 55;
        sigma = 1.2;
        rad = 5;
    case 4
        Img = imread('4.pgm');
        stat = 5;
        sigma = 1.5;
        rad = 5;
    case 5
        Img = imread('5.bmp');
        stat = 10;
        sigma = 1.5;
        rad = 5;
end

if size(Img,3)>1
    Img = rgb2gray(Img);
end


[row,col] = size(Img);
phi = ones(row,col);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% you can add different types of noise.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Img = imnoise(Img,'speckle',0.03);

phi(stat:row-stat,stat:col-stat) = -1;

u = - phi;
figure, subplot(2,2,1);imshow(Img);hold on;
[c, h] = contour(u, [0 0], 'r','LineWidth',2);
title('Initial contour');
hold off;
Img_a=Img;
Img = double(Img);
tic;
G = fspecial('gaussian',rad, sigma);

Iter = 100;
% the weights of weighted global and local region-based SPF
wg = 1;
wl = 0.1;
subplot(2,2,2);
for n = 1:Iter
    [ux, uy] = gradient(u);
    mu = sqrt(ux.^2 + uy.^2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Weighted global region-based SPF (GRSPF)    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    c1 = sum(sum(Img.*(u<0)))/(sum(sum(u<0)));
    c2 = sum(sum(Img.*(u>=0)))/(sum(sum(u>=0)));
    m1 = median(Img(u<0));
    m2 = median(Img(u>=0));    
    wg1= sum(sum((u>0)))/(row*col);
    wg2= sum(sum((u<0)))/(row*col);

    grspf=Img-((wg1*c1^2+wg2*m1^2-c2^2)/(wg1*2*c1+wg2*2*m1-2*c2));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Weighted local region-based SPF (GLSPF)    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nImg = conv2(Img,G,'same');
    u1=conv2(double(u>0),G,'same');
    s1=conv2(Img.*(u>0),G,'same'); 
    u2=conv2(double(u<0),G,'same');
    s2=conv2(double(Img.*(u<0)),G,'same'); 
    
    wl1= sum(sum((nImg-s1)>0))/(row*col);
    wl2= sum(sum((nImg-s2)>0))/(row*col);
  
    f1 = sum(s1(:))./sum(u1(:));
    f2 = sum(s2(:))./sum(u2(:));
    lrspf = nImg - (wl1*f1 + wl2*f2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Weighted hybrid region-based SPF (WHSPF)    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    gmax=max(abs(grspf(:)));
    lmax=max(abs(lrspf(:)));
    if gmax>lmax
        grspf=grspf/gmax*lmax;
    else 
        lrspf=lrspf/lmax*gmax;
    end
    hrspf = wg.*grspf+wl.*lrspf;    
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Computing the level set function (LSF)    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    u = u + abs(c1+m1-2*c2)*hrspf.*mu;
    
    if mod(n,10)==0
    imshow(Img_a); hold on;
    [c, h] = contour(u, [0 0], 'r');
    iterNum = [num2str(n), 'iterations'];
    title(iterNum);
    pause(0.02);
    end
    u = (u >= 0) - ( u< 0);% the selective step.
    objPos=bwareaopen(u>0,20);
    objPos=imclose(objPos,strel('disk',3)); objNeg = ~objPos;  
    u = objPos - objNeg;
    u = conv2(u, G, 'same');
end
imshow(Img_a);hold on;
[c, h] = contour(u, [0 0], 'r','LineWidth',2);
subplot(2,2,3);
toc
seg = u>0;
imshow(seg);
