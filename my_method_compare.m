% for public release
clear
close all;
clc;

PatchSize=8;
showflag=0;

addpath(genpath('sparsefusion'));
addpath(genpath('ksvdbox'))


[imagename1, imagepath1]=uigetfile('LytroDataset\*.jpg;*.bmp;*.png;*.tif;*.tiff;*.pgm;*.gif','Please choose the first input image');
image_input1=imread(strcat(imagepath1,imagename1));    
[imagename2, imagepath2]=uigetfile('LytroDataset\*.jpg;*.bmp;*.png;*.tif;*.tiff;*.pgm;*.gif','Please choose the second input image');
image_input2=imread(strcat(imagepath2,imagename2));    



if size(image_input1)~=size(image_input2)
    error('two images are not the same size.');
end


img1=double(image_input1);
img2=double(image_input2);

% img11=impad(img1,PatchSize);
% img22=impad(img2,PatchSize);

overlap = 6;                    
epsilon=0.1;
%%

load('sparsefusion\mydic\mylearnd.mat')

    


timeprifix=datestr(datetime);
timeprifix(isspace(timeprifix)) = [];
timeprifix([end-2 end-5])=[];
%   compare method dchwt

%%
if size(img1,3)==1   %for gray images
    im_nsct_sr = lp_sr_fuse(img1,img2,level,3,3,D,overlap,epsilon);  % method of lapalacian and sr
    imwrite(uint8(im_nsct_sr),['myresults/' timeprifix 'im_lp_sr.bmp']);
    
    [imgf,tscore,fullscore]=mysparse_fusion(img1,img2,D,overlap,epsilon);
    [mysc,scerode,maxregion]=xingtaixue(fullscore);
    myif=img1.*mysc+(1-mysc).*img2;
    imwrite(uint8(myif),['myresults/' timeprifix 'myif.bmp']);
    
else                 %for color images
    imgf=zeros(size(img1));
    mysc=zeros(size(img1));
    im_nsct_sr=zeros(size(img1));
    im_csr=zeros(size(img1));
    im_sr=zeros(size(img1));
   %*********************************************************************** 
    for i=1:3
    im_sr(:,:,i)=sparse_fusion(img1(:,:,i),img2(:,:,i),D,overlap,epsilon);
    end
    imwrite(uint8(im_sr),['myresults/' timeprifix 'im_mwg.bmp']);
    fprintf('base_SR done\n')
  %************************************************************************  
    if_dchwt=mydchwt(img1,img2);
    imdchwt=uint8(if_dchwt);
    imwrite(imdchwt,['myresults/' timeprifix 'imdchwt.bmp']);
    fprintf('DCHWT done\n')
   %***********************************************************************
    gffinput=cell(1,2);
    gffinput{1,1}=img1;
    gffinput{1,2}=img2;
    im_gff=gff(gffinput,0);
    imwrite(im_gff,['myresults/' timeprifix 'imgff.bmp']);
    fprintf('GFF done\n')
   %***********************************************************************
    for i=1:3
      im_nsct_sr(:,:,i) = nsct_sr_fuse(img1(:,:,i),img2(:,:,i),[2,3,3,4],D,overlap,epsilon);  % method of lapalacian and sr
    end
    imwrite(uint8(im_nsct_sr),['myresults/' timeprifix 'im_nsct_sr.bmp']);
    fprintf('NSCT_SR done\n')
   %***********************************************************************
    load('csrdic.mat')
    for i=1:3
        im_csr(:,:,i) = CSR_Fusion(img1(:,:,i),img2(:,:,i),D_CSR,0.01,1);  % method of csr
    end
    imwrite(uint8(im_csr),['myresults/' timeprifix 'im_csr.bmp']);
    fprintf('CSR done\n')
   %***********************************************************************
    for i=1:3   
        [fullscore]=my_if(img1(:,:,i),img2(:,:,i),D,overlap,epsilon);
        [mysc(:,:,i),scerode,maxregion]=xingtaixue(fullscore); 
    end
    myif=img1.*mysc+(1-mysc).*img2;
    imwrite(uint8(myif),['myresults/' timeprifix 'myif.bmp']);
    fprintf('proposed method done\n')
   %***********************************************************************
    
%     model_name = 'model/cnnmodel.mat';
%     for i=1:3
%     im_cnn(:,:,i)=CNN_Fusion(img1(:,:,i),img2(:,:,i),model_name);
%     end
%     imwrite(uint8(im_cnn),['myresults/' timeprifix 'im_cnn.bmp']);

    
end


    




 if(size(img1,3)==1)
     MetricY=metricYang(img1,img2,imgf);
     MetricY1=metricYang(img1,img2,myif);
     MetricXy=metricXydeas(img1,img2,imgf);
     MetricXy1=metricXydeas(img1,img2,myif);
 else
     grim1=rgb2gray(uint8(img1));
     grim2=rgb2gray(uint8(img2));
     
     grimf2=rgb2gray(uint8(imdchwt));
     grimf6=rgb2gray(uint8(myif));
     grimf4=rgb2gray(uint8(im_nsct_sr));
     grimf5=rgb2gray(uint8(im_csr));
     grimf1=rgb2gray(uint8(im_sr));
     grimf3=rgb2gray(im_gff);
     
     qabfdchwt=metricXydeas(grim1,grim2,grimf2);
     qabfmyif=metricXydeas(grim1,grim2,grimf6);
     qabfnsctsr=metricXydeas(grim1,grim2,grimf4);
     qabfcsr=metricXydeas(grim1,grim2,grimf5);
     qabfsr=metricXydeas(grim1,grim2,grimf1);
     qabfgff=metricXydeas(grim1,grim2,grimf3);
     
     qnmidchwt=metricMI(grim1,grim2,grimf2);
     qnmimyif=metricMI(grim1,grim2,grimf6);
     qnminsctsr=metricMI(grim1,grim2,grimf4);
     qnmicsr=metricMI(grim1,grim2,grimf5);
     qnmisr=metricMI(grim1,grim2,grimf1);
     qnmigff=metricMI(grim1,grim2,grimf3);
     
     qhvsdchwt=metricChenBlum(grim1,grim2,grimf2);
     qhvsmyif=metricChenBlum(grim1,grim2,grimf6);
     qhvsnsctsr=metricChenBlum(grim1,grim2,grimf4);
     qhvscsr=metricChenBlum(grim1,grim2,grimf5);
     qhvssr=metricChenBlum(grim1,grim2,grimf1);
     qhvsgff=metricChenBlum(grim1,grim2,grimf3);
     
     viffdchwt=VIFF_Public(grim1,grim2,grimf2);
     viffmyif=VIFF_Public(grim1,grim2,grimf6);
     viffnsctsr=VIFF_Public(grim1,grim2,grimf4);
     viffcsr=VIFF_Public(grim1,grim2,grimf5);
     viffsr=VIFF_Public(grim1,grim2,grimf1);
     viffgff=VIFF_Public(grim1,grim2,grimf3);
     

 end
 
 save(['myresults/' timeprifix '.mat'])
 
format

fprintf([' qnmisr' ' qnmidchwt' ' qnmigff' ' qnminsctsr' ' qnmicsr' ' qnmimyif','\n']);
disp([qnmisr qnmidchwt  qnmigff qnminsctsr qnmicsr qnmimyif])


fprintf([' qabfsr' ' qabfdchwt' ' qabfgff' ' qabfnsctsr' ' qabfcsr' ' qabfmyif' '\n']);
disp([qabfsr qabfdchwt  qabfgff qabfnsctsr qabfcsr qabfmyif])

fprintf([' qhvssr' ' qhvsdchwt' ' qhvsgff' ' qhvsnsctsr' ' qhvscsr' ' qhvsmyif' '\n']);
disp([qhvssr qhvsdchwt  qhvsgff qhvsnsctsr qhvscsr qhvsmyif])

fprintf([' viffsr' ' viffdchwt' ' viffgff' ' viffnsctsr' ' viffcsr' ' viffmyif' '\n']);
disp([viffsr viffdchwt  viffgff viffnsctsr viffcsr viffmyif])

 
if(showflag)
figure;
subplot(2,3,1),imshow(uint8(im_sr));
subplot(2,3,2),imshow(uint8(imdchwt));
subplot(2,3,3),imshow(uint8(im_gff));
subplot(2,3,4),imshow(uint8(im_nsct_sr));
subplot(2,3,5),imshow(uint8(im_csr));
subplot(2,3,6),imshow(uint8(myif));
img1=uint8(img1);
img2=uint8(img2);
img3=uint8(im_sr);
img4=uint8(imdchwt);
img5=uint8(im_gff);
img6=uint8(im_nsct_sr);
img7=uint8(im_csr);
img8=uint8(myif);
vmar=30;
hmar=60;
marg=255*ones(size(img1,1),30,3);
marg2=255*ones(hmar,size(img1,2)*4+vmar*3,3);
imall=[img1 marg img2 marg img3 marg img4; marg2; img5 marg img6 marg img7 marg img8];
imshow(imall)
% for diving man magnified local fusion 
var=10;har=20;
r1=320;r2=470;
c1=190;c2=390;
marg=255*ones(r2-r1+1,var,3);
marg2=255*ones(har,(c2-c1+1)*4+var*3,3);
magim=[img1(320:470,190:390,:) marg img2(320:470,190:390,:) marg img3(320:470,190:390,:) marg img4(320:470,190:390,:);
    marg2; img5(320:470,190:390,:) marg img6(320:470,190:390,:) marg img7(320:470,190:390,:) marg img8(320:470,190:390,:)];
%magim=[img3(320:470,190:390,:) marg img4(320:470,190:390,:) marg img5(320:470,190:390,:);marg2 ;img6(320:470,190:390,:) marg img7(320:470,190:390,:) marg img8(320:470,190:390,:)];
figure,imshow(magim)
end

