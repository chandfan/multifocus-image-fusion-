%% public release. for paper:mulatifocus image fusion
function [scomap,sc1,myll]=xingtaixue(fullscore)
sc1=1-fullscore>.5;
sc1=imfill(sc1,'holes');
% imshow(sc1);
% [h,w]=size(sc1);
% [xx,yy]=meshgrid(1:h,1:w);
% mesh(xx,yy,double(sc1));
% scoreraw=zeros(h*w,3);
% for i=1:h
%     for j=1:w
%         scoreraw(i*j,:)=[i j double(sc1(i,j))];
%         
%         
%     end
% end
% ind=kmeans(scoreraw,2);
% cluscore=reshape(ind,h,w);
% [scx, scy]=find(sc1==1);

sc1=double(sc1);
se = strel('disk',3);        
sc1 = imerode(sc1,se);
% imshow(sc1)


[ll,lnum]=bwlabel(sc1,4);
max=0;
indmax=0;
for k=1:lnum
    [y,x]=find(ll==k);
    nsize=length(y);
    if(nsize>max)
        max=nsize;
        indmax=k;
    end
end
if(indmax==0)
    return
end

myll=ll==indmax;
% figure,imshow(myll);


I2 = imfill(myll,'holes');
% figure,imshow(I2);
scomap=I2;
