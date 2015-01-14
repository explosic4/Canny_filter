%本程序用于实验canny滤波算法
%作者：大指挥官
%程序状态：完结，但是效果不如MATLAB自带的Canny滤波器效果好。

pic=imread('e:\test-pic\car1.jpg');
pic=rgb2gray(pic);
%pic=imresize(pic,0.5);
pic_s=pic;
px=size(pic,1);
py=size(pic,2);
% figure(1);
% imshow(pic);
% title('原图');
%%%%%%%%%%%%%%%%%%
%构建高斯滤波算子
%%%%%%%%%%%%%%%%%%
sigma=1;
Gx=floor(3*sigma+0.5);
x=-1*Gx:1*Gx;
y=x;
[x,y]=meshgrid(x,y);
z=1/(2*pi*sigma^2)*exp(-1*(x.^2+y.^2)/(2*sigma.^2));
    %高斯算子归一化
z_sum=sum(sum(z));
z=z/z_sum;
% figure(2);
% mesh(z);
pic=uint8(conv2(double(pic),z));
pic=pic(Gx+1:px+Gx,Gx+1:py+Gx);
pic=uint8(pic);
figure(3);
imshow(pic);
title(strcat('sigma=',num2str(sigma)));

%%%%%%%%%%%%%%%%%%%%%%%
%计算两个方向上的一阶差分
%%%%%%%%%%%%%%%%%%%%%%%
A=-1/2*[-1 1;-1 1];
B=-1/2*[1 1;-1 -1];
pic_x=conv2(double(pic),A);
pic_x=pic_x(2:px+1,2:py+1);
pic_y=conv2(double(pic),B);
pic_y=pic_y(2:px+1,2:py+1);
% figure(4);
% imshow(pic_x,[]);
% title('x方向差分');
% figure(5);
% imshow(pic_y,[]);
% title('y方向差分');

%%%%%%%%%%%%%%%%%%%%%%%
%计算梯度的幅值与方向
%%%%%%%%%%%%%%%%%%%%%%%
gradamp=sqrt(double(pic_x.^2)+double(pic_y.^2));
gradori=atan(double(pic_y)./double(pic_x));
figure(6);
imshow(gradamp,[]);
title('梯度幅值');

%%%%%%%%%%%%%%%%%%%%%%%
%非极大值抑制
%%%%%%%%%%%%%%%%%%%%%%%
%此处是以梯度图为模板，比较梯度图在梯度方向上的局部最大值，从而细化边缘
picre=zeros(size(gradamp));
picre(1,:)=0;
picre(px,:)=0;
picre(:,1)=0;
picre(:,py)=0;
%g1=0;g2=0;g3=0;g4=0;%这四个变量用于插值计算
t1=0;t2=0;          %这两个变量保存插好的值
h1=0;h2=0;          %插值计算时的两个端点的权重
for i=2:px-1
    for j=2:py-1
        if gradamp(i,j)~=0
            theta=gradori(i,j);
            if theta<0
                theta=theta+pi;
            end
            %  g1    g2
            %        c
            %        g3    g4
            if (theta>=pi/2)&&(theta<pi*3/4)
                h1=tan(theta-pi/2);
                h2=1-h1;
                t1=h1*gradamp(i-1,j-1)+h2*gradamp(i-1,j);
                t2=h1*gradamp(i+1,j+1)+h2*gradamp(i+1,j);
            end
            %  g1
            %  g2    c    g3
            %             g4
            if (theta>=pi*3/4)&&(theta<=pi)
                h1=tan(pi-theta);
                h2=1-h1;
                t1=h1*gradamp(i-1,j-1)+h2*gradamp(i,j-1);
                t2=h1*gradamp(i+1,j+1)+h2*gradamp(i,j+1);
            end
            %        g2    g1
            %        c
            %  g4    g3
            if (theta>=pi/4)&&(theta<pi/2)
                h1=tan(pi/2-theta);
                h2=1-h1;
                t1=h1*gradamp(i-1,j+1)+h2*gradamp(i-1,j);
                t2=h1*gradamp(i+1,j-1)+h2*gradamp(i+1,j);
            end
            %             g1
            %  g3    c    g2
            %  g4
            if (theta>=0)&&(theta<pi/4)
                h1=tan(theta);
                h2=1-h1;
                t1=h1*gradamp(i-1,j+1)+h2*gradamp(i,j+1);
                t2=h1*gradamp(i+1,j-1)+h2*gradamp(i,j-1);
            end
            if (gradamp(i,j)>t1) && (gradamp(i,j)>t2)
                picre(i,j)=128;
            end
        end
    end
end
figure(7);
imshow(picre,[]);

%%%%%%%%%%%%%%%%
%双阈值算法
%%%%%%%%%%%%%%%%
%第一步得到两个阈值
%第二步在picre中，把大于高阈值的，设置为255，把低于低阈值的设为0。
%       这样在picre中，255的值代表真边缘，128的值代表假边缘
%第三步，把所有255的点邻域内值为128的点，设置成255。
%第四步，把孤立的值为128的点删去。
%计算阈值
gradamp(picre==0)=0;
gradmax=max(max(gradamp));
gradall=size(find(gradamp~=0),1);%梯度图中，不为0梯度的数目。
level=0.790;
bizhi=0.4;

t=0:0.1:gradmax;
gradhist=hist(gradamp(:),t);
gradhist(1)=0;
%figure(8);
%stem(t,gradhist);
    %找到一个梯度值，所有大于这个梯度的像素数占梯度图中不为0总像素数的1-High
gradlevel=floor(gradall*level+0.5);
gradsum=0;
for i=1:size(t,2)
    gradsum=gradsum+gradhist(i);
    if(gradsum>gradlevel)
        break;
    end
end
threshold_h=(i-1)/10;
threshold_l=bizhi*threshold_h;
flag=find(gradamp>threshold_h);
picre(flag)=255;
flag1=find(gradamp<threshold_l);
picre(flag1)=0;
picre1=picre;
for i=2:px-1
    for j=2:py-1
        if picre(i,j)==255
            picre(picre(i-1:i+1,j-1:j+1)==128)=255;
        end
    end
end
while picre~=picre1
    for i=2:px-1
        for j=2:py-1
            if picre(i,j)==255
                picre(picre(i-1:i+1,j-1:j+1)==128)=255;
            end
        end
    end
end
picre(picre==128)=0;

figure(8);
imshow(picre,[]);
title('本程序实现滤波');


figure(10);
imshow(edge(pic,'canny'));
title('matlab函数');
