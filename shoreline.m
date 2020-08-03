function  shoreline( h_cd,h_name,SL_name )
%需要确定的参数:h0;OA0;L0,dx0,dx;
h0=0;OA0=70;L0=500;dx0=25;dx=5;

%导入h
h_name=strcat(h_cd,h_name);
load(h_name);
h=h;


%去掉矩阵外围的NaN值
h=h(2:size(h,1)-1,2:size(h,2)-1);

%去掉河道L0
idx=ceil(L0/dx0-0.5)+1;
h=h(:,idx:size(h,2));

%二维插值
% [x0,y0]=meshgrid(1:size(h,2),1:size(h,1));
% [x,y]=meshgrid(1:0.2:size(h,2),1:0.2:size(h,1));  %加密5倍
% h=interp2(x0,y0,h,x,y,'linear');  %二维线性插值

%根据水深将分析区域分为陆域和水域
Lbw=h>h0; % 陆地
Wbw=~Lbw;  % 水域

%求水域最大连通域，作为新的水域并计算新的陆域（去除被陆地包围的水域）
imlabel=bwlabel(Wbw);
stats=regionprops(imlabel,'Area');
area=cat(1,stats.Area);
index=find(area==max(area));
W_out=ismember(imlabel,index);  %新的水域，用于后续计算
ind_W=find(W_out==1);  %水域点的index
Lbw=~W_out;  %新的陆域，用于后续计算
ind_L=find(Lbw==1);  %陆域点的index

%计算包围陆域的最小凸多边形
[Ly,Lx]=find(Lbw==1);  %陆域点的坐标，行为纵坐标y，列为横坐标x
dt=DelaunayTri(Lx,Ly);
k=convexHull(dt);
Lcx=Lx(k);Lcy=Ly(k);  %求得的凸多边形的坐标
ind_Lc=sub2ind(size(Lbw),Lcy,Lcx);  %将凸多边形点的index

%计算位于凸多边形外的水域，作为open water
[Wy,Wx]=find(W_out==1);  %水域点的坐标
ind_in=inpolygon(Wx,Wy,Lcx,Lcy);  %水域坐标点位于凸多边形内的index
ind_W_out_in=sub2ind(size(W_out),Wy(ind_in),Wx(ind_in));  %位于凸多边形内水域点的index
ind_W_out_open=sub2ind(size(W_out),Wy(~ind_in),Wx(~ind_in));  %位于凸多边形外的水域点（open water）的index，

%求水域的外（内）边界，作为水陆交界
se=strel('square',3);
% I=imdilate(W_out,se)-W_out;  %水域外边界，水陆交界点，所有点属于陆域
I=W_out-imerode(W_out,se);  %水域内边界，水陆交界点，所有点属于水域
[Iy,Ix]=find(I==1);  %
ind_I=find(I==1);  %

%求图像边界上的陆域点
Lex1=find(Lx==1);Lex2=find(Lx==size(Lbw,2));  %左右边界
Ley1=find(Ly==1);Ley2=find(Ly==size(Lbw,1));  %上下边界

%求Q set和T set
Qx=[Wx(ind_in);Ix];Qy=[Wy(ind_in);Iy]; %Q的坐标
ind_Q=sub2ind(size(W_out),Qy,Qx);  %Q点的index
Tx=[Ix;Lx(Lex1);Lx(Lex2);Lx(Ley1);Lx(Ley2)];  %T的横坐标（列）
Ty=[Iy;Ly(Lex1);Ly(Lex2);Ly(Ley1);Ly(Ley2)];  %T的纵坐标（行）
ind_T=sub2ind(size(W_out),Ty,Tx);  %T点的index

%在T set上加上一列bank
Tx_bank=-1*ones(size(Lbw,1),1);
Ty_bank=(1:length(Tx_bank))';
Tx=[Tx;Tx_bank];
Ty=[Ty;Ty_bank];

% 求Q set中的点的open angle
OA=zeros(size(W_out));
for i=1:length(Qx);
    OA_temp=zeros(1,length(Tx));
    for j=1:length(Tx);
        Q0x=Qx(i);Q0y=Qy(i);
        T0x=Tx(j);T0y=Ty(j);
        if Q0x==T0x
            if Q0y<T0y
                OA_temp(j)=90;
            elseif Q0y>T0y
                 OA_temp(j)=270;
            else
                OA_temp(j)=0;
            end 
        elseif Q0y==T0y
                if Q0x<T0x
                    OA_temp(j)=0;
                else
                    OA_temp(j)=180;
                end 
        else
            OA_tan=atan((Q0y-T0y)/(Q0x-T0x))*180/pi;
            if OA_tan>0
                if Q0x<T0x
                    OA_temp(j)=OA_tan;
                else
                     OA_temp(j)=OA_tan+180;
                end
            elseif OA_tan<0
                if Q0x<T0x
                    OA_temp(j)=OA_tan+360;
                else
                     OA_temp(j)=OA_tan+180;
                end
            end
        end
    end
    OAs=sort(OA_temp);
    dOA=zeros(1,length(OAs));
    for k=1:length(OAs)-1
        dOA(k)=OAs(k+1)-OAs(k);
    end
    dOA(length(OAs))=360-OAs(length(OAs))+OAs(1);
    dOA=sort(dOA,'descend');
    OA(Qy(i),Qx(i))=dOA(1)+dOA(2)+dOA(3);
end     

%求临界角度对应的OA二值图
OAc=OA;
OAc(ind_W_out_open)=180;  %open water的OA定义为180
OAc(OAc>180)=180;  %大于180的点定为180
OAbw=OAc;
OAbw=OAbw>=OA0;  %二值图

%求大于OA0区域的最大连通域，去掉边界处可能出现的孤立区域
imlabel=bwlabel(OAbw);
stats=regionprops(imlabel,'Area');
area=cat(1,stats.Area);
index=find(area==max(area));
OAbw=ismember(imlabel,index);  %

% %看情况是否需要进行高斯模糊
% [SLs SL_h]=contour(OAbw,[1 1]);  %做二值图的等值线图
% SLs_num=SLs(2,1);  %获得第一条等值线的点数
% if SLs_num~=length(SLs)-1;  %根据第一条等值线点数是等于所有点数；判断是否有多条等值线
%     G=fspecial('gaussian',[5 5],2);
%     OAbw=imfilter(OAbw,G,'same');  %高斯模糊
% end
% 
% %重求最大连通域，去掉可能出现的孤立区域
% imlabel=bwlabel(OAbw);
% stats=regionprops(imlabel,'Area');
% area=cat(1,stats.Area);
% index=find(area==max(area));
% OAbw=ismember(imlabel,index);

%重作等值线
[SLs SLs_h]=contour(OAbw,[1 1]);  
SL=contourdata(SLs); %获得等值线信息的函数

%将点数最多的等值线作为岸线
SL_num=cat(1,SL.num);
index=find(SL_num==max(SL_num));
SLx=SL(index).xdata;
SLy=SL(index).ydata;
L_SL=0;
for i=1:length(SLx)-1
    dL=sqrt((SLx(i)-SLx(i+1))^2+(SLy(i)-SLy(i+1))^2);
    L_SL=L_SL+dL;
end

    ind_SL=sub2ind(size(OAbw),SLy,SLx);
    P_SL=zeros(size(OAbw));
    P_s=ones(1,length(ind_SL));
    P_SL(ind_SL)=P_s;
    save(SL_name,'SL','ind_SL','P_SL','OAbw');

end

