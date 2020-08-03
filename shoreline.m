function  shoreline( h_cd,h_name,SL_name )
%��Ҫȷ���Ĳ���:h0;OA0;L0,dx0,dx;
h0=0;OA0=70;L0=500;dx0=25;dx=5;

%����h
h_name=strcat(h_cd,h_name);
load(h_name);
h=h;


%ȥ��������Χ��NaNֵ
h=h(2:size(h,1)-1,2:size(h,2)-1);

%ȥ���ӵ�L0
idx=ceil(L0/dx0-0.5)+1;
h=h(:,idx:size(h,2));

%��ά��ֵ
% [x0,y0]=meshgrid(1:size(h,2),1:size(h,1));
% [x,y]=meshgrid(1:0.2:size(h,2),1:0.2:size(h,1));  %����5��
% h=interp2(x0,y0,h,x,y,'linear');  %��ά���Բ�ֵ

%����ˮ����������Ϊ½���ˮ��
Lbw=h>h0; % ½��
Wbw=~Lbw;  % ˮ��

%��ˮ�������ͨ����Ϊ�µ�ˮ�򲢼����µ�½��ȥ����½�ذ�Χ��ˮ��
imlabel=bwlabel(Wbw);
stats=regionprops(imlabel,'Area');
area=cat(1,stats.Area);
index=find(area==max(area));
W_out=ismember(imlabel,index);  %�µ�ˮ�����ں�������
ind_W=find(W_out==1);  %ˮ����index
Lbw=~W_out;  %�µ�½�����ں�������
ind_L=find(Lbw==1);  %½����index

%�����Χ½�����С͹�����
[Ly,Lx]=find(Lbw==1);  %½�������꣬��Ϊ������y����Ϊ������x
dt=DelaunayTri(Lx,Ly);
k=convexHull(dt);
Lcx=Lx(k);Lcy=Ly(k);  %��õ�͹����ε�����
ind_Lc=sub2ind(size(Lbw),Lcy,Lcx);  %��͹����ε��index

%����λ��͹��������ˮ����Ϊopen water
[Wy,Wx]=find(W_out==1);  %ˮ��������
ind_in=inpolygon(Wx,Wy,Lcx,Lcy);  %ˮ�������λ��͹������ڵ�index
ind_W_out_in=sub2ind(size(W_out),Wy(ind_in),Wx(ind_in));  %λ��͹�������ˮ����index
ind_W_out_open=sub2ind(size(W_out),Wy(~ind_in),Wx(~ind_in));  %λ��͹��������ˮ��㣨open water����index��

%��ˮ����⣨�ڣ��߽磬��Ϊˮ½����
se=strel('square',3);
% I=imdilate(W_out,se)-W_out;  %ˮ����߽磬ˮ½����㣬���е�����½��
I=W_out-imerode(W_out,se);  %ˮ���ڱ߽磬ˮ½����㣬���е�����ˮ��
[Iy,Ix]=find(I==1);  %
ind_I=find(I==1);  %

%��ͼ��߽��ϵ�½���
Lex1=find(Lx==1);Lex2=find(Lx==size(Lbw,2));  %���ұ߽�
Ley1=find(Ly==1);Ley2=find(Ly==size(Lbw,1));  %���±߽�

%��Q set��T set
Qx=[Wx(ind_in);Ix];Qy=[Wy(ind_in);Iy]; %Q������
ind_Q=sub2ind(size(W_out),Qy,Qx);  %Q���index
Tx=[Ix;Lx(Lex1);Lx(Lex2);Lx(Ley1);Lx(Ley2)];  %T�ĺ����꣨�У�
Ty=[Iy;Ly(Lex1);Ly(Lex2);Ly(Ley1);Ly(Ley2)];  %T�������꣨�У�
ind_T=sub2ind(size(W_out),Ty,Tx);  %T���index

%��T set�ϼ���һ��bank
Tx_bank=-1*ones(size(Lbw,1),1);
Ty_bank=(1:length(Tx_bank))';
Tx=[Tx;Tx_bank];
Ty=[Ty;Ty_bank];

% ��Q set�еĵ��open angle
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

%���ٽ�Ƕȶ�Ӧ��OA��ֵͼ
OAc=OA;
OAc(ind_W_out_open)=180;  %open water��OA����Ϊ180
OAc(OAc>180)=180;  %����180�ĵ㶨Ϊ180
OAbw=OAc;
OAbw=OAbw>=OA0;  %��ֵͼ

%�����OA0����������ͨ��ȥ���߽紦���ܳ��ֵĹ�������
imlabel=bwlabel(OAbw);
stats=regionprops(imlabel,'Area');
area=cat(1,stats.Area);
index=find(area==max(area));
OAbw=ismember(imlabel,index);  %

% %������Ƿ���Ҫ���и�˹ģ��
% [SLs SL_h]=contour(OAbw,[1 1]);  %����ֵͼ�ĵ�ֵ��ͼ
% SLs_num=SLs(2,1);  %��õ�һ����ֵ�ߵĵ���
% if SLs_num~=length(SLs)-1;  %���ݵ�һ����ֵ�ߵ����ǵ������е������ж��Ƿ��ж�����ֵ��
%     G=fspecial('gaussian',[5 5],2);
%     OAbw=imfilter(OAbw,G,'same');  %��˹ģ��
% end
% 
% %���������ͨ��ȥ�����ܳ��ֵĹ�������
% imlabel=bwlabel(OAbw);
% stats=regionprops(imlabel,'Area');
% area=cat(1,stats.Area);
% index=find(area==max(area));
% OAbw=ismember(imlabel,index);

%������ֵ��
[SLs SLs_h]=contour(OAbw,[1 1]);  
SL=contourdata(SLs); %��õ�ֵ����Ϣ�ĺ���

%���������ĵ�ֵ����Ϊ����
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

