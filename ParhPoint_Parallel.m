    clc
    clear
       
    L = 3.95;  %车长/m
    W = 1.97;  %车宽/m
    l = 2.48;   %轴距/m
    lf = 0.8;    %前悬/m
    lr = 0.67;   %后悬/m
    delta_f = 0.524;  %前轮最大转角/rad
    omiga_f = 0.524;  %前轮最大转角转速/rad/s 
    Rmin = 4.3;  %最小转弯半径/m
    h = 4;      %道路宽度/m
    Lp = 7;   %车位长度/m
    Wp = 2.2;    %车位宽度/m 

    xM4 = 0.8;
    yM4 = 1.1;    
    xM0 = 9;
    yM0 = 3.8;    
    
    % 起始点坐标转换，以终止点为坐标原点
    xM1 = xM0-xM4;%8.2
    yM1 = yM0-yM4;
    
    if xM0>(9.8-((4.3-yM0)*(9.8-8.4)/(4.3-3.3)))
        a0 = yM1/(9.8-((4.3-yM0)*(9.8-8.2)/(4.3-3.3))-xM4);
        b0 = -yM1/(2*pi);
        c0 = 2*pi/(9.8-((4.3-yM0)*(9.8-8.2)/(4.3-3.3))-xM4);      
    else 
        a0 = yM1/xM1;
        b0 = -yM1/(2*pi);
        c0 = 2*pi/xM1;              
    end
    
    k=1;
    
    for xm01 = 2:0.01:2.56
        xmr01(k,1) = xm01;
        ymr01(k,1) = (4.3^2-(xm01-2.554)^2)^(1/2)-3.11;
        phi01(k,1) = -(2*xm01 - 1277/250)/(2*(1849/100 - (xm01 - 1277/500)^2)^(1/2));
        tho01(k,1) = 1/Rmin;
        delta_fr01(k,1) = 0.524;
        omiga_fr01(k,1) = 0.524;
    k = k+1;
    end
    omiga_fr = [];
    tho = [];
    k = 1;
    for xm = 1.2:0.01:xM1
        xmr(k,1) = xm+xM4;
        ym(k,1) = a0*xm+b0*sin(c0*xm);
        ymr(k,1) = ym(k,1)+yM4;
        phi(k,1) = atan(a0+b0*c0*cos(c0*xm));
        d_ym(k,1) = a0+b0*c0*cos(c0*xm);
        d2_ym(k,1) = -b0*c0*c0*sin(c0*xm);
        
        if xmr(k,1) >= (9.8-((4.3-yM0)*(9.8-8.2)/(4.3-3.3)))
            ymr(k,1) = yM0;
            phi(k,1) = 0;  
            d_ym(k,1) = 0;
            d2_ym(k,1) = 0;
        end    
        
        tho(k,1) = d2_ym(k,1)/((1+d_ym(k,1)^2)^(3/2));
        delta_fr(k,1) = atan(2.48*(d2_ym(k,1))/((1+d_ym(k,1)^2)^(3/2)));
        omiga_fr(k,1) = -((62*b0*c0^3*cos(c0*xm))/(25*((a0 + b0*c0*cos(c0*xm))^2 + 1)^(3/2)) + (186*b0^2*c0^4*sin(c0*xm)^2*(a0 + b0*c0*cos(c0*xm)))/(25*((a0 + b0*c0*cos(c0*xm))^2 + 1)^(5/2)))/(((a0 + b0*c0*cos(c0*xm))^2 + 1)^(1/2)*((3844*b0^2*c0^4*sin(c0*xm)^2)/(625*((a0 + b0*c0*cos(c0*xm))^2 + 1)^3) + 1));
        if xmr(k,1) >= (9.8-((4.3-yM0)*(9.8-8.2)/(4.3-3.3)))
            omiga_fr(k,1) = 0;
        end    
        k = k+1;
    end
%% config
% 1
for i =1:length(xmr)
dx_1(i) = xmr(end-i+1)';
dy_1(i) = ymr(end-i+1)';
end
for i=1:length(dx_1)-1
    dyaw_1(i+1)=atan((dy_1(i+1)-dy_1(i))/(dx_1(i+1)-dx_1(i)));
    dtheta_1(i)=dyaw_1(i+1)-dyaw_1(i);
    ds_1(i)=sqrt((dy_1(i+1)-dy_1(i))^2+(dx_1(i+1)-dx_1(i))^2);
    dkappa_1(i)=dtheta_1(i)/ds_1(i);
end
dkappa_1=[dkappa_1,dkappa_1(end)];
dir_1=2*ones(size(dx_1)); 

% 2
% dx_2=[dx_1(end),dx_1(end)];
% dy_2=[dy_1(end),dy_1(end)];
dx_2 = xmr01';
dy_2 = ymr01';
for i=1:length(dx_2)-1
    dyaw_2(i+1)=atan((dy_2(i+1)-dy_2(i))/(dx_2(i+1)-dx_2(i)));
    dtheta_2(i)=dyaw_2(i+1)-dyaw_2(i);
    ds_2(i)=sqrt((dy_2(i+1)-dy_2(i))^2+(dx_2(i+1)-dx_2(i))^2);
    dkappa_2(i)=dtheta_2(i)/ds_2(i);
end
dkappa_2=[dkappa_2,dkappa_2(end)];
dkappa_2(1)=dkappa_1(end);
dir_2=1*ones(size(dx_2)); 

% 3
dx_3=[dx_2(end),dx_2(end)];
dy_3=[dy_2(end),dy_2(end)];
for i=1:length(dx_3)-1
    dyaw_3(i+1)=atan((dy_3(i+1)-dy_3(i))/(dx_3(i+1)-dx_3(i)));
    dtheta_3(i)=dyaw_3(i+1)-dyaw_3(i);
    ds_3(i)=sqrt((dy_3(i+1)-dy_3(i))^2+(dx_3(i+1)-dx_3(i))^2);
    dkappa_3(i)=dtheta_3(i)/ds_3(i); 
end
dkappa_3=[dkappa_3,dkappa_3(end)];
dkappa_3(1)=dkappa_2(end);
dir_3=2*ones(size(dx_3));

mv_x = 1.8;
mv_y = 2.2;

for i=1:length(dx_1)
    dx_1(i) = dx_1(i) - mv_x;
    dy_1(i) = dy_1(i) - mv_y;
end
for i=1:length(dx_2)
    dx_2(i) = dx_2(i) - mv_x;
    dy_2(i) = dy_2(i) - mv_y;
end
for i=1:length(dx_3)
    dx_3(i) = dx_3(i) - mv_x;
    dy_3(i) = dy_3(i) - mv_y;
end

path_x=[dx_1,dx_2,dx_3];
path_y=[dy_1,dy_2,dy_3];

plot(path_x,path_y,'r--');


%%
% 绘图
% 车辆各点路径图
figure
box off
set(0,'defaultfigurecolor','w')
hold on
color01 = [1,0,1];
color02 = [0,1,0];
color03 = [1,1,0];
color04 = [0,1,1];
color05 = [0,0,1];
 plot(xmr,ymr,'color',color05,'linewidth',2)
 plot(xmr01,ymr01,'color',color05,'linewidth',2)
%车位
Car=[-1,1.2,1.2,7,7,14;2.2,2.2,0,0,2.2,2.2];%车子的相对坐标
plot(Car(1,:),Car(2,:),'-r','LineWidth',2)
% plot([1.2,7],[2.2,2.2],'--r','LineWidth',2)
%边界线
plot([-1,13.9],[6.2,6.2],'--r','LineWidth',2)

%范围、线宽、图底颜色
box off
axis([-1 14 -1 7]); 

%字体大小
set(gca,'FontSize',24);

%坐标轴粗细
set(gca,'LineWidth',2)

%轴的名称
xlabel('X(m)');% x轴名称
ylabel('Y(m)'); 

%带箭头的XY轴
annotation('arrow',[0.13 0.13],[0.06 0.99],'LineWidth',2,'HeadStyle','plain','HeadLength',18,'HeadWidth',8);
annotation('arrow',[0.1 0.95],[0.126 0.126],'LineWidth',2,'HeadStyle','plain','HeadLength',18,'HeadWidth',8);

%修改XY轴的刻度
set(gca,'XTickLabel',{0:2:14})
set(gca,'YTickLabel',{'','0','1','2','3','4','5','6','7'})




