% carsim中的输出端为X坐标（m)、Y坐标(m)、横摆角(°)、质心处的纵向速度（km/h）、刹车、挡位、方向盘转角(°）
function [sys,x0,str,ts] = MPCAPA(t,x,u,flag)
  switch flag %sys输出根据flag的不同而不同
     case 0
      [sys,x0,str,ts] = mdlInitializeSizes; % Initialization 初始化
      
     case 2
      sys = mdlUpdates(t,x,u); % Update discrete states 更新离散状态
      
     case 3
      sys = mdlOutputs(t,x,u); % Calculate outputs 计算输出
    
     case {1,4,9} % Unused flags 滞空
      sys = [];
      
     otherwise
      error(['unhandled flag = ',num2str(flag)]); % Error handling
    end
% End of dsfunc.


%==============================================================
% Initialization
%=============================================================
function [sys,x0,str,ts] = mdlInitializeSizes

% Call simsizes for a sizes structure, fill it in, and convert it 
% to a sizes array.

sizes = simsizes; 
sizes.NumContStates  = 0; % 连续状态量
sizes.NumDiscStates  = 3; % 离散状态量
sizes.NumOutputs     = 2; % 输出量的个数
sizes.NumInputs      = 8; % 输入量的个数
sizes.DirFeedthrough = 1; % Matrix D is non-empty. 在输出函数中直接调用carsim的输入量u(1)~u(3)，使用1。
sizes.NumSampleTimes = 1; % 采样时间的个数
sys = simsizes(sizes); % 设置完后赋值给sys输出
x0 =[0;0;0]; % 状态量初始化，横纵坐标以及横摆角都为0  
global U; % 设定U为全局变量
U=[0;0]; % 初始的控制量为0,0
% Initialize the discrete states.
str = [];             % Set str to an empty matrix.  str为保留参数，mathworks公司还没想好怎么用它，一般在初始化中将其滞空即可
ts  = [0.01 0];       % ts为1*2维度的向量，采样周期+偏移量sample time: [period, offset]，仿真开始0s后每0.01秒运行一次
%End of mdlInitializeSizes
		      
%==============================================================
% Update the discrete states
%==============================================================
function sys = mdlUpdates(t,x,u) % flag=2 表示此时要计算下一个离散状态，即x(k+1)
  
sys = x;
%End of mdlUpdate.

%==============================================================
% Calculate outputs
%==============================================================
function sys = mdlOutputs(t,x,u) % flag=3 表示此时要计算输出，如果sys=[],则表示没有输出
    global a b u_piao; % a,b控制量输出方程中的系数矩阵；u_piao矩阵为全局矩阵,为推导公式中的delt_u_piao(k)
    global U; % U全局控制量误差，2*1 维度矩阵
    global kesi; % kesi为全局，新的状态向量，[k时刻的状态量误差；k-1时刻的控制量误差],5*1矩阵
    tic
    Nx=3;%状态量的个数3
    Nu=2;%控制量的个数2
    Np=60;%预测时域，t+60*0.01 （一般为20~30）
    Nc=30;%控制时域 （一般为预测时域的10%~20%）
    Row=10;%松弛因子权重
    fprintf('Update start, t=%6.3f\n',t)

    % u(3)为车辆横摆角
    t_d =u(3)*pi/180; % t_d为横摆角，CarSim输出的为角度，角度转换为弧度
    
    r(1)=u(4);% 参考路径点x坐标值
    r(2)=u(5);% 参考路径点y坐标值
    r(3)=u(6);% 参考路径点中偏航角yaw值 
    vd1=u(7); % 路径规划信息中的参考速度 （跟踪不同的轨迹有不用的参考速度，规划给的）
    vd2=u(8); % 路径规划信息中的参考方向盘转角 （跟踪不同的轨迹有不用的参考角度，规划给的）

    kesi=zeros(Nx+Nu,1); % 构造新的状态量矩阵 ，5*1维矩阵

    % u(1)为车辆x坐标 u(2)为车辆y坐标 t_d=u(3)为车辆横摆角
    kesi(1)=u(1)-r(1);%u(1)==X(1) 将状态量误差放到相应的位置上，第一行X的误差。kesi(1)为坐标x的误差量。k时刻的实际的x-参考的x
    kesi(2)=u(2)-r(2);%u(2)==X(2) 第二行Y的误差。kesi(2)为坐标y的误差量。k时刻的实际的y-参考的y
    kesi(3)=t_d-r(3); %u(3)==X(3) 第三行是横摆角的误差。kesi(3)为坐标yaw的误差量。k时刻的实际的横摆角-参考的横摆角

    kesi(4)=U(1); % 速度误差放到第五行 我们上一时刻已经算出来的控制量的误差量，k-1时刻的控制量速度
    kesi(5)=U(2); % 前轮偏角误差放到第五行 上一时刻的，k-1时刻的控制量前轮转角
    fprintf('Update start, u(1)=%4.2f\n',U(1))
    fprintf('Update start, u(2)=%4.2f\n',U(2))

    T=0.01; % 采样时间为0.01s，即10毫秒
  % C级轿车 L = 2.77622;
    L = 2.66; % 轴距为2.66m

  % 矩阵初始化   
    u_piao=zeros(Nx,Nu); % 构造3*2维度的矩阵来存放控制量误差
    Q=100*eye(Nx*Np,Nx*Np); % 权重矩阵Q为90*90的单位矩阵   
    R=5*eye(Nu*Nc); % 权重矩阵R为2*2的单位矩阵
    a=[1    0   -vd1*sin(t_d)*T;
       0    1   vd1*cos(t_d)*T;
       0    0   1;]; % a为我们线性离散化后的第一个系数矩阵
    b=[cos(t_d)*T   0;
       sin(t_d)*T   0;
       tan(vd2)*T/L   vd1*T/(cos(vd2)^2);];% b为我们线性离散化后的第二个系数矩阵
    % tan(vd2)*T/L vd1*T/(cos(vd2)^2)
    A_cell=cell(2,2); % 构建2*2的元胞数组
    B_cell=cell(2,1); % 构建2*1的元胞数组
    A_cell{1,1}=a; % 将a矩阵放到A_cell的第一行第一个位置
    A_cell{1,2}=b; % 将b矩阵放到A_cell的第一行第二个位置
    A_cell{2,1}=zeros(Nu,Nx); % 将2*3的零矩阵放到A_cell第二行的第一个位置
    A_cell{2,2}=eye(Nu); % 将2*2的单位阵放到A_cell第二行的第二个位置
    B_cell{1,1}=b; % 将b矩阵放到B_cell的第一行
    B_cell{2,1}=eye(Nu); % 将2*2的单位阵放到B_cell第二行
    A=cell2mat(A_cell); % 这里的A就是我们在推导下一时刻的状态空间时候的A 元胞数组转换成矩阵cell2mat
    B=cell2mat(B_cell); % 这里的B就是我们在推导下一时刻的状态空间时候的B
    C=[1 0 0 0 0;0 1 0 0 0;0 0 1 0 0;];%这个C矩阵是我们输出方程yita的系数矩阵，因为我们可能不需要把每个状态量都输出所以就可以通过设置C矩阵来输出我们想要的状态量，在这里我们输出的是X的误差、Y的误差以及横摆角的误差
    PHI_cell=cell(Np,1);%这个PHI是我们通过总结规律得到的等式右边的第一个系数矩阵，30*1维度
    THETA_cell=cell(Np,Nc);%这里的THETA为我们通过总结规律的到的等式右边的第二个系数矩阵，60*30维度，具体请详见我们CSDN的推导过程
   
    % 通过循环来给第一个系数矩阵PHI赋值
    for j=1:1:Np 
        PHI_cell{j,1}=C*A^j;
        %通过循环来给第二个系数矩阵THETA赋值
        for k=1:1:Nc 
            if k<=j
                THETA_cell{j,k}=C*A^(j-k)*B; % C为3*5矩阵；A为5*5；B为5*2，所以C*A*B为3*2矩阵
            else 
                THETA_cell{j,k}=zeros(Nx,Nu);
            end
        end
    end

    PHI=cell2mat(PHI_cell);%size(PHI)=[Nx*Np Nx+Nu] 180*5维度
    THETA=cell2mat(THETA_cell);%size(THETA)=[Nx*Np Nu*(Nc+1)]

    H_cell=cell(2,2);%这里的H为我们二次规划中的H矩阵，以下来构造二次规划中的H矩阵
    H_cell{1,1}=THETA'*Q*THETA+R;
    H_cell{1,2}=zeros(Nu*Nc,1); %60*1维度
    H_cell{2,1}=zeros(1,Nu*Nc); %1*60维度
    H_cell{2,2}=Row; %H矩阵的右下角的元素就只有一个就是我们的松弛因子
    H=2*cell2mat(H_cell); %由于松弛因子的影响，最终的H矩阵为61*61 %书上这里写错了！应该是2*

    error=PHI*kesi;%这里的error就是我们所设的E矩阵
    f_cell=cell(1,2);%f为二次规划的第二个向量，下面我们来构造它
    f_cell{1,1}=2*error'*Q*THETA;
    f_cell{1,2}=0;
    f=(cell2mat(f_cell))'; %这里的公式为什么书上是f=-cell2mat(f_cell);
    % f=f';
        
 %% 以下为约束生成区域
 %不等式约束
    A_t=zeros(Nc,Nc);%见falcone论文 P181
    for p=1:1:Nc % 第p行
        for q=1:1:Nc
            if q<=p % 列数只有小于或等于行数的时候才赋值为1，构成下三角矩阵
                A_t(p,q)=1;
            else  
                A_t(p,q)=0;
            end
        end 
    end 
    A_I=kron(A_t,eye(Nu));%对应于falcone论文约束处理的矩阵A,求克罗内克积
    Ut=kron(ones(Nc,1),U);%此处感觉论文里的克罗内科积有问题,暂时交换顺序
    
    umin=[-0.1;-0.08];%维数与控制变量的个数相同
    umax=[0.1;0.08];
    delta_umin=[-0.02;-0.002;];%delta_umin=[0.05;-0.0082;];原代码有错，速度变化下界没加负号
    delta_umax=[0.02;0.002];
    Umin=kron(ones(Nc,1),umin);
    Umax=kron(ones(Nc,1),umax);

  % 二次规划不等式约束 Ax<=b                                                         
    A_cons_cell={A_I zeros(Nu*Nc,1);-A_I zeros(Nu*Nc,1)};
    b_cons_cell={Umax-Ut;-Umin+Ut};
    A_cons=cell2mat(A_cons_cell);%（求解方程）状态量不等式约束增益矩阵，转换为绝对值的取值范围
    b_cons=cell2mat(b_cons_cell);%（求解方程）状态量不等式约束的取值
 
  % 二次规划上下限约束 自变量上下限约束 lb<= X <=ub
    M=10; % 松弛因子上限
    delta_Umin=kron(ones(Nc,1),delta_umin);
    delta_Umax=kron(ones(Nc,1),delta_umax);  
    lb=[delta_Umin;0];%（求解方程）状态量下界，包含控制时域内控制增量和松弛因子
    ub=[delta_Umax;M];%（求解方程）状态量上界，包含控制时域内控制增量和松弛因子
    
    %% 开始求解过程
    %options = optimset('Algorithm','active-set'); %新版quadprog不能用有效集法，这里选用内点法 
    options = optimset('Algorithm','interior-point-convex'); 
    %满足 lb ≤ x ≤ ub 的限制条件下求解上述问题。输入 lb 和 ub 是由双精度值组成的向量，这些限制适用于每个 x 分量。如果不存在等式，请设置 Aeq = [] 和 beq = []。
    [X,fval,exitflag]=quadprog(H,f,A_cons,b_cons,[],[],lb,ub,[],options);% 使用 options 中指定的优化选项求解上述问题。使用 optimoptions 创建 options。如果您不提供初始点，请设置 x0 = []。
    %% 计算输出    
    u_piao(1)=X(1);% X为二次规划求出来的控制量序列,取第一行X(1),控制时域内最优的速度控制增量 delt_u_piao(1)
    u_piao(2)=X(2);% X为二次规划求出来的控制量序列,取第二行X(2),最优前轮转角控制增量 delt_u_piao(2)
    U(1)=kesi(4)+u_piao(1);% kesi(4)为速度误差（k-1时刻）,与速度误差增量相加得到我们下一时刻（k时刻）的kesi矩阵的第四行 U(1)=u_piao(k)(1)
    U(2)=kesi(5)+u_piao(2);% kesi(5)为k-1时刻前轮偏角误差，道理同上，与前轮偏角误差的增量相加得到我们下一时刻的kesi矩阵的第五行 U(2)=u_piao(k)(2)
    u_real(1)=U(1)+vd1;% 与参考速度vd1相加才是真正的速度控制量 vd1是ur(k)(1)
    u_real(2)=U(2)+vd2;% 与参考转角vd2相加才是真正的转角控制量 vd2是ur(k)(2)
    sys= u_real;% 真正的控制量输出（2*1矩阵）
    toc % tic与toc结合使用，可以测量运行经过的时间！
% End of mdlOutputs.