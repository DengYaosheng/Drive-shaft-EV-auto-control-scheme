% carsim�е������ΪX���꣨m)��Y����(m)����ڽ�(��)�����Ĵ��������ٶȣ�km/h����ɲ������λ��������ת��(�㣩
function [sys,x0,str,ts] = MPCAPA(t,x,u,flag)
  switch flag %sys�������flag�Ĳ�ͬ����ͬ
     case 0
      [sys,x0,str,ts] = mdlInitializeSizes; % Initialization ��ʼ��
      
     case 2
      sys = mdlUpdates(t,x,u); % Update discrete states ������ɢ״̬
      
     case 3
      sys = mdlOutputs(t,x,u); % Calculate outputs �������
    
     case {1,4,9} % Unused flags �Ϳ�
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
sizes.NumContStates  = 0; % ����״̬��
sizes.NumDiscStates  = 3; % ��ɢ״̬��
sizes.NumOutputs     = 2; % ������ĸ���
sizes.NumInputs      = 8; % �������ĸ���
sizes.DirFeedthrough = 1; % Matrix D is non-empty. �����������ֱ�ӵ���carsim��������u(1)~u(3)��ʹ��1��
sizes.NumSampleTimes = 1; % ����ʱ��ĸ���
sys = simsizes(sizes); % �������ֵ��sys���
x0 =[0;0;0]; % ״̬����ʼ�������������Լ���ڽǶ�Ϊ0  
global U; % �趨UΪȫ�ֱ���
U=[0;0]; % ��ʼ�Ŀ�����Ϊ0,0
% Initialize the discrete states.
str = [];             % Set str to an empty matrix.  strΪ����������mathworks��˾��û�����ô������һ���ڳ�ʼ���н����Ϳռ���
ts  = [0.01 0];       % tsΪ1*2ά�ȵ���������������+ƫ����sample time: [period, offset]�����濪ʼ0s��ÿ0.01������һ��
%End of mdlInitializeSizes
		      
%==============================================================
% Update the discrete states
%==============================================================
function sys = mdlUpdates(t,x,u) % flag=2 ��ʾ��ʱҪ������һ����ɢ״̬����x(k+1)
  
sys = x;
%End of mdlUpdate.

%==============================================================
% Calculate outputs
%==============================================================
function sys = mdlOutputs(t,x,u) % flag=3 ��ʾ��ʱҪ������������sys=[],���ʾû�����
    global a b u_piao; % a,b��������������е�ϵ������u_piao����Ϊȫ�־���,Ϊ�Ƶ���ʽ�е�delt_u_piao(k)
    global U; % Uȫ�ֿ�������2*1 ά�Ⱦ���
    global kesi; % kesiΪȫ�֣��µ�״̬������[kʱ�̵�״̬����k-1ʱ�̵Ŀ��������],5*1����
    tic
    Nx=3;%״̬���ĸ���3
    Nu=2;%�������ĸ���2
    Np=60;%Ԥ��ʱ��t+60*0.01 ��һ��Ϊ20~30��
    Nc=30;%����ʱ�� ��һ��ΪԤ��ʱ���10%~20%��
    Row=10;%�ɳ�����Ȩ��
    fprintf('Update start, t=%6.3f\n',t)

    % u(3)Ϊ������ڽ�
    t_d =u(3)*pi/180; % t_dΪ��ڽǣ�CarSim�����Ϊ�Ƕȣ��Ƕ�ת��Ϊ����
    
    r(1)=u(4);% �ο�·����x����ֵ
    r(2)=u(5);% �ο�·����y����ֵ
    r(3)=u(6);% �ο�·������ƫ����yawֵ 
    vd1=u(7); % ·���滮��Ϣ�еĲο��ٶ� �����ٲ�ͬ�Ĺ켣�в��õĲο��ٶȣ��滮���ģ�
    vd2=u(8); % ·���滮��Ϣ�еĲο�������ת�� �����ٲ�ͬ�Ĺ켣�в��õĲο��Ƕȣ��滮���ģ�

    kesi=zeros(Nx+Nu,1); % �����µ�״̬������ ��5*1ά����

    % u(1)Ϊ����x���� u(2)Ϊ����y���� t_d=u(3)Ϊ������ڽ�
    kesi(1)=u(1)-r(1);%u(1)==X(1) ��״̬�����ŵ���Ӧ��λ���ϣ���һ��X����kesi(1)Ϊ����x���������kʱ�̵�ʵ�ʵ�x-�ο���x
    kesi(2)=u(2)-r(2);%u(2)==X(2) �ڶ���Y����kesi(2)Ϊ����y���������kʱ�̵�ʵ�ʵ�y-�ο���y
    kesi(3)=t_d-r(3); %u(3)==X(3) �������Ǻ�ڽǵ���kesi(3)Ϊ����yaw���������kʱ�̵�ʵ�ʵĺ�ڽ�-�ο��ĺ�ڽ�

    kesi(4)=U(1); % �ٶ����ŵ������� ������һʱ���Ѿ�������Ŀ��������������k-1ʱ�̵Ŀ������ٶ�
    kesi(5)=U(2); % ǰ��ƫ�����ŵ������� ��һʱ�̵ģ�k-1ʱ�̵Ŀ�����ǰ��ת��
    fprintf('Update start, u(1)=%4.2f\n',U(1))
    fprintf('Update start, u(2)=%4.2f\n',U(2))

    T=0.01; % ����ʱ��Ϊ0.01s����10����
  % C���γ� L = 2.77622;
    L = 2.66; % ���Ϊ2.66m

  % �����ʼ��   
    u_piao=zeros(Nx,Nu); % ����3*2ά�ȵľ�������ſ��������
    Q=100*eye(Nx*Np,Nx*Np); % Ȩ�ؾ���QΪ90*90�ĵ�λ����   
    R=5*eye(Nu*Nc); % Ȩ�ؾ���RΪ2*2�ĵ�λ����
    a=[1    0   -vd1*sin(t_d)*T;
       0    1   vd1*cos(t_d)*T;
       0    0   1;]; % aΪ����������ɢ����ĵ�һ��ϵ������
    b=[cos(t_d)*T   0;
       sin(t_d)*T   0;
       tan(vd2)*T/L   vd1*T/(cos(vd2)^2);];% bΪ����������ɢ����ĵڶ���ϵ������
    % tan(vd2)*T/L vd1*T/(cos(vd2)^2)
    A_cell=cell(2,2); % ����2*2��Ԫ������
    B_cell=cell(2,1); % ����2*1��Ԫ������
    A_cell{1,1}=a; % ��a����ŵ�A_cell�ĵ�һ�е�һ��λ��
    A_cell{1,2}=b; % ��b����ŵ�A_cell�ĵ�һ�еڶ���λ��
    A_cell{2,1}=zeros(Nu,Nx); % ��2*3�������ŵ�A_cell�ڶ��еĵ�һ��λ��
    A_cell{2,2}=eye(Nu); % ��2*2�ĵ�λ��ŵ�A_cell�ڶ��еĵڶ���λ��
    B_cell{1,1}=b; % ��b����ŵ�B_cell�ĵ�һ��
    B_cell{2,1}=eye(Nu); % ��2*2�ĵ�λ��ŵ�B_cell�ڶ���
    A=cell2mat(A_cell); % �����A�����������Ƶ���һʱ�̵�״̬�ռ�ʱ���A Ԫ������ת���ɾ���cell2mat
    B=cell2mat(B_cell); % �����B�����������Ƶ���һʱ�̵�״̬�ռ�ʱ���B
    C=[1 0 0 0 0;0 1 0 0 0;0 0 1 0 0;];%���C�����������������yita��ϵ��������Ϊ���ǿ��ܲ���Ҫ��ÿ��״̬����������ԾͿ���ͨ������C���������������Ҫ��״̬���������������������X����Y������Լ���ڽǵ����
    PHI_cell=cell(Np,1);%���PHI������ͨ���ܽ���ɵõ��ĵ�ʽ�ұߵĵ�һ��ϵ������30*1ά��
    THETA_cell=cell(Np,Nc);%�����THETAΪ����ͨ���ܽ���ɵĵ��ĵ�ʽ�ұߵĵڶ���ϵ������60*30ά�ȣ��������������CSDN���Ƶ�����
   
    % ͨ��ѭ��������һ��ϵ������PHI��ֵ
    for j=1:1:Np 
        PHI_cell{j,1}=C*A^j;
        %ͨ��ѭ�������ڶ���ϵ������THETA��ֵ
        for k=1:1:Nc 
            if k<=j
                THETA_cell{j,k}=C*A^(j-k)*B; % CΪ3*5����AΪ5*5��BΪ5*2������C*A*BΪ3*2����
            else 
                THETA_cell{j,k}=zeros(Nx,Nu);
            end
        end
    end

    PHI=cell2mat(PHI_cell);%size(PHI)=[Nx*Np Nx+Nu] 180*5ά��
    THETA=cell2mat(THETA_cell);%size(THETA)=[Nx*Np Nu*(Nc+1)]

    H_cell=cell(2,2);%�����HΪ���Ƕ��ι滮�е�H����������������ι滮�е�H����
    H_cell{1,1}=THETA'*Q*THETA+R;
    H_cell{1,2}=zeros(Nu*Nc,1); %60*1ά��
    H_cell{2,1}=zeros(1,Nu*Nc); %1*60ά��
    H_cell{2,2}=Row; %H��������½ǵ�Ԫ�ؾ�ֻ��һ���������ǵ��ɳ�����
    H=2*cell2mat(H_cell); %�����ɳ����ӵ�Ӱ�죬���յ�H����Ϊ61*61 %��������д���ˣ�Ӧ����2*

    error=PHI*kesi;%�����error�������������E����
    f_cell=cell(1,2);%fΪ���ι滮�ĵڶ�������������������������
    f_cell{1,1}=2*error'*Q*THETA;
    f_cell{1,2}=0;
    f=(cell2mat(f_cell))'; %����Ĺ�ʽΪʲô������f=-cell2mat(f_cell);
    % f=f';
        
 %% ����ΪԼ����������
 %����ʽԼ��
    A_t=zeros(Nc,Nc);%��falcone���� P181
    for p=1:1:Nc % ��p��
        for q=1:1:Nc
            if q<=p % ����ֻ��С�ڻ����������ʱ��Ÿ�ֵΪ1�����������Ǿ���
                A_t(p,q)=1;
            else  
                A_t(p,q)=0;
            end
        end 
    end 
    A_I=kron(A_t,eye(Nu));%��Ӧ��falcone����Լ������ľ���A,������ڿ˻�
    Ut=kron(ones(Nc,1),U);%�˴��о�������Ŀ����ڿƻ�������,��ʱ����˳��
    
    umin=[-0.1;-0.08];%ά������Ʊ����ĸ�����ͬ
    umax=[0.1;0.08];
    delta_umin=[-0.02;-0.002;];%delta_umin=[0.05;-0.0082;];ԭ�����д��ٶȱ仯�½�û�Ӹ���
    delta_umax=[0.02;0.002];
    Umin=kron(ones(Nc,1),umin);
    Umax=kron(ones(Nc,1),umax);

  % ���ι滮����ʽԼ�� Ax<=b                                                         
    A_cons_cell={A_I zeros(Nu*Nc,1);-A_I zeros(Nu*Nc,1)};
    b_cons_cell={Umax-Ut;-Umin+Ut};
    A_cons=cell2mat(A_cons_cell);%����ⷽ�̣�״̬������ʽԼ���������ת��Ϊ����ֵ��ȡֵ��Χ
    b_cons=cell2mat(b_cons_cell);%����ⷽ�̣�״̬������ʽԼ����ȡֵ
 
  % ���ι滮������Լ�� �Ա���������Լ�� lb<= X <=ub
    M=10; % �ɳ���������
    delta_Umin=kron(ones(Nc,1),delta_umin);
    delta_Umax=kron(ones(Nc,1),delta_umax);  
    lb=[delta_Umin;0];%����ⷽ�̣�״̬���½磬��������ʱ���ڿ����������ɳ�����
    ub=[delta_Umax;M];%����ⷽ�̣�״̬���Ͻ磬��������ʱ���ڿ����������ɳ�����
    
    %% ��ʼ������
    %options = optimset('Algorithm','active-set'); %�°�quadprog��������Ч����������ѡ���ڵ㷨 
    options = optimset('Algorithm','interior-point-convex'); 
    %���� lb�0�2�܁0�2x�0�2�܁0�2ub ����������������������⡣���� lb �� ub ����˫����ֵ��ɵ���������Щ����������ÿ�� x ��������������ڵ�ʽ�������� Aeq�0�2=�0�2[] �� beq�0�2=�0�2[]��
    [X,fval,exitflag]=quadprog(H,f,A_cons,b_cons,[],[],lb,ub,[],options);% ʹ�� options ��ָ�����Ż�ѡ������������⡣ʹ�� optimoptions ���� options����������ṩ��ʼ�㣬������ x0�0�2=�0�2[]��
    %% �������    
    u_piao(1)=X(1);% XΪ���ι滮������Ŀ���������,ȡ��һ��X(1),����ʱ�������ŵ��ٶȿ������� delt_u_piao(1)
    u_piao(2)=X(2);% XΪ���ι滮������Ŀ���������,ȡ�ڶ���X(2),����ǰ��ת�ǿ������� delt_u_piao(2)
    U(1)=kesi(4)+u_piao(1);% kesi(4)Ϊ�ٶ���k-1ʱ�̣�,���ٶ����������ӵõ�������һʱ�̣�kʱ�̣���kesi����ĵ����� U(1)=u_piao(k)(1)
    U(2)=kesi(5)+u_piao(2);% kesi(5)Ϊk-1ʱ��ǰ��ƫ��������ͬ�ϣ���ǰ��ƫ������������ӵõ�������һʱ�̵�kesi����ĵ����� U(2)=u_piao(k)(2)
    u_real(1)=U(1)+vd1;% ��ο��ٶ�vd1��Ӳ����������ٶȿ����� vd1��ur(k)(1)
    u_real(2)=U(2)+vd2;% ��ο�ת��vd2��Ӳ���������ת�ǿ����� vd2��ur(k)(2)
    sys= u_real;% �����Ŀ����������2*1����
    toc % tic��toc���ʹ�ã����Բ������о�����ʱ�䣡
% End of mdlOutputs.