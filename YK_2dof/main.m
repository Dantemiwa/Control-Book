%如果出现了Ew < 0 eta为无穷或者无法解析,而其他指标正常
%此时未必就稳定不了，只是会出现静差。比如将w0设置成1的时候就会出现这种状况

clear;
clc
s = tf('s');
Ps = (s - 1) / (s*(s-2));%名义传递函数
Ps_star = ctranspose(Ps);%共轭，Ps(-s)
sigma = 0.1;%设置方差
mu = 0;
disp("pertubation of P(s)")
a = normrnd(mu,sigma,2,1)
delPs = (a(1) * s - a(2)) / (s * (s-2));
k = 1;
Q = 1;
Ft = 1;
Gd = 0;
Gm = 0;
miu = 1;
Gs = CancelPoleZero(sigma^2 * (1 - s^2) / (s^2 * (s^2 - 4)),4);
w0 = 0.0126;%0.01;%设置前置滤波器用来降低模型的影响
cf_flag =1;%是否使用滤波器
if cf_flag == 1
    Zu = 1 / (s+w0);
    Zu_star = ctranspose(Zu);
    Zw = 0;
    Zw_star = ctranspose(Zw);
else
    Zu = 0;
    Zu_star = ctranspose(Zu);
    Zw = 0;
    Zw_star = ctranspose(Zw);
end
%设计二自由度伺服系统来跟踪一个谱密度为 1 /s^2的输入

Ou = 1 / s;
Ou_star = -Ou;
Gu = -Ou * Ou;

%判断 P(s)可控性,不可控则证明存在hidden modes,不符合论文中的assumption
rank_Ps  =  rank(ctrb(ss(Ps)));
if rank_Ps == 0
    error("System Uncontrollable")
else
    disp("System Controllable")
end
%获取系统分母与分子
A1 = tf(cell2mat(Ps.den),1);%SISO下为分母
B1 = tf(cell2mat(Ps.num),1);%SISO下为分子

A1_star = ctranspose(A1);%A1*
B1_star = ctranspose(B1);

%求u的最优输入
Au2 = (A1_star * (Ps_star*Ps + k*Q) * A1);
Au = ConjFactor(Au2);%Wiener Hopf分解，用来求次acceptable的次最优解
Au_star = ctranspose(Au);
Gamma = (Au_star^-1 * A1_star * Ps_star * Ou);
Gamma_p = minreal(PartFrac(Gamma));
%接论文表达式(55),此处为theorm 2的结论之一，这样取的时候Eu花费的能量最少
Ru_min = minreal((A1 * Au^-1 * Gamma_p * Ou^-1));%minreal:最小实现
if cf_flag == 1
    Ru = Ru_min *  w0 / (s+w0);
else
    Ru = Ru_min;
end
%求w的最优输入
I = minreal(A1*Ft * (Gd + miu * Gs) * Ps*A1);
I_star = ctranspose(I);
G = minreal(Gm + Ft * (Gd + miu* Gs) * ctranspose(Ft)); 
O = minreal(ConjFactor(A1*G*A1_star));
%Our = 0.01 * (s+1) / (s * (s + 2));
O_star = ctranspose(O);
%Rw_min = minreal(Aur^-1 * PartFrac(Aur_star^-1 * Ir_star * Our_star^-1) * Our^-1)
%Rw_min = minreal(Aur^-1 * PartFrac((Aur_star^-1) * Ir_star * (Our_star^-1)) * Our^-1);
%Y = 1/ (2*B1) ;
us = Ou;
rs = minreal(Ru * us);%minreal加上，保证积分不出问题
ys = minreal(rs * Ps);
es = (us - ys);
Y = 1;
Rw_min = A1 * Au^-1 * ...
    (PartFrac(Au_star^-1 * I_star * O_star^-1) + PartFracNeg(Au*A1^-1*Y*O))...
    *O^-1*A1;

%测试：单自由度控制器
%使用双自由度控制器时屏蔽掉
% Ru = s*(s-2)*((11+4*sqrt(7))*s - 2) / ((s + 2)*(s^2+sqrt(7)*s+1));
% Ru_min = Ru;
% Rw_min = Ru;%削取一个自由度进行对比

%小数前两位一样的极点给约掉
%[zw,pw,kw] = tf2zp(cell2mat(Rw_min.num),cell2mat(Rw_min.den));
Rw_min = minreal(CancelPoleZero(Rw_min,2));
% zw1 = round(zw,2);pw1 = round(pw,2);kw1 = round(kw,2);
% Rw_min = minreal(zpk(zw1,pw1,kw1));%w为输入的闭环
% Rw_min = tf(Rw_min);
Rw = Rw_min;
%计算代价 Eu Ew Et Es eta等


Ge = (es * ctranspose(es));
%Ge这里重点标记，
Ge = (tf(round(cell2mat(Ge.num),6),round(cell2mat(Ge.den),6)));%控制精度，忽略极小项，以防积分崩溃
Gr = (rs * ctranspose(rs));
%注意，以下minreal很多都是不能去掉的，去掉了以后，非最小实现会可能会导致积分出错。
G_deles = ((1 - Ps * Rw * Ft) *  (Gs * Gr)  * ctranspose((1 - Ps * Rw * Ft)));%系统敏感度计算
G_deles = (CancelPoleZero(G_deles,3));%相近的零极点相消
G_deles_sym = tf_to_sym(G_deles);
fe = tf_to_sym(Ge);%计算代价Et
fr =  tf_to_sym(Gr);%计算代价Es
fu = (1 - Ps*Ru_min) * Gu * ctranspose((1 - Ps*Ru_min)) + Ru_min * Gu * ctranspose(Ru_min);
fu = tf(round(cell2mat(fu.num),6),round(cell2mat(fu.den),6));%控制精度，忽略极小项，以防积分崩溃
fu_sym = tf_to_sym(fu);

Eu_min = vpa(1/(2*pi*1i) * int(fu_sym,[-inf*1i,inf*1i]),5)%计算模型的最小代价Eu
if cf_flag == 0
    f_Zu = 0;
    f_Zu_star = 0;
    Eu = Eu_min;
else
    f_Zu = tf_to_sym(Zu);
    f_Zu_star = tf_to_sym(Zu_star);
    Eu = Eu_min + vpa(1/(1i*2*pi) * int(f_Zu*f_Zu_star,[-inf*1i,inf*1i]),3)%Zu带来的代价增量，接式75
end

Et = vpa(1/(2*pi*1i) * int(fe,[-inf*1i,inf*1i]) ,3)
Es =  vpa(1/(2*pi*1i) * int(fr,[-inf*1i,inf*1i]) ,3)
Ew = Et + Es - Eu
%vpa(1/(2*pi*i) * int(tf_to_sym((1 - Ps * Ru)*Gu*ctranspose(1-Ps*Ru)),[-inf*i,inf*i]) ,3)
delEt = vpa(1/(2*pi*1i) * int(G_deles_sym,[-inf*1i,inf*1i]) ,3)
eta = sqrt(delEt / Et)


%计算Cw and Cu
%小数前两位一样的极点给约掉
X = minreal(A1^-1 * (1  - B1*Y)); %
H1 = minreal(A1^-1 * Ru_min);
K1 = minreal((A1^-1 * Rw_min - Y) * A1^-1);
Cu_min = minreal((X - K1*B1)^-1 * (H1));
Cu_min = CancelPoleZero(Cu_min,2);

%Cw = minreal((X - K1*B1)^-1  * - (Y + K1 * A1));
Cw = Rw_min/(1-Rw_min*Ft*Ps);
Cw = CancelPoleZero(Cw,2)
if cf_flag == 1
    Cf =  w0 / (s + w0);
else
    Cf = 1;
end

Cu = Cu_min * Cf
Tc = [Cu  ,-Cw];
Fs = minreal(Cw/Cu_min);
sys = (feedback(Cu_min * minreal(1*delPs + Ps),Fs*Ft,-1));
sys = minreal(Cf*sys);%要用最小实现，不用最小实现可能会在数值上崩掉
figure(1)
step(sys)
figure(2)
nyquist(sys)%无滤波


function tf_G1=sym_to_tf(G1)%符号函数转换为传递函数形式  
[num,den]=numden(G1);%提取符号表达式分子和分母  
Num=sym2poly(num);%返回多项式项式系数  
Den=sym2poly(den);  
tf_G1 = tf(Num,Den);  
end  

function sym_G1=tf_to_sym(G1)%传递函数转换为符号函数形式  
syms s  
[num,den]=tfdata(G1);  
Num=poly2sym(num,s);%把系数系数转换为多项式  
Den=poly2sym(den,s);%  
sym_G1=Num/Den;  
end 

function Gout = ConjFactor(Gin)%功率谱共轭分解为F(s) * F(-s)，返回 Re s <= 0
    syms s
    assume(s,'positive');
    Gin = minreal(Gin);
    Num = round((cell2mat(Gin.num)),12);%精度
    Den = round((cell2mat(Gin.den)),12);
    fun_num = (tf_to_sym(tf(Num,1)));
    fun_den = (tf_to_sym(tf(Den,1)));
    m_num = factor(fun_num,s,'FactorMode','real');%factor分解出来的是一系列互质的因数，要提取其中s <= 0的部分
    m_den = factor(fun_den,s,'FactorMode','real');
    index = 1;
    mult_num = 1;
    mult_den = 1;
    for i = 1 : length(m_num)
        if isSymType(m_num(i),'real')
            if double(m_num(i)) <0
                mult_num = sqrt(-double(m_num(i)));
            end
            
        elseif isa(m_num(i),'sym')
            G = sym_to_tf(m_num(i));
            z = zero(G);
            if (isempty(z) == 0 && z <= 0) %假如是零点
                mult_num = mult_num * G;
            end
        end
        %aur(index,:) = sym2poly(m(i) * m(i + 2));%刚因式分解出来之后全是一次的，需要再合起来
        index = index + 1;
    end
    index = 1;
    zeroFlag = 0;
    for i = 1 : length(m_den)
        if isSymType(m_den(i),'real')
            if double(m_den(i)) <0
                mult_den = sqrt(-double(m_num(i)));
            end
        elseif isa(m_den(i),'sym')
            G = sym_to_tf(m_den(i));
            z = zero(G);
            if (isempty(z) == 0 && z < 0) %假如是零点
                mult_den = mult_den * G;
            elseif (isempty(z) == 0 && z == 0 && zeroFlag == 0)
                zeroFlag = 1;
                mult_den = mult_den * G;
            % 偶数个重复0极点的时候使用，其实应该开头判断下各极点的阶数比较好，此处偷懒了
            elseif (isempty(z) == 0 && z == 0 && zeroFlag == 1) 
                zeroFlag = 0;
            end
            
        end
        %aur(index,:) = sym2poly(m(i) * m(i + 2));%刚因式分解出来之后全是一次的，需要再合起来
        index = index + 1;
    end    
    %此处可用零点来判断传回去的部分是否为s > 0
    Gout = mult_num / mult_den;%tf(aur(2,:),1);
end

function Gout = PartFrac(Gin) %部分分式法，并返回 Re s <= 0部分 （也就是极点小于0）
    Gin = minreal(Gin);
    [r,p,k] = residue(cell2mat(Gin.num),cell2mat(Gin.den));%留数法 r / (s - p) + k
    if isempty(k) == 1
        k = 0;
    end
    G = 0;
    s = tf('s');
    r = round(r,3);%有科学计数法表示的数时，不能用vpa来取精度
    for i = 1 : length(p)
        if p(i) <= 0
          G = G + r(i) / (s - p(i));    
        end
    end
    Gout = G;
end

function Gout = PartFracNeg(Gin) %部分分式法，并返回 Re s <= 0部分 （也就是极点小于0）
    Gin = minreal(Gin);
    [r,p,k] = residue(cell2mat(Gin.num),cell2mat(Gin.den));%留数法 r / (s - p) + k
    if isempty(k) == 1
        k = 0;
    end
    G = 0;
    s = tf('s');
    r = round(r,3);%有科学计数法表示的数时，不能用vpa来取精度
    for i = 1 : length(p)
        if p(i) > 0
          G = G + r(i) / (s - p(i));    
        end
    end
    Gout = G;
end

function Gout = CancelPoleZero(Gin,number) %小数点后第num位的零极点相消
    z = zero(Gin);
    p = pole(Gin);
    G = zpk(Gin);
     Num = cell2mat(Gin.num);
     Den = cell2mat(Gin.den);
     Num = round(Num,6);
     Den = round(Den,6);

    while Num(1) == 0
        Num(1) = [];
    end
    while Den(1) == 0
        Den(1) = [];
    end
    [z,p,k] = tf2zp(Num,Den);
    
    z = round(z,number);
    p = round(p,number);
    k = G.k;
    Gout = tf(minreal(zpk(z,p,k)));
end
