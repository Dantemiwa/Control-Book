%���������Ew < 0 etaΪ��������޷�����,������ָ������
%��ʱδ�ؾ��ȶ����ˣ�ֻ�ǻ���־�����罫w0���ó�1��ʱ��ͻ��������״��

clear;
clc
s = tf('s');
Ps = (s - 1) / (s*(s-2));%���崫�ݺ���
Ps_star = ctranspose(Ps);%���Ps(-s)
sigma = 0.1;%���÷���
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
w0 = 0.0126;%0.01;%����ǰ���˲�����������ģ�͵�Ӱ��
cf_flag =1;%�Ƿ�ʹ���˲���
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
%��ƶ����ɶ��ŷ�ϵͳ������һ�����ܶ�Ϊ 1 /s^2������

Ou = 1 / s;
Ou_star = -Ou;
Gu = -Ou * Ou;

%�ж� P(s)�ɿ���,���ɿ���֤������hidden modes,�����������е�assumption
rank_Ps  =  rank(ctrb(ss(Ps)));
if rank_Ps == 0
    error("System Uncontrollable")
else
    disp("System Controllable")
end
%��ȡϵͳ��ĸ�����
A1 = tf(cell2mat(Ps.den),1);%SISO��Ϊ��ĸ
B1 = tf(cell2mat(Ps.num),1);%SISO��Ϊ����

A1_star = ctranspose(A1);%A1*
B1_star = ctranspose(B1);

%��u����������
Au2 = (A1_star * (Ps_star*Ps + k*Q) * A1);
Au = ConjFactor(Au2);%Wiener Hopf�ֽ⣬�������acceptable�Ĵ����Ž�
Au_star = ctranspose(Au);
Gamma = (Au_star^-1 * A1_star * Ps_star * Ou);
Gamma_p = minreal(PartFrac(Gamma));
%�����ı��ʽ(55),�˴�Ϊtheorm 2�Ľ���֮һ������ȡ��ʱ��Eu���ѵ���������
Ru_min = minreal((A1 * Au^-1 * Gamma_p * Ou^-1));%minreal:��Сʵ��
if cf_flag == 1
    Ru = Ru_min *  w0 / (s+w0);
else
    Ru = Ru_min;
end
%��w����������
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
rs = minreal(Ru * us);%minreal���ϣ���֤���ֲ�������
ys = minreal(rs * Ps);
es = (us - ys);
Y = 1;
Rw_min = A1 * Au^-1 * ...
    (PartFrac(Au_star^-1 * I_star * O_star^-1) + PartFracNeg(Au*A1^-1*Y*O))...
    *O^-1*A1;

%���ԣ������ɶȿ�����
%ʹ��˫���ɶȿ�����ʱ���ε�
% Ru = s*(s-2)*((11+4*sqrt(7))*s - 2) / ((s + 2)*(s^2+sqrt(7)*s+1));
% Ru_min = Ru;
% Rw_min = Ru;%��ȡһ�����ɶȽ��жԱ�

%С��ǰ��λһ���ļ����Լ��
%[zw,pw,kw] = tf2zp(cell2mat(Rw_min.num),cell2mat(Rw_min.den));
Rw_min = minreal(CancelPoleZero(Rw_min,2));
% zw1 = round(zw,2);pw1 = round(pw,2);kw1 = round(kw,2);
% Rw_min = minreal(zpk(zw1,pw1,kw1));%wΪ����ıջ�
% Rw_min = tf(Rw_min);
Rw = Rw_min;
%������� Eu Ew Et Es eta��


Ge = (es * ctranspose(es));
%Ge�����ص��ǣ�
Ge = (tf(round(cell2mat(Ge.num),6),round(cell2mat(Ge.den),6)));%���ƾ��ȣ����Լ�С��Է����ֱ���
Gr = (rs * ctranspose(rs));
%ע�⣬����minreal�ܶ඼�ǲ���ȥ���ģ�ȥ�����Ժ󣬷���Сʵ�ֻ���ܻᵼ�»��ֳ���
G_deles = ((1 - Ps * Rw * Ft) *  (Gs * Gr)  * ctranspose((1 - Ps * Rw * Ft)));%ϵͳ���жȼ���
G_deles = (CancelPoleZero(G_deles,3));%������㼫������
G_deles_sym = tf_to_sym(G_deles);
fe = tf_to_sym(Ge);%�������Et
fr =  tf_to_sym(Gr);%�������Es
fu = (1 - Ps*Ru_min) * Gu * ctranspose((1 - Ps*Ru_min)) + Ru_min * Gu * ctranspose(Ru_min);
fu = tf(round(cell2mat(fu.num),6),round(cell2mat(fu.den),6));%���ƾ��ȣ����Լ�С��Է����ֱ���
fu_sym = tf_to_sym(fu);

Eu_min = vpa(1/(2*pi*1i) * int(fu_sym,[-inf*1i,inf*1i]),5)%����ģ�͵���С����Eu
if cf_flag == 0
    f_Zu = 0;
    f_Zu_star = 0;
    Eu = Eu_min;
else
    f_Zu = tf_to_sym(Zu);
    f_Zu_star = tf_to_sym(Zu_star);
    Eu = Eu_min + vpa(1/(1i*2*pi) * int(f_Zu*f_Zu_star,[-inf*1i,inf*1i]),3)%Zu�����Ĵ�����������ʽ75
end

Et = vpa(1/(2*pi*1i) * int(fe,[-inf*1i,inf*1i]) ,3)
Es =  vpa(1/(2*pi*1i) * int(fr,[-inf*1i,inf*1i]) ,3)
Ew = Et + Es - Eu
%vpa(1/(2*pi*i) * int(tf_to_sym((1 - Ps * Ru)*Gu*ctranspose(1-Ps*Ru)),[-inf*i,inf*i]) ,3)
delEt = vpa(1/(2*pi*1i) * int(G_deles_sym,[-inf*1i,inf*1i]) ,3)
eta = sqrt(delEt / Et)


%����Cw and Cu
%С��ǰ��λһ���ļ����Լ��
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
sys = minreal(Cf*sys);%Ҫ����Сʵ�֣�������Сʵ�ֿ��ܻ�����ֵ�ϱ���
figure(1)
step(sys)
figure(2)
nyquist(sys)%���˲�


function tf_G1=sym_to_tf(G1)%���ź���ת��Ϊ���ݺ�����ʽ  
[num,den]=numden(G1);%��ȡ���ű��ʽ���Ӻͷ�ĸ  
Num=sym2poly(num);%���ض���ʽ��ʽϵ��  
Den=sym2poly(den);  
tf_G1 = tf(Num,Den);  
end  

function sym_G1=tf_to_sym(G1)%���ݺ���ת��Ϊ���ź�����ʽ  
syms s  
[num,den]=tfdata(G1);  
Num=poly2sym(num,s);%��ϵ��ϵ��ת��Ϊ����ʽ  
Den=poly2sym(den,s);%  
sym_G1=Num/Den;  
end 

function Gout = ConjFactor(Gin)%�����׹���ֽ�ΪF(s) * F(-s)������ Re s <= 0
    syms s
    assume(s,'positive');
    Gin = minreal(Gin);
    Num = round((cell2mat(Gin.num)),12);%����
    Den = round((cell2mat(Gin.den)),12);
    fun_num = (tf_to_sym(tf(Num,1)));
    fun_den = (tf_to_sym(tf(Den,1)));
    m_num = factor(fun_num,s,'FactorMode','real');%factor�ֽ��������һϵ�л��ʵ�������Ҫ��ȡ����s <= 0�Ĳ���
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
            if (isempty(z) == 0 && z <= 0) %���������
                mult_num = mult_num * G;
            end
        end
        %aur(index,:) = sym2poly(m(i) * m(i + 2));%����ʽ�ֽ����֮��ȫ��һ�εģ���Ҫ�ٺ�����
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
            if (isempty(z) == 0 && z < 0) %���������
                mult_den = mult_den * G;
            elseif (isempty(z) == 0 && z == 0 && zeroFlag == 0)
                zeroFlag = 1;
                mult_den = mult_den * G;
            % ż�����ظ�0�����ʱ��ʹ�ã���ʵӦ�ÿ�ͷ�ж��¸�����Ľ����ȽϺã��˴�͵����
            elseif (isempty(z) == 0 && z == 0 && zeroFlag == 1) 
                zeroFlag = 0;
            end
            
        end
        %aur(index,:) = sym2poly(m(i) * m(i + 2));%����ʽ�ֽ����֮��ȫ��һ�εģ���Ҫ�ٺ�����
        index = index + 1;
    end    
    %�˴�����������жϴ���ȥ�Ĳ����Ƿ�Ϊs > 0
    Gout = mult_num / mult_den;%tf(aur(2,:),1);
end

function Gout = PartFrac(Gin) %���ַ�ʽ���������� Re s <= 0���� ��Ҳ���Ǽ���С��0��
    Gin = minreal(Gin);
    [r,p,k] = residue(cell2mat(Gin.num),cell2mat(Gin.den));%������ r / (s - p) + k
    if isempty(k) == 1
        k = 0;
    end
    G = 0;
    s = tf('s');
    r = round(r,3);%�п�ѧ��������ʾ����ʱ��������vpa��ȡ����
    for i = 1 : length(p)
        if p(i) <= 0
          G = G + r(i) / (s - p(i));    
        end
    end
    Gout = G;
end

function Gout = PartFracNeg(Gin) %���ַ�ʽ���������� Re s <= 0���� ��Ҳ���Ǽ���С��0��
    Gin = minreal(Gin);
    [r,p,k] = residue(cell2mat(Gin.num),cell2mat(Gin.den));%������ r / (s - p) + k
    if isempty(k) == 1
        k = 0;
    end
    G = 0;
    s = tf('s');
    r = round(r,3);%�п�ѧ��������ʾ����ʱ��������vpa��ȡ����
    for i = 1 : length(p)
        if p(i) > 0
          G = G + r(i) / (s - p(i));    
        end
    end
    Gout = G;
end

function Gout = CancelPoleZero(Gin,number) %С������numλ���㼫������
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
