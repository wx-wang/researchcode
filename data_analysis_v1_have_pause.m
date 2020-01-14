%实验数据处理
%%
% 数据导入
clear all;clc
syms y;
syms t;
%data read
dread = xlsread('C:\Users\xingqy\Desktop\AC-22-TS\0,25seg_pause\AC22TS (7-m-lang).xlsx',2);%attention:which sheet?
pause = 0.25;%rest time
N = dread(:,2);
time = dread(:,3);
Horizontal_stress = dread(:,4);
Horizontal_strain = dread(:,5);
data = [N, time, Horizontal_stress, Horizontal_strain];
%%
%data analysis
len = length(N); %length of data
flag = 1; %flag
flag1 = 0;
flag2 = 1;
flag3 = 1;
flag4 = 0;
for i = 1:len-1
    Nd(i-flag1,flag) = N(i); %flag N 将不同循环数据分开
    timed(i-flag1,flag) = time(i);%将不同循环时间分开
    Horizontal_stress_d(i-flag1,flag) = Horizontal_stress(i); %对应应力
    Horizontal_strain_d(i-flag1,flag) = Horizontal_strain(i); %对应应变
    time_d(i-flag1,flag) = time(i); %对应时间

    if data(i+1) - data(i) ~= 0 %判断N次数，确定是否在当前循环
        if data(i+1) - data(i) == 1 %in a circle
            y11 = Horizontal_stress_d(flag2:i-flag1,flag);
            y12 = Horizontal_strain_d(flag2:i-flag1,flag);
            t1 = time_d(flag2:i-flag1,flag);
            for j = 1: length(t1) %排除pause time，对剩下的点进行拟合
                if abs(t1(length(t1))-t1(j)) > pause
                    flag4 = flag4+1;
                else
                end
            end
            y_stress = y11(1:flag4);
            y_strain = y12(1:flag4);
            t_c = t1(1:flag4);
            
            %计算每个循环下，数据点三角函数拟合参数
            delta_ss(flag3,:) = lsqcurvefit(@(a,t_c)a(1)+a(2)*sin(2*pi*10*t_c+a(3))...
            ,ones(1,3),t_c,y_stress);%三角函数拟合
            delta_ee(flag3,:) = lsqcurvefit(@(a,t_c)a(1)+a(2)*sin(2*pi*10*t_c+a(3))...
            ,ones(1,3),t_c,y_strain);
            delta_ss1(flag3,1) = abs(2* delta_ss(flag3,2)); %确定应力幅值
            delta_eps1(flag3,1) = abs(2* delta_ee(flag3,2));%确定应变幅值
            theta1(flag3,1) = (delta_ee(flag3,3) - delta_ss(flag3,3))*180/pi; %相位角拟合

            

            flag3 = flag3+1;
            flag2 = i-flag1+1;
            flag4 = 0;

        elseif data(i+1) - data(i) > 1
            y11 = Horizontal_stress_d(flag2:i-flag1,flag);
            y12 = Horizontal_strain_d(flag2:i-flag1,flag);
            t1 = time_d(flag2:i-flag1,flag);
            for j = 1: length(t1) %排除pause time，对剩下的点进行拟合
                if abs(t1(length(t1))-t1(j)) > pause
                    flag4 = flag4+1;
                else
                end
            end
            y_stress = y11(1:flag4);
            y_strain = y12(1:flag4);
            t_c = t1(1:flag4);

            delta_ss(flag3,:) = lsqcurvefit(@(a,t_c)a(1)+a(2)*sin(2*pi*10*t_c+a(3))...
            ,ones(1,3),t_c,y_stress);
            delta_ee(flag3,:) = lsqcurvefit(@(a,t_c)a(1)+a(2)*sin(2*pi*10*t_c+a(3))...
            ,ones(1,3),t_c,y_strain);

            delta_ss1(flag3,1) = abs(2* delta_ss(flag3,2));
            delta_eps1(flag3,1) = abs(2* delta_ee(flag3,2)); 
            theta1(flag3,1) = (delta_ee(flag3,3) - delta_ss(flag3,3))*180/pi;
            ini_ss(flag3,1)=delta_ss(flag3,1);
            ini_eps(flag3,1)=delta_ee(flag3,1);
            %对一段时间内循环参数进行平均
            Delta_S(flag,1) = abs(mean(delta_ss1));
            Delta_EPS(flag,1) = abs(mean(delta_eps1));
            Theta(flag,1) = (mean(theta1));
            INI_S(flag,1)=mean(ini_ss);
            INI_E(flag,1)=mean(ini_eps);
            p(flag,:) = [Nd(1,flag),Nd(i-flag1,flag)];
            flag = flag+1;
            flag1 = i;
            flag2 = 1;
            flag3 = 1;
            flag4 = 0;
            delta_ss1 = 0;
            delta_eps1 = 0;

        end

    end

end
%%
% 数据拟合
posion_ratio = 0.15+0.35/(1+exp(3.1849-0.04233*(9/5*20+32)));%泊松比
E = (1+3*posion_ratio)*Delta_S.*(Delta_EPS.^-1);
p1 = p(:,1)+p(:,2);
ER = 0.5* E .* p1;
x1 = p(:,1);
poly = polyfit(x1,ER,6);
y = poly(1)*t^6+poly(2)*t^5+poly(3)*t^4+poly(4)*t^3+poly(5)*t^2+poly(6)*t^1+poly(7);
z = subs(y,t,[0:1:p(flag-1,1)]);
[ER_max,N_max] = max(z);
%%
% 数据导出
filename = '0.25sec_7-m-lang.xlsx';
sheet1 = 'delta_stress_parameter';
sheet2 = 'delta_strain_parameter';
title = {'N_0 [-]', 'N_end [-]', 'Delta S [MPa]', 'Delta Eps [-]', '|E*| [MPa]', 'Energy ratio (|E*|N) [MPa]'};
title1 = {'N_maxER'};

xlswrite(filename, N_max, 'sheet1', 'A2');
xlswrite(filename, title1, 'sheet1', 'A1');
xlswrite(filename, title, 'summary1', 'A1');
xlswrite(filename, Delta_S, 'summary1', 'C2');
xlswrite(filename, Delta_EPS, 'summary1', 'D2');
xlswrite(filename, E, 'summary1', 'E2');
xlswrite(filename, ER, 'summary1', 'F2');
xlswrite(filename, p, 'summary1', 'A2');
xlswrite(filename, INI_S, sheet1, 'A2');
xlswrite(filename, INI_E, sheet2, 'A2');