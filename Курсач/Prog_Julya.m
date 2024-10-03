clc, clear
%Исходные данные вариант 3
lam_0=[1;0;0;0];
om_0=[3.0e-05;0;0];
t0=0;
tk=50;
Om_0=[0;0;0];
t1=0;
dt=0.01;
i=1;
om=om_0;
Om=Om_0;
Lam=lam_0;
dhdt=[0;0;0];
psi_i=0;
%Данные орбиты вариант 3
r_nach=[-435.428056*10^3;-4043.868745*10^3;5697.150865*10^3]; %Начальное значение вектора r
V_nach=[4.407121*10^3;-6.180556*10^3;-4.050160*10^3]; %Начальное значение вектора V
r=r_nach;
V=V_nach;
Lam_tr = [1; 0; 0; 0];
t1_mass = zeros(701,1);
rx_mass = zeros(701,1); ry_mass = zeros(701,1); rz_mass = zeros(701,1);
Vx_mass = zeros(701,1); Vy_mass = zeros(701,1); Vz_mass = zeros(701,1);
Lam1_mass = zeros(701,1); Lam2_mass = zeros(701,1); Lam3_mass = zeros(701,1); Lam4_mass = zeros(701,1);
omx_mass = zeros(701,1); omy_mass = zeros(701,1); omz_mass = zeros(701,1);
Omx_mass = zeros(701,1); Omy_mass = zeros(701,1); Omz_mass = zeros(701,1);
Lam_tr_mas1 = zeros(701,1); Lam_tr_mas2 = zeros(701,1); Lam_tr_mas3 = zeros(701,1);
Fi_mas = zeros(701,1); Fi_mas_tr = zeros(701,1);
Vozm_mom_mas1 = zeros(701,1); Vozm_mom_mas2 = zeros(701,1); Vozm_mom_mas3 = zeros(701,1);
iter=0;
while (t1 <= 700)
    Vmoment = Ras_vozm_mom(r, Lam);
    if mod(i-1,100)==0
        t1_mass(iter+1) = t1;
        rx_mass(iter+1) = r(1,1); ry_mass(iter+1) = r(2,1); rz_mass(iter+1) = r(3,1);
        Vx_mass(iter+1) = V(1,1); Vy_mass(iter+1) = V(2,1); Vz_mass(iter+1) = V(3,1);
        Lam1_mass(iter+1) = Lam(2,1); Lam2_mass(iter+1) = Lam(3,1); Lam3_mass(iter+1) = Lam(4,1);
        omx_mass(iter+1) = om(1,1); omy_mass(iter+1) = om(2,1); omz_mass(iter+1) = om(3,1);
        Omx_mass(iter+1) = Om(1,1); Omy_mass(iter+1) = Om(2,1); Omz_mass(iter+1) = Om(3,1);
        Lam_tr_mas1(iter+1) = Lam_tr(2,1); Lam_tr_mas2(iter+1) = Lam_tr(3,1); Lam_tr_mas3(iter+1) = Lam_tr(4,1);
        Fi_mas(iter+1) = 2*acos(Lam(1,1)); Fi_mas_tr(iter+1) = 2*acos(Lam_tr(1,1));
        Vozm_mom_mas1(iter+1) = Vmoment(1,1); Vozm_mom_mas2(iter+1) = Vmoment(2,1); Vozm_mom_mas3(iter+1) = Vmoment(3,1);
        iter = iter+1;
    end
    i = i+1;
    [r, V, Lam, om, Om] = Intregr_Runge(t1, dt, r, V, om, Om, Lam, dhdt);
    if (mod(i,10)==0)
        if t1<=tk
            [dhdt, psi_i, Q, Q_3, Lam_tr] = Ras_dhdt(t1, Lam, om, Om, psi_i);
        end
        if t1>tk
            [dhdt] = Ras_dhdt_end(t1, Lam,  om, Om, Lam_tr);
        end
    end
    t1 = t1+dt;
end
Table=table(t1_mass, rx_mass, ry_mass, rz_mass, Vx_mass, Vy_mass, Vz_mass,...
Lam1_mass, Lam2_mass, Lam3_mass, omx_mass, omy_mass, omz_mass,...
Omx_mass, Omy_mass, Omz_mass, Lam_tr_mas1, Lam_tr_mas2, Lam_tr_mas3,...
Fi_mas, Fi_mas_tr)
fileID=fopen('Результаты.txt','w');%открытие файла
writetable(Table,'Результаты.txt')%запись таблицы в файл
fclose(fileID);%закрытие файла

function [a_grav_g] = Ras_agrav(t1, r)
    A = -19089.451590;
    B = 8640184.812866;
    C = 0.093104;
    D = -6.2*10^-6;
    JD0 = 2451545;
    JDD = 36525;
    DS2R = 7.272205216643039903848712*10^-5;
    mu_z = 3.986004418*10^14;%гравитацинный парметр Земли
    J2 = 1.08262668355*10^-3;%коэффициент второй зональной гармоники
    R = 6.378137*10^6;%средний радиус Земли в метрах
    JD=2459761.375+(t1/(24*3600));
    t=(JD-JD0)/JDD;
    f=86400 * mod(JD,1.0);
    alpha=DS2R*((A+(B+(C+D*t)*t)*t)+f);
    alpha=mod(alpha, 2*pi);
    if alpha<0
       alpha=alpha + 2*pi;
    end

    Mgr=[cos(alpha) sin(alpha) 0;
        -sin(alpha) cos(alpha) 0;
         0 0 1];

    rg=Mgr * r;

    per=[(1-5*(r(3,1)/norm(rg))^2)*(r(1,1)/(norm(rg)));
         (1-5*(r(3,1)/norm(rg))^2)*(r(2,1)/(norm(rg)));
         (3-5*(r(3,1)/norm(rg))^2)*(r(3,1)/(norm(rg)))];

    a_grav_g = -1.5*J2*(mu_z/((norm(rg))^2))*((R/(norm(rg)))^2)*per;
    a_grav = Mgr' * a_grav_g;
end

function [a_sun] = Ras_asun(t1, r)
    muS=1.32712440018*10^20;%гравитационный параметр Солнца
    sumOm = 282.940;%сумма долготы восходящего узла и аргумент перицентра
    eps = deg2rad(23.43929111);%наклонение плоскости эклиптики
    JD=2459761.375+(t1/(24*3600));
    T=(JD - 2451545.0)/36525.0;
    M=(357.5226 + 35999.049 * T);
    lamS=deg2rad(sumOm) + deg2rad(M) + deg2rad((6892/3600)) * sin(deg2rad(M)) + deg2rad((72/3600)) * sin(2*deg2rad(M));
    rSun=(149.619 - (2.499 * cos(deg2rad(M))) - 0.021 * cos(2*deg2rad(M))) * 10^6;

    rs = rSun * ([cos(lamS);
                  sin(lamS)*cos(eps);
                  sin(lamS)*sin(eps)])*1000;

    a_sun = muS * ((rs - r)/(norm(rs - r)^3) - rs/(norm(rs)^3));
end

function [Vmoment] = Ras_vozm_mom(r, Lam)
    mu_z = 3.986004418*10^14;
    J=[0.7 0 0;0 1 0;0 0 1];
    r_norm = r*(1/norm(r));
    r0=[0; r_norm(1,1); r_norm(2,1); r_norm(3,1)];

    Lam_=[Lam(1,1); -Lam(2,1); -Lam(3,1); -Lam(4,1)];

    Lam_l=[Lam_(1,1) -Lam_(2,1) -Lam_(3,1) -Lam_(4,1);
           Lam_(2,1) Lam_(1,1) -Lam_(4,1) Lam_(3,1);
           Lam_(3,1) Lam_(4,1) Lam_(1,1) -Lam_(2,1);
           Lam_(4,1) -Lam_(3,1) Lam_(2,1) Lam_(1,1)];

    Lam_r0=Lam_l*r0;

    Lam_l_r0=[Lam_r0(1,1) -Lam_r0(2,1) -Lam_r0(3,1) -Lam_r0(4,1);
              Lam_r0(2,1) Lam_r0(1,1) -Lam_r0(4,1) Lam_r0(3,1);
              Lam_r0(3,1) Lam_r0(4,1) Lam_r0(1,1) -Lam_r0(2,1);
              Lam_r0(4,1) -Lam_r0(3,1) Lam_r0(2,1) Lam_r0(1,1)];

    Lam_l_r0_lam=Lam_l_r0*Lam;
    Eta=[Lam_l_r0_lam(2,1);
         Lam_l_r0_lam(3,1);
         Lam_l_r0_lam(4,1)];
    R=norm(r);
    J_Eta=J*Eta;
    Eta_cross=cross(Eta, J_Eta);

    Vmoment=(3*mu_z/(2*(R^3)))*(Eta_cross);
end

function [dhdt, psi_i, Q, Q_3, Lam_tr] = Ras_dhdt(t1, Lam, om, Om, psi_i)
    dPsi=deg2rad(0.02);
    e=[0;0;1];
    kp=1;
    kd=1.5;
    om_max=629;
    h_max=0.005;
    psi_i=psi_i + dPsi; 
    Lam_tr=[cos(0.5*psi_i);
           e(1,1)*sin(0.5*psi_i);
           e(2,1)*sin(0.5*psi_i);
           e(3,1)*sin(0.5*psi_i)];

    Lam_tr_=[Lam_tr(1,1) Lam_tr(2,1) Lam_tr(3,1) Lam_tr(4,1);
            -Lam_tr(2,1) Lam_tr(1,1) Lam_tr(4,1) -Lam_tr(3,1);
            -Lam_tr(3,1) -Lam_tr(4,1) Lam_tr(1,1) Lam_tr(2,1);
            -Lam_tr(4,1) Lam_tr(3,1) -Lam_tr(2,1) Lam_tr(1,1)];

    Q=Lam_tr_*Lam;
    Q_3=[Q(2,1);Q(3,1);Q(4,1)];
    dhdt=kp*sign(Q(1,1))*Q_3+kd*om;

    i=1;
    while i <= 3
        if abs(-dhdt(i))>=h_max
            dhdt(i)=-h_max*sign(-dhdt(i));
        end
        if abs(Om(i))>=om_max
            dhdt(i)=0;
        end
        i=i+1;
    end
    end

    function [dhdt] = Ras_dhdt_end(t1, Lam,  om, Om, Lam_tr)
    kp=1;
    kd=1.5;
    om_max=629;
    h_max=0.005;

    Lam_tr_=[Lam_tr(1,1) Lam_tr(2,1) Lam_tr(3,1) Lam_tr(4,1);
            -Lam_tr(2,1) Lam_tr(1,1) Lam_tr(4,1) -Lam_tr(3,1);
            -Lam_tr(3,1) -Lam_tr(4,1) Lam_tr(1,1) Lam_tr(2,1);
            -Lam_tr(4,1) Lam_tr(3,1) -Lam_tr(2,1) Lam_tr(1,1)];

    Q=Lam_tr_*Lam;
    Q_3=[Q(2,1);Q(3,1);Q(4,1)];
    dhdt=kp*sign(Q(1,1))*Q_3+kd*om;

    i=1;
    while i <= 3
        if abs(-dhdt(i))>=h_max
            dhdt(i)=-h_max*sign(-dhdt(i));
        end
        if abs(Om(i))>=om_max
            dhdt(i)=0;
        end
        i=i+1;
    end
end

function [r, V, Lam, om, Om] = Intregr_Runge(t1, dt, r, V, om, Om, Lam, dhdt)
    mu_z = 3.986004418*10^14;%гравитацинный парметр Земли
    om_max=629;
    h_max=0.005;
    I=h_max/om_max;
    J=[0.7 0 0;0 1 0;0 0 1];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a_grav=Ras_agrav(t1, r);
    a_sun=Ras_asun(t1, r);
    k11=(-mu_z/((norm(r))^3)*r) +a_sun +a_grav; 
    k12=V;
    Vmoment=Ras_vozm_mom(r, Lam);
    k13=dhdt*(1/I);
    J_I=J*om+I*Om;

    k14=inv(J)*(Vmoment-cross(om, J_I)-dhdt);

    Lamom_0=[0; om(1,1); om(2,1); om(3,1)];
    Lamom=[Lam(1,1) -Lam(2,1) -Lam(3,1) -Lam(4,1);
           Lam(2,1) Lam(1,1) -Lam(4,1) Lam(3,1);
           Lam(3,1) Lam(4,1) Lam(1,1) -Lam(2,1);
           Lam(4,1) -Lam(3,1) Lam(2,1) Lam(1,1)];

    dLam=Lamom*Lamom_0;
    k15=0.5*dLam;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a_grav=Ras_agrav(t1+0.5*dt, r+0.5*dt*k12);
    a_sun=Ras_asun(t1+0.5*dt, r+0.5*dt*k12);
    k21=(-mu_z/((norm(r+0.5*dt*k12))^3)*(r+0.5*dt*k12)) +a_sun +a_grav;  
    k22=V+0.5*dt*k11;
    Vmoment=Ras_vozm_mom(r+0.5*k12, Lam+0.5*dt*k15);
    k23=(1/I)*dhdt;
    J_I=J*(om+0.5*dt*k14)+I*(Om+0.5*dt*k13);
    k24=inv(J)*(Vmoment-cross(om+0.5*k14, J_I)-dhdt);  

    Lamom_0=[0; om(1,1)+0.5*dt*k14(1,1); om(2,1)+0.5*dt*k14(2,1); om(3,1)+0.5*dt*k14(3,1)];
    Lamom=[Lam(1,1)+0.5*dt*k15(1,1) -Lam(2,1)+0.5*dt*k15(2,1) -Lam(3,1)+0.5*dt*k15(3,1) -Lam(4,1)+0.5*dt*k15(4,1);
           Lam(2,1)+0.5*dt*k15(2,1) Lam(1,1)+0.5*dt*k15(1,1) -Lam(4,1)+0.5*dt*k15(4,1) Lam(3,1)+0.5*dt*k15(3,1);
           Lam(3,1)+0.5*dt*k15(3,1) Lam(4,1)+0.5*dt*k15(4,1) Lam(1,1)+0.5*dt*k15(1,1) -Lam(2,1)+0.5*dt*k15(2,1);
           Lam(4,1)+0.5*dt*k15(4,1) -Lam(3,1)+0.5*dt*k15(3,1) Lam(2,1)+0.5*dt*k15(2,1) Lam(1,1)+0.5*dt*k15(1,1)];

    dLam=Lamom*Lamom_0;
    k25=0.5*dLam;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a_grav=Ras_agrav(t1+0.5*dt, r+0.5*dt*k22);
    a_sun=Ras_asun(t1+0.5*dt, r+0.5*dt*k22);
    k31=(-mu_z/((norm(r+0.5*dt*k22))^3)*(r+0.5*dt*k22))+a_sun +a_grav;  
    k32=V+0.5*dt*k21;
    Vmoment=Ras_vozm_mom(r+0.5*k22, Lam+0.5*dt*k25);
    k33=(1/I)*dhdt;
    J_I=J*(om+0.5*dt*k24)+I*(Om+0.5*dt*k23);
    k34=inv(J)*(Vmoment-cross(om+0.5*k24, J_I)-dhdt);  

    Lamom_0=[0; om(1,1)+0.5*dt*k24(1,1); om(2,1)+0.5*dt*k24(2,1); om(3,1)+0.5*dt*k24(3,1)];
    Lamom=[Lam(1,1)+0.5*dt*k25(1,1) -Lam(2,1)+0.5*dt*k25(2,1) -Lam(3,1)+0.5*dt*k25(3,1) -Lam(4,1)+0.5*dt*k25(4,1);
           Lam(2,1)+0.5*dt*k25(2,1) Lam(1,1)+0.5*dt*k25(1,1) -Lam(4,1)+0.5*dt*k25(4,1) Lam(3,1)+0.5*dt*k25(3,1);
           Lam(3,1)+0.5*dt*k25(3,1) Lam(4,1)+0.5*dt*k25(4,1) Lam(1,1)+0.5*dt*k25(1,1) -Lam(2,1)+0.5*dt*k25(2,1);
           Lam(4,1)+0.5*dt*k25(4,1) -Lam(3,1)+0.5*dt*k25(3,1) Lam(2,1)+0.5*dt*k25(2,1) Lam(1,1)+0.5*dt*k25(1,1)];

    dLam=Lamom*Lamom_0;
    k35=0.5*dLam;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    a_grav=Ras_agrav(t1+dt, r+dt*k32);
    a_sun=Ras_asun(t1+dt, r+dt*k32);
    k41=(-mu_z/((norm(r+dt*k32))^3)*(r+dt*k32)) +a_sun +a_grav;  
    k42=V+dt*k31;
    Vmoment=Ras_vozm_mom(r+k32, Lam+dt*k35);
    k43=(1/I)*dhdt;
    J_I=J*(om+dt*k34)+I*(Om+dt*k33);
    k44=inv(J)*(Vmoment-cross(om+k34, J_I)-dhdt);  

    Lamom_0=[0; om(1,1)+dt*k34(1,1); om(2,1)+dt*k34(2,1); om(3,1)+dt*k34(3,1)];
    Lamom=[Lam(1,1)+dt*k35(1,1) -Lam(2,1)+dt*k35(2,1) -Lam(3,1)+dt*k35(3,1) -Lam(4,1)+dt*k35(4,1);
           Lam(2,1)+dt*k35(2,1) Lam(1,1)+dt*k35(1,1) -Lam(4,1)+dt*k35(4,1) Lam(3,1)+dt*k35(3,1);
           Lam(3,1)+dt*k35(3,1) Lam(4,1)+dt*k35(4,1) Lam(1,1)+dt*k35(1,1) -Lam(2,1)+dt*k35(2,1);
           Lam(4,1)+dt*k35(4,1) -Lam(3,1)+dt*k35(3,1) Lam(2,1)+dt*k35(2,1) Lam(1,1)+dt*k35(1,1)];

    dLam=Lamom*Lamom_0;
    k45=0.5*dLam;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V=V+(dt/6)*(k11+2*k21+2*k31+k41);
    r=r+(dt/6)*(k12+2*k22+2*k32+k42);
    Om=Om+(dt/6)*(k13+2*k23+2*k33+k43);
    om=om+(dt/6)*(k14+2*k24+2*k34+k44);
    Lam=Lam+(dt/6)*(k15+2*k25+2*k35+k45);
    sch=1;
    if (Lam(1,1)^2 + Lam(2,1)^2 + Lam(3,1)^2 + Lam(4,1)^2)^0.5~=1
        while sch<=4
            Lam(sch,1)=Lam(sch,1)/((Lam(1,1)^2+Lam(2,1))^2+Lam(3,1)^2+Lam(4,1)^2)^0.5;
            sch=sch+1;
        end
        sch=1;
    end
end
