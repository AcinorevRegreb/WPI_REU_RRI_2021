%% This contains our in host function

function dydt = Our_In_Host_func2(t,y,p)

U=y(1);
I=y(2);
V_C=y(3);
F=y(4);
V_P=y(5);
A=y(6);
A_hat=y(7);
T8=y(8);
T8_hat=y(9);
B=y(10);
B_hat=y(11);
T4=y(12);
T4_hat=y(13);
S_P=y(14);
L_P=y(15);
Y=y(16);

dydt=zeros(16,1);

%Organ Compartment
dydt(1) = -(1/7).*p.beta.*U .*V_C;       
dydt(2) = (1/7).*p.beta.*U .*V_C - p.d_C.*I - p.d_I.*I.*T8_hat;
dydt(3) = p.p_v.*(I/(1+(p.eps.*F))) + p.alpha_pc.*V_P - p.alpha_cp.*V_C ...
    - p.d_Y.*V_C.*Y - p.E*V_C - p.avg_lambda_M*V_C;

%IFN
dydt(4) = p.g_F + (2/7).*p.p_F.*(I)-p.d_F.*F;

%Plasma viral load
dydt(5) = -p.alpha_pc.*V_P + p.alpha_cp.*V_C - p.d_Y.*V_P.*Y;

%Naive APC
dydt(6) = -p.g_A.*((A-p.A_gen)/(1+abs(A-p.A_gen))) - (1/7).*p.t_A .* A .*V_C;
%Active APC
dydt(7) = (1/7).*p.t_A .* (A.*V_C) - p.d_A.* A_hat;

%Naive CD8
dydt(8) = -p.g_t8.* ((T8-p.t8_gen)/(1+abs(T8-p.t8_gen))) - p.t_t .* A_hat .*T8;
%Active CD8
dydt(9) = p.t_t.*A_hat.*T8 - p.d_t8.*T8_hat;

%Naive B cells
dydt(10) = -p.g_B.* ((B-p.B_gen)/(1+abs(B-p.B_gen))) - p.t_B .* A_hat .*B;
%Active B Cells
dydt(11) = p.t_B.*A_hat.*B - (p.t_S + p.t_L).*T4_hat.*B_hat;

%Naive CD4
dydt(12) = -p.g_t4.* ((T4-p.t4_gen)/(1+abs(T4-p.t4_gen))) - p.t_t .* A_hat .*T4;
%Activated CD4
dydt(13) = p.t_t.*A_hat.*T4 - p.d_t4.*T4_hat;

%Long and short lived plasma cells
dydt(14) = p.t_S.*T4_hat.*B_hat - p.d_S.*S_P;
dydt(15) = p.t_L.*T4_hat.*B_hat - p.d_L.*L_P;

%Antibodies
dydt(16) = p.g_S .* S_P + p.g_L.*L_P - p.lambda_Y.*Y.*(V_P+V_C)-p.c_Y.*Y;

end