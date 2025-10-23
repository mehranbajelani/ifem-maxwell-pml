% Example: customize plot appearance in showrate script

N1 = 1./h1; N2 = 1./h2; N3 =1./h3;
if (nargin<=2)
    k1 = 1; opt1 = '-*'; % line style
end

r1 = -showrate(N1,err1,k1,opt1);
hold on
r2 = -showrate(N2,err2,k2,'--sr'); % red dashed line with square markers
r3 = -showrate(N3,err3,k3,':dg');  % green dotted line with diamond markers

% Update title with LaTeX interpreter for nicer font
title(['\bf Rate of convergence'], 'FontSize',16, 'Interpreter','latex');

% Update legend with custom labels
h_leg = legend(str1,'$C_1 h^1$','$C_2 h^2$','$C_3 h^2$', ...
    'Location','southwest');
set(h_leg,'FontSize',12,'Interpreter','latex');

% Axis labels
xlabel('$\\log(1/h)$','FontSize',14,'Interpreter','latex');
ylabel('Error','FontSize',14,'Interpreter','latex');

% Grid and background
set(gca,'FontSize',12);
grid on