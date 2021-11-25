% 
clear;

Delta_1=0e-2;

gx=0.00;
gy=0.3;
%  0.2
gz=0.;
mu=0.05;

sigma_x=[0,1;1,0];sigma_y=sqrt(-1)*[0,-1;1,0];sigma_z=[1,0;0,-1];
pair_mat=[0,1;-1,0];
k_all_range=1;

% sigma_xz=kron(sigma_z,sigma_x);sigma_yz=kron(sigma_z,sigma_y);sigma_zz=kron(sigma_z,sigma_z);
% sigma_x0=kron(eye(2),sigma_x);sigma_z0=kron(eye(2),sigma_z);sigma_y0=kron(eye(2),sigma_y);
% pair_mat_1=zeros(4);pair_mat_1(1,2)=1;pair_mat_1(2,1)=-1;
% pair_mat_2=zeros(4);pair_mat_2(3,4)=1;pair_mat_2(4,3)=-1;
% pair_mat_1_tri=diag([1,1,0,0],0);
% pair_mat_2_tri=diag([0,0,1,1],0);
% % to be checked
% pair_mat_inter=zeros(4);pair_mat_inter(1,4)=1;pair_mat_inter(2,3)=-1;pair_mat_inter(3:4,1:2)=-pair_mat_inter(1:2,3:4);
% pair_mat_inter_tri=diag([1,1],2)-diag([1,1],-2);
% % to be checked
% pair_mat_inter=zeros(4);pair_mat_inter(1,4)=1;pair_mat_inter(2,3)=-1;pair_mat_inter(3:4,1:2)=pair_mat_inter(1:2,3:4);
Ham_SC=zeros(4,4);
n_kx=201;n_ky=101;
kx_list=linspace(-k_all_range,k_all_range,n_kx);
ky_list=linspace(-k_all_range,k_all_range,n_ky);


E=zeros(n_kx,n_ky,4);
for i_kx=1:n_kx
    kx=kx_list(i_kx);
    for i_ky=1:n_ky
        ky=ky_list(i_ky);
        Ham_SC(1:2,1:2)=Ham_normal_Top( kx,ky,gx,gy,gz )-mu*eye(2);
        Ham_SC(1:2,3:4)=Delta_1*pair_mat;
        Ham_SC(3:4,3:4)=-conj(Ham_normal_Top( -kx,-ky,gx,gy,gz))+mu*eye(2);
        Ham_SC(3:4,1:2)=Ham_SC(1:2,3:4)';
        E(i_kx,i_ky,:)=eig(Ham_SC);
    end
end


figure
mesh(kx_list,ky_list,E(:,:,2)');hold on;mesh(kx_list,ky_list,E(:,:,3)');
Ene_x=zeros(n_kx);x_x=kx_list;y_x=zeros(n_kx,1);plot3(x_x,y_x,Ene_x,'k');
Ene_y=zeros(n_ky);x_y=zeros(n_ky,1);y_y=ky_list;plot3(x_y,y_y,Ene_y,'k');
% xlabel('$k_x/\sqrt{\frac{v}{\lambda}}$', 'FontName', 'Times New Roman','FontSize',18,'Color','k', 'Interpreter', 'LaTeX');
% ylabel('$k_y/\sqrt{\frac{v}{\lambda}}$', 'FontName', 'Times New Roman','FontSize',18,'Color','k', 'Interpreter', 'LaTeX');zlim([-0.2,0.2]);
% zlabel('$E/\epsilon^*$', 'FontName', 'Times New Roman','FontSize',18,'Color','k', 'Interpreter', 'LaTeX');zlim([-0.2,0.2]);
xlabel('$k_x$', 'FontName', 'Times New Roman','FontSize',18,'Color','k', 'Interpreter', 'LaTeX');
ylabel('$k_y$', 'FontName', 'Times New Roman','FontSize',18,'Color','k', 'Interpreter', 'LaTeX');zlim([-0.2,0.2]);
zlabel('$E$', 'FontName', 'Times New Roman','FontSize',18,'Color','k', 'Interpreter', 'LaTeX');zlim([-0.2,0.2]);

xlim([-1.0,1.0]);ylim([-1.0,1.0]);
%titleNam=sprintf('Direct Gap D_{sin 1}:%.2f,D_{sin 2}:%.2f,g_y:%.2f,g_z: %.3f,mu:%.3f',Delta_1,Delta_2,gy,gz,mu);
%title(titleNam)

% 2-D
Ene_line=zeros(n_kx,4);
for i_kx=1:n_kx
    kx=kx_list(i_kx);
        ky=0;
        Ham_SC(1:2,1:2)=Ham_normal_Top( kx,ky,gx,gy,gz )-mu*eye(2);
        Ham_SC(1:2,3:4)=Delta_1*pair_mat;
        Ham_SC(3:4,3:4)=-conj(Ham_normal_Top( -kx,-ky,gx,gy,gz))+mu*eye(2);
        Ham_SC(3:4,1:2)=Ham_SC(1:2,3:4)';
        Ene_line(i_kx,:)=eig(Ham_SC);
end
figure
plot(x_x,Ene_line,'linewidth',2);hold on;plot(x_x,zeros(numel(x_x),1),'--');
xlabel('$k_x$', 'FontName', 'Times New Roman','FontSize',18,'Color','k', 'Interpreter', 'LaTeX');
ylabel('$E$', 'FontName', 'Times New Roman','FontSize',18,'Color','k', 'Interpreter', 'LaTeX');
% ylim([-0.5,0.5])
figure
plot(x_x,Ene_line,'linewidth',2);hold on;plot(x_x,zeros(numel(x_x),1),'--');
xlabel('$k_x$', 'FontName', 'Times New Roman','FontSize',18,'Color','k', 'Interpreter', 'LaTeX');
ylabel('$E$', 'FontName', 'Times New Roman','FontSize',18,'Color','k', 'Interpreter', 'LaTeX');
ylim([-0.5,0.5])
