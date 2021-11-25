clear;

Delta_1=10e-2;

gx=0.0;
%  gy=-0.5;gz=-0.1;
gy=0.22;gz=-0.;

% gy=0.22;gz=0;
mu=0.2;
dE=6e-3;
k_all_range=0.8;spin_length=3;

sigma_x=[0,1;1,0];sigma_y=sqrt(-1)*[0,-1;1,0];sigma_z=[1,0;0,-1];
pair_mat=[0,-1;1,0];

Ham_SC=zeros(4,4);
n_kx=321;n_ky=321;
kx_list=linspace(-k_all_range,k_all_range,n_kx);
ky_list=linspace(-k_all_range,k_all_range,n_ky);


Ene=zeros(n_kx,n_ky,8);
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
contour(kx_list,ky_list,E(:,:,2)',[0,1]);hold on;
contour(kx_list,ky_list,E(:,:,3)',[0,1]);




%plot figure
index_1=find(abs(E(:,:,2))<dE);
i_ky_FS_1=floor(index_1./n_kx)+1;
i_kx_FS_1=mod(index_1,n_kx);
figure
plot(kx_list(i_kx_FS_1),ky_list(i_ky_FS_1),'.');hold on;

index_2=find(abs(E(:,:,3))<dE);
i_ky_FS_2=floor(index_2./n_kx)+1;
i_kx_FS_2=mod(index_2,n_kx);
plot(kx_list(i_kx_FS_2),ky_list(i_ky_FS_2),'.');
xlim([-k_all_range,k_all_range]);ylim([-k_all_range,k_all_range]);
xlabel('k_x');ylabel('k_y');

i_kx_FS=[i_kx_FS_1;i_kx_FS_2];
i_ky_FS=[i_ky_FS_1;i_ky_FS_2];
kx_FS_list=kx_list(i_kx_FS);
ky_FS_list=ky_list(i_ky_FS);
% plot(kx_FS_list,ky_FS_list,'.');
% xlim([-k_all_range,k_all_range]);ylim([-k_all_range,k_all_range]);
% xlabel('k_x');ylabel('k_y');
n_FS=length(kx_FS_list);
n_cut_FS=length(i_kx_FS_1);
eig_FS=zeros(4,1);
% c_k\dagger \sigma c_k-h_k\dagger sigma^*h_k
sigma_spin_z=diag([1,-1,-1,1],0);
sigma_spin_x=diag([1,0,-1],1);sigma_spin_x=sigma_spin_x+sigma_spin_x';
sigma_spin_y=-sqrt(-1)*diag([1,0,1],1);sigma_spin_y=sigma_spin_y+sigma_spin_y';
sigma_charge_z=-diag([1,1,-1,-1],0);
Spin_Po=zeros(n_FS,3);
PH=zeros(n_FS,1);

for i_FS=1:n_FS
    kx=kx_FS_list(i_FS);
    ky=ky_FS_list(i_FS);
    Ham_SC(1:2,1:2)=Ham_normal_Top( kx,ky,gx,gy,gz )-mu*eye(2);
    Ham_SC(1:2,3:4)=Delta_1*pair_mat;
    Ham_SC(3:4,3:4)=-conj(Ham_normal_Top( -kx,-ky,gx,gy,gz))+mu*eye(2);
    Ham_SC(3:4,1:2)=Ham_SC(1:2,3:4)';
    [Vec,~]=eig(Ham_SC);
    if i_FS<=n_cut_FS
        eig_FS=Vec(:,2);
    else
        eig_FS=Vec(:,3);
    end
    PH(i_FS)=eig_FS'*sigma_charge_z*eig_FS;
    Spin_Po(i_FS,3)=eig_FS'*sigma_spin_z*eig_FS;
    Spin_Po(i_FS,2)=eig_FS'*sigma_spin_y*eig_FS;
    Spin_Po(i_FS,1)=eig_FS'*sigma_spin_x*eig_FS;
end
figure
quiver(kx_FS_list(1:3:n_FS)',ky_FS_list(1:3:n_FS)',0.3*Spin_Po((1:3:n_FS),1),0.3*Spin_Po((1:3:n_FS),2),spin_length);hold on;
% scatter3(kx_FS_list(1:n_FS)',ky_FS_list(1:n_FS)',PH,[],PH,'.');colorbar;colormap(jet);
% xlabel('$k_x/\sqrt{\frac{v_F}{\lambda}}$', 'FontName', 'Times New Roman','FontSize',18,'Color','k', 'Interpreter', 'LaTeX');
% ylabel('$k_y/\sqrt{\frac{v_F}{\lambda}}$', 'FontName', 'Times New Roman', 'FontSize',18,'Color','k', 'Interpreter', 'LaTeX');
scatter3(kx_FS_list(1:n_FS)',ky_FS_list(1:n_FS)',PH,[],PH,'.');colorbar;colormap(jet);
xlabel('$k_x$', 'FontName', 'Times New Roman','FontSize',18,'Color','k', 'Interpreter', 'LaTeX');
ylabel('$k_y$', 'FontName', 'Times New Roman', 'FontSize',18,'Color','k', 'Interpreter', 'LaTeX');

