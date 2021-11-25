clear;
% coding parameter
nk=111;nE=101;n_mat=4;
n_kx=nk;n_ky=nk;

k_all_range=1.;
dk=2*k_all_range/(nk-1);
kx_list=linspace(-k_all_range,k_all_range,nk);
ky_list=linspace(-k_all_range,k_all_range,nk);

E_min=-0.5;E_max=0.5;
E_list=linspace(E_min,E_max,nE);
eta=1e-2;
 

% real tuning parameter
Delta_1=10e-2;
gx=0.00;
gy=0.15;
%  0.2
gz=0.;
mu=0.2;
sigma_x=[0,1;1,0];sigma_y=sqrt(-1)*[0,-1;1,0];sigma_z=[1,0;0,-1];
pair_mat=[0,1;-1,0];
 
Ham_SC=zeros(n_mat,n_mat);

Ene_list=zeros(n_kx,n_ky,4);
GR=zeros(nE,n_kx,n_ky);
GA=zeros(nE,n_kx,n_ky);
A=zeros(nE,n_kx,n_ky);
for i_E=1:nE
    for i_kx=1:n_kx
    kx=kx_list(i_kx);
        for i_ky=1:n_ky
            ky=ky_list(i_ky);
            
            Ham_SC(1:2,1:2)=Ham_normal_Top(kx,ky,gx,gy,gz)-mu*eye(2);
            Ham_SC(1:2,3:4)=Delta_1*pair_mat;
            Ham_SC(3:4,3:4)=-conj(Ham_normal_Top( -kx,-ky,gx,gy,gz))+mu*eye(2);
            Ham_SC(3:4,1:2)=Ham_SC(1:2,3:4)';
            Ene_list(i_kx,i_ky,:)=eig(Ham_SC);
            GR(i_E,i_kx,i_ky)=trace(1/2/pi*pinv((E_list(i_E)+1j*eta)*eye(n_mat)-Ham_SC));
            GA(i_E,i_kx,i_ky)=trace(1/2/pi*pinv((E_list(i_E)-1j*eta)*eye(n_mat)-Ham_SC));
        end
    end
end
figure(2)
A=1j*(GR-GA);
Aky=reshape(A(:,:,round(nk/2)),[nE,n_kx]);
subplot(1,2,1)
surf(kx_list,E_list,real(Aky/2/pi));shading('interp');view(0,90);colorbar;
 
GRk=trapz(dk,trapz(dk,GR,3),2);
GAk=trapz(dk,trapz(dk,GA,3),2);
Ak=1j*(GRk-GAk);
subplot(1,2,2)
plot(real(Ak)/2/pi,E_list,'.-');
xlabel('DOS');ylabel('E');xlim([0,1])

