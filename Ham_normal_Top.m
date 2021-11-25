function [ Ham ] = Ham_normal_Top( kx,ky,gx,gy,gz)
%HAM_NORMAL 此处显示有关此函数的摘要
%   此处显示详细说明
sigma_x=[0,1;1,0];sigma_y=sqrt(-1)*[0,-1;1,0];sigma_z=[1,0;0,-1];
% sigma_xz=kron(sigma_z,sigma_x);sigma_yz=kron(sigma_z,sigma_y);sigma_zz=kron(sigma_z,sigma_z);
% sigma_x0=kron(eye(2),sigma_x);sigma_z0=kron(eye(2),sigma_z);sigma_y0=kron(eye(2),sigma_y);
k_plus=kx+sqrt(-1)*ky;k_minus=kx-sqrt(-1)*ky;
Ham=-(kx*sigma_y-ky*sigma_x)+0*(kx*kx*eye(2)+ky*ky*eye(2))+(gz)*sigma_z+gx*sigma_x+gy*sigma_y;
%1/2*(k_plus^3+k_minus^3)
end

