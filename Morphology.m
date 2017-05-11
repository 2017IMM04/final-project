%{
This is a test code 
This code will generate an "annealed" morphology
%}

%Dimensions
N_x = 60;
N_y = 60;
N_z = 30;

%Physical properties
J = 10000;
k = 1.38e-23;
T = 298;

%Steps
N_steps = 10000000;

%initialization of the block
a = solid(N_x, N_y, N_z);

for i = 1:N_steps
   %generate two random sets of coordinates
   x_1 = ceil(rand()*N_x);
   y_1 = ceil(rand()*N_y);
   z_1 = ceil(rand()*N_z);
   x_2 = ceil(rand()*N_x);
   y_2 = ceil(rand()*N_y);
   z_2 = ceil(rand()*N_z);
   if a.site(x_1, y_1, z_1) ~= a.site(x_2, y_2, z_2)
   %calculate site interaction energy of 1 and 2
   site_1 = a.site(x_1+1,y_1,z_1)+a.site(x_1-1,y_1,z_1)+a.site(x_1,y_1+1,z_1)+a.site(x_1,y_1-1,z_1)+a.site(x_1,y_1,z_1+1)+a.site(x_1,y_1,z_1-1);
   site_2 = a.site(x_2+1,y_2,z_2)+a.site(x_2-1,y_2,z_2)+a.site(x_2,y_2+1,z_2)+a.site(x_2,y_2-1,z_2)+a.site(x_2,y_2,z_2+2)+a.site(x_2,y_2,z_2-1);
   diff = -J*(a.site(x_1,y_1,z_1)*site_2+a.site(x_2,y_2,z_2)*site_1-a.site(x_1,y_1,z_1)*site_1-a.site(x_2,y_2,z_2)*site_2);
   if diff <= 0 %if exchange lowers the energy, accept 
       temp = a.data_matrix(x_1,y_1,z_1);
       a.data_matrix(x_1,y_1,z_1) = a.data_matrix(x_2,y_2,z_2);
       a.data_matrix(x_2,y_2,z_2) = temp;
   else %accept with boltzmann probability 
       temp_rand = rand();
       if exp(-diff/k/T) > temp_rand
           temp = a.data_matrix(x_1,y_1,z_1);
           a.data_matrix(x_1,y_1,z_1) = a.data_matrix(x_2,y_2,z_2);
           a.data_matrix(x_2,y_2,z_2) = temp;
       end
   end
   end
end

%Graphics
X=[];
Y=[];
Z=[];
Color = [];
for x = 1:N_x
    for y = 1:N_y
        for z = 1:N_z
            X=[X,x];
            Y=[Y,y];
            Z=[Z,z];
            Color=[Color,a.data_matrix(x,y,z)];
        end
    end
end
scatter3(X,Y,Z,10,Color);
