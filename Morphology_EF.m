%{
This is a test code 
This code will generate an "annealed" morphology while under an external
potential
%}

%Dimensions
N_x = 60;
N_y = 60;
N_z = 30;

%Physical constants
k = 1.38e-23;
J = 1/k;
T = 298;
scale_factor = 0.05;

%Steps
N_steps = 1000000;

%initialization of the block
%a = solid(N_x, N_y, N_z);
input_morphology = solid_rand(N_x, N_y, N_z);

for i = 1:N_steps
   %generate two random sets of coordinates
   x_1 = ceil(rand()*N_x);
   y_1 = ceil(rand()*N_y);
   z_1 = ceil(rand()*N_z);
   x_2 = ceil(rand()*N_x);
   y_2 = ceil(rand()*N_y);
   z_2 = ceil(rand()*N_z);
   if input_morphology.site(x_1, y_1, z_1) ~= input_morphology.site(x_2, y_2, z_2)
   %calculate site interaction energy of 1 and 2
   %site_1 = a.site(x_1+1,y_1,z_1)+a.site(x_1-1,y_1,z_1)+a.site(x_1,y_1+1,z_1)+a.site(x_1,y_1-1,z_1)+a.site(x_1,y_1,z_1+1)+a.site(x_1,y_1,z_1-1);
   %site_2 = a.site(x_2+1,y_2,z_2)+a.site(x_2-1,y_2,z_2)+a.site(x_2,y_2+1,z_2)+a.site(x_2,y_2-1,z_2)+a.site(x_2,y_2,z_2+1)+a.site(x_2,y_2,z_2-1);
   site_1 = scale_factor*(z_1-15)+input_morphology.site(x_1+1,y_1,z_1)+input_morphology.site(x_1-1,y_1,z_1)+input_morphology.site(x_1,y_1+1,z_1)+input_morphology.site(x_1,y_1-1,z_1)+input_morphology.site(x_1,y_1,z_1+1)+input_morphology.site(x_1,y_1,z_1-1)+(input_morphology.site(x_1+1,y_1+1,z_1)+input_morphology.site(x_1+1,y_1-1,z_1)+input_morphology.site(x_1-1,y_1+1,z_1)+input_morphology.site(x_1-1,y_1-1,z_1)+input_morphology.site(x_1+1,y_1,z_1+1)+input_morphology.site(x_1+1,y_1,z_1-1)+input_morphology.site(x_1-1,y_1,z_1+1)+input_morphology.site(x_1-1,y_1,z_1-1)+input_morphology.site(x_1,y_1+1,z_1+1)+input_morphology.site(x_1,y_1+1,z_1-1)+input_morphology.site(x_1,y_1-1,z_1+1)+input_morphology.site(x_1,y_1-1,z_1-1))./sqrt(2);
   %site_1 = site_1+(a.site(x_1+1,y_1+1,z_1)+a.site(x_1+1,y_1-1,z_1)+a.site(x_1-1,y_1+1,z_1)+a.site(x_1-1,y_1-1,z_1)+a.site(x_1+1,y_1,z_1+1)+a.site(x_1+1,y_1,z_1-1)+a.site(x_1-1,y_1,z_1+1)+a.site(x_1-1,y_1,z_1-1)+a.site(x_1,y_1+1,z_1+1)+a.site(x_1,y_1+1,z_1-1)+a.site(x_1,y_1-1,z_1+1)+a.site(x_1,y_1-1,z_1-1))./sqrt(2);
   site_2 = scale_factor*(z_2-15)+input_morphology.site(x_2+1,y_2,z_2)+input_morphology.site(x_2-1,y_2,z_2)+input_morphology.site(x_2,y_2+1,z_2)+input_morphology.site(x_2,y_2-1,z_2)+input_morphology.site(x_2,y_2,z_2+1)+input_morphology.site(x_2,y_2,z_2-1)+(input_morphology.site(x_2+1,y_2+1,z_2)+input_morphology.site(x_2+1,y_2-1,z_2)+input_morphology.site(x_2-1,y_2+1,z_2)+input_morphology.site(x_2-1,y_2-1,z_2)+input_morphology.site(x_2+1,y_2,z_2+1)+input_morphology.site(x_2+1,y_2,z_2-1)+input_morphology.site(x_2-1,y_2,z_2+1)+input_morphology.site(x_2-1,y_2,z_2-1)+input_morphology.site(x_2,y_2+1,z_2+1)+input_morphology.site(x_2,y_2+1,z_2-1)+input_morphology.site(x_2,y_2-1,z_2+1)+input_morphology.site(x_2,y_2-1,z_2-1))./sqrt(2);
   %site_2 = site_2+(a.site(x_2+1,y_2+1,z_2)+a.site(x_2+1,y_2-1,z_2)+a.site(x_2-1,y_2+1,z_2)+a.site(x_2-1,y_2-1,z_2)+a.site(x_2+1,y_2,z_2+1)+a.site(x_2+1,y_2,z_2-1)+a.site(x_2-1,y_2,z_2+1)+a.site(x_2-1,y_2,z_2-1)+a.site(x_2,y_2+1,z_2+1)+a.site(x_2,y_2+1,z_2-1)+a.site(x_2,y_2-1,z_2+1)+a.site(x_2,y_2-1,z_2-1))./sqrt(2);
   diff = -J*(input_morphology.site(x_1,y_1,z_1)*site_2+input_morphology.site(x_2,y_2,z_2)*site_1-input_morphology.site(x_1,y_1,z_1)*site_1-input_morphology.site(x_2,y_2,z_2)*site_2);
   if diff <= 0 %if exchange lowers the energy, accept 
       temp = input_morphology.data_matrix(x_1,y_1,z_1);
       input_morphology.data_matrix(x_1,y_1,z_1) = input_morphology.data_matrix(x_2,y_2,z_2);
       input_morphology.data_matrix(x_2,y_2,z_2) = temp;
   else %accept with boltzmann probability 
       temp_rand = rand();
       if exp(-diff) > temp_rand
           temp = input_morphology.data_matrix(x_1,y_1,z_1);
           input_morphology.data_matrix(x_1,y_1,z_1) = input_morphology.data_matrix(x_2,y_2,z_2);
           input_morphology.data_matrix(x_2,y_2,z_2) = temp;
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
            Color=[Color,input_morphology.data_matrix(x,y,z)];
        end
    end
end
scatter3(X,Y,Z,10,Color);

%save('input_file.mat','input_morphology.datamatrix');
