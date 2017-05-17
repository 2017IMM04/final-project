classdef spatial
    
   properties
       N_x;
       N_y;
       N_z;
       particle_type_matrix;%1
       potential_matrix;%2
       morphology_matrix;%3
   end
   
   methods
       function obj = spatial(input_morphology, x_val, y_val, z_val)
       obj.N_x = x_val;
       obj.N_y = y_val;
       obj.N_z = z_val;
       obj.morphology_matrix = input_morphology;%morphology is from input
       obj.potential_matrix = zeros(x_val, y_val, z_val);%potential initialized to zero
       obj.particle_type_matrix = zeros(x_val, y_val, z_val);%particle type initialized to none
       end
       
       %function data_at gives convienient way to write in periodic
       %conditions: for example if the data is stored in A, to obtain the potential of sites next to a
       %particle at x,y,z, use A.data_at(2,x,y,z,0,0,1),
       %A.data_at(2,x,y,z,0,0,-1) etc.
       function out = data_at(obj, num, start_x, start_y, start_z, val_1, val_2, val_3)
           if num == 1
               out = obj.particle_matrix(mod(start_x+val_1-1, obj.N_x)+1, mod(start_y+val_2-1, obj.N_y)+1, mod(start_z+val_3-1, obj.N_z)+1);
           elseif num == 2
               out = obj.potential_matrix(mod(start_x+val_1-1, obj.N_x)+1, mod(start_y+val_2-1, obj.N_y)+1, mod(start_z+val_3-1, obj.N_z)+1);
           else
               out = obj.morphology_matrix(mod(start_x+val_1-1, obj.N_x)+1, mod(start_y+val_2-1, obj.N_y)+1, mod(start_z+val_3-1, obj.N_z)+1);
           end
           
       end
       
   end
    
    
    
end