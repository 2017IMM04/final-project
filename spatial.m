classdef spatial < handle
    %the handle superclass is used so that we can modify the contents
    %directly (i.e. it is a "pointer", not a copy)
   properties
       N_x;
       N_y;
       N_z;
       work_function_diff;
       dielectric = 3.5*8.8542e-12;
       electron_charge = 1.60217e-19;
       particle_type_matrix;%1
       potential_matrix;%2
       morphology_matrix;%3
   end
   
   methods (Access = public)
       
       %instantiate an instance of this class as follows:
       % this_block = spatial(random_morph, 60, 60, 30, 2.1)
       %for a block with random_morph morphology, 2.1V potential difference
       %function obj = spatial(input_morphology, x_val, y_val, z_val, pot_diff)
       function obj = spatial(x_val, y_val, z_val, pot_diff)
       obj.N_x = x_val;
       obj.N_y = y_val;
       obj.N_z = z_val;
       obj.work_function_diff = pot_diff;
       %obj.morphology_matrix = input_morphology;%morphology is from input
       
       obj.particle_type_matrix = zeros(x_val, y_val, z_val);%particle type initialized to none
       obj.potential_matrix = zeros(x_val, y_val, z_val);
       for x = 1:x_val
           for y = 1:y_val
               for z = 1:z_val
                   obj.potential_matrix(x,y,z)= (z-(1+z_val)/2)/z_val*obj.work_function_diff;
               end
           end
       end
       %potential initialized to linear potential, with z=30 positive
       end
       
       %function that gives update functionality
       function update(obj, type, pos_x, pos_y, pos_z)
           if type ~= obj.particle_type_matrix(pos_x, pos_y, pos_z) %if changed type
               if obj.particle_type_matrix(pos_x, pos_y, pos_z) == 1 || obj.particle_type_matrix(pos_x, pos_y, pos_z) == -1
                   %if it is a particle type that generates an electric field
                   %delete the effect of this particle on the potential
                   for x = -3:3
                       for y = -3:3
                           for z = -3:3 %for each x,y,z neighbor
                               r_square = x^2+y^2+z^2;
                               if r_square <= 9%if satisfies cutoff
                                   original = obj.data_at(2,pos_x,pos_y,pos_z,x,y,z);
                                   if r_square == 0%takes care of self-interaction
                                       final = original - obj.particle_type_matrix(pos_x,pos_y,pos_z)*obj.electron_charge/(4*pi*obj.dielectric*3e-9);
                                   else
                                       final = original - obj.particle_type_matrix(pos_x,pos_y,pos_z)*obj.electron_charge/(4*pi*obj.dielectric*sqrt(r_square)*3e-9);
                                   end
                                   obj.modify_data_at(final, 2,pos_x,pos_y,pos_z,x,y,z);
                               end
                           end
                       end
                   end
               end
               obj.particle_type_matrix(pos_x, pos_y, pos_z) = type;%deletion of effects done, change type
               if type == 1 || type == -1%if type creates potential, update potential
                   for x = -3:3
                       for y = -3:3
                           for z = -3:3 %for each x,y,z neighbor
                               r_square = x^2+y^2+z^2;
                               if r_square <= 9%if satisfies cutoff
                                   original = obj.data_at(2,pos_x,pos_y,pos_z,x,y,z);
                                   if r_square == 0%takes care of self-interaction
                                       final = original + type*obj.electron_charge/(4*pi*obj.dielectric*3e-9);
                                   else
                                       final = original + type*obj.electron_charge/(4*pi*obj.dielectric*sqrt(r_square)*3e-9);
                                   end
                                   obj.modify_data_at(final, 2,pos_x,pos_y,pos_z,x,y,z);
                               end
                           end
                       end
                   end
               end
           end
       end
   
       %function data_at gives convienient way to write in periodic
       %conditions: for example if the data is stored in A, to obtain the potential of sites next to a
       %particle at x,y,z, use A.data_at(2,x,y,z,0,0,1),
       %A.data_at(2,x,y,z,0,0,-1) etc.
       function out = data_at(obj, num, start_x, start_y, start_z, val_1, val_2, val_3)
            if num == 1
                out = obj.particle_type_matrix(mod(start_x+val_1-1, obj.N_x)+1, mod(start_y+val_2-1, obj.N_y)+1, mod(start_z+val_3-1, obj.N_z)+1);
            elseif num == 2
                out = obj.potential_matrix(mod(start_x+val_1-1, obj.N_x)+1, mod(start_y+val_2-1, obj.N_y)+1, mod(start_z+val_3-1, obj.N_z)+1);
            else
                out = obj.morphology_matrix(mod(start_x+val_1-1, obj.N_x)+1, mod(start_y+val_2-1, obj.N_y)+1, mod(start_z+val_3-1, obj.N_z)+1);
            end
       end
   
   end
   
   methods (Access = protected)
       %similar to data_at, but write only, no return value
       function modify_data_at(obj, input, num, start_x, start_y, start_z, val_1, val_2, val_3)
           if num == 1
               obj.particle_type_matrix(mod(start_x+val_1-1, obj.N_x)+1, mod(start_y+val_2-1, obj.N_y)+1, mod(start_z+val_3-1, obj.N_z)+1)=input;
           elseif num == 2
               obj.potential_matrix(mod(start_x+val_1-1, obj.N_x)+1, mod(start_y+val_2-1, obj.N_y)+1, mod(start_z+val_3-1, obj.N_z)+1)=input;
           else
               obj.morphology_matrix(mod(start_x+val_1-1, obj.N_x)+1, mod(start_y+val_2-1, obj.N_y)+1, mod(start_z+val_3-1, obj.N_z)+1)=input;
           end
       end 
    end
end