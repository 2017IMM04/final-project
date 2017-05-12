classdef solid_rand
   properties
   %shorthand numbers for site properties
   %site_p = -1;
   %site_n = 1;
   %dimensions of the solid
   N_x;
   N_y;
   N_z;
   %matrix to save the solid site data -1 for hole, 1 for electron conductor
   data_matrix;
   
   end
   
   methods
       %constructor
       function obj = solid_rand(val_1, val_2, val_3)
            obj.N_x = val_1;
            obj.N_y = val_2;
            obj.N_z = val_3;
            obj.data_matrix = zeros(val_1, val_2, val_3);
            for x = 1:val_1 
                for y = 1:val_2
                    for z = 1:val_3
                        temp = rand();
                        if temp > 0.5
                            obj.data_matrix(x,y,z) = 1; % n type
                        else
                            obj.data_matrix(x,y,z) = -1; % p type
                        end
                    end
                end
            end
       end
       
       %solid site function (takes care of periodic conditions
       function out = site(obj, val_1, val_2, val_3)
           out = obj.data_matrix(mod(val_1-1, obj.N_x)+1, mod(val_2-1, obj.N_y)+1, mod(val_3-1, obj.N_z)+1);
       end
   end      
end