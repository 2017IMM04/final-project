classdef solid_checkered
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
   situ;
   
   end
   
   methods
       %constructor
       function obj = solid_checkered(val_1, val_2, val_3, check_size)
           %check_size means the edge width of the check (unit = lattice)
            obj.N_x = val_1;
            obj.N_y = val_2;
            obj.N_z = val_3;
            obj.data_matrix = zeros(val_1, val_2, val_3);
            obj.situ = zeros(val_1, val_2, val_3);
            for x = 1:val_1 
                for y = 1:val_2
                    if rem((floor((x-1)/check_size)+floor((y-1)/check_size)),2) == 0
                        obj.data_matrix(x,y,:) = 1; % n type
                    else
                        obj.data_matrix(x,y,:) = -1; % p type
                    end
                end
            end
       end
   end
end


           
       
           
       
