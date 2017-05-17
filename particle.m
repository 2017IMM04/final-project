classdef particle
    properties
        type;
        %type according to charges 1 for hole, 2 for exciton, -1 for
        %electron, 0 for nothing
        pos_x;
        pos_y;
        pos_z;
        
    end
    
    methods
        %constructor, 
        %for example, define a particle by using:
        % A = particle(1, 1, 2, 3) for a hole at coordinates(1,2,3)   
        function obj = particle(type_val, x_val, y_val, z_val)
            obj.type = type_val;
            obj.pos_x = x_val;
            obj.pos_y = y_val;
            obj.pos_z = z_val;
        end
    end
        
end