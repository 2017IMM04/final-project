%calculate pair distribution function
%input type(1 for p-type, -1 for n-type) and morphology matrix
function [p_distribution,n_distribution]=pair_distribution(type,morphology)
[N_x,N_y,N_z]=size(morphology);
max_dist=N_x.^2+N_y.^2+N_z.^2;
p_distribution=zeros(max_dist,1);
n_distribution=zeros(max_dist,1);
if type == 1   
    for x=1:N_x
        for y=1:N_y
            for z=1:N_z
               if morphology(x,y,z) == -1  %(p-type = -1 in solid_rand)
                   for d=0:N_x         %searching range of pair distribution function: x=N_x, y=N_y, z=N_z
                       for e=0:N_y
                           for f=0:N_z
                              sqdist = d.^2 + e.^2 + f.^2;
                               if morphology(mod(x+d-1,N_x)+1,mod(y+e-1,N_y)+1,mod(z+f-1,N_z)+1) == 1  %search for n-type
                                   if x ==d && y == e && z == f
                                       continue
                                   else                                   
                                       p_distribution(sqdist)=p_distribution(sqdist)+1;
                                   end
                               else
                                   continue
                               end
                           end
                       end
                   end    
               else
                   continue
               end
            end
        end
    end
    n_distribution=0;
    for u=1:max_dist        %removing 0 from the discrete distribution
        if p_distribution(u) == 0
            p_distribution(u) = nan;
        else
            continue
        end
    end
    dist=sqrt(1:max_dist)';
    scatter(dist,p_distribution,'filled');
elseif type == -1
    for x=1:N_x
        for y=1:N_y
            for z=1:N_z
               if morphology(x,y,z) == 1  %(n-type = 1 in solid_rand)
                   for d=0:N_x         %searching range of pair distribution function: x=N_x, y=N_y, z=N_z
                       for e=0:N_y
                           for f=0:N_z
                              sqdist = d.^2 + e.^2 + f.^2;
                               if morphology(mod(x+d-1,N_x)+1,mod(y+e-1,N_y)+1,mod(z+f-1,N_z)+1) == -1  %search for p-type
                                   if x ==d && y == e && z == f
                                       continue
                                   else
                                       n_distribution(sqdist)=n_distribution(sqdist)+1;
                                   end
                               else
                                   continue
                               end
                           end
                       end
                   end    
               else
                   continue
               end
            end
        end
    end
    p_distribution=0;
    for u=1:max_dist        %removing 0 from the discrete distribution
        if n_distribution(u) == 0
            n_distribution(u) = nan;
        else
            continue
        end
    end
    dist=sqrt(1:max_dist)';
    scatter(dist,n_distribution,'filled');
else
    disp('Wrong description for morphology type')
end


