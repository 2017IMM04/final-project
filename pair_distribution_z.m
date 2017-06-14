%calculate pair_distribution parallel to electrode area
%1 for p-type search and -1 for n-type search
function [p_distribution,n_distribution]=pair_distribution_z(type,morphology)
[N_x,N_y,N_z]=size(morphology);
length=[N_x,N_y,N_z];
max_dist=N_z^2;
p_distribution=zeros(max_dist,1);
n_distribution=zeros(max_dist,1);

if type == 1  %(1 for p-type search)
    for x=1:N_x
        for y=1:N_y
            for z=1:N_z
                if morphology(x,y,z) == -1 %-1 for p-type in solid_rand
                    for f=-N_z:N_z
                         sqdist = f^2;
                        if sqdist <= (min(length))^2
                            if morphology(x,y,mod(z+f-1,N_z)+1) == 1  %search for n-type
                               if f == 0
                                   continue
                               else                                    
                                   p_distribution(sqdist)=p_distribution(sqdist)+1; 
                               end
                            elseif morphology(x,y,mod(z+f-1,N_z)+1) == -1  %-1 if morphlolgy is p-type
                               if f == 0
                                   continue
                               else
                                   p_distribution(sqdist)=p_distribution(sqdist)-1;
                               end
                           end
                       else
                            continue
                        end
                    end
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

elseif type == -1 %(-1 for n-type search)
    for x=1:N_x
        for y=1:N_y
            for z=1:N_z
                if morphology(x,y,z) == 1 %1 for n-type in solid_rand
                    for f=-N_z:N_z
                         sqdist = f^2;
                        if sqdist <= (min(length))^2
                            if morphology(x,y,mod(z+f-1,N_z)+1) == -1  %search for p-type
                               if f == 0
                                   continue
                               else                                    
                                   n_distribution(sqdist)=n_distribution(sqdist)+1; 
                               end
                            elseif morphology(x,y,mod(z+f-1,N_z)+1) == 1  %1 if morphlolgy is p-type
                               if f == 0
                                   continue
                               else
                                   n_distribution(sqdist)=n_distribution(sqdist)-1;
                               end
                           end
                       else
                            continue
                        end
                    end
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
    disp('Wrong type description');
end
end