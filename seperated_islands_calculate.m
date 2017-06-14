%morphology analyzer for finding each useless cluster 

load('10000000steps_5.mat'); %input morphology data 
number_of_cluster = 0; 
cluster = [] ; 

A=a;

%first delete the useful semiconductor from cathode (z = 30)
for i= 1:A.N_x 
    for j = 1:A.N_y
        if A.data_matrix(i,j,30) == -1 
          A = delector(A,i,j,30,-1);  
        end
        if A.data_matrix(i,j,1) == 1 
          A = delector(A,i,j,1,1);
        end
    end
end

%then delete the useful semiconductor from anode (z = 30)


%analyze the cluster 
for i= 1:A.N_x
    for j = 1:A.N_y 
        for k = 1:A.N_z
            if A.data_matrix(i,j,k) == -1
                number_of_cluster = number_of_cluster +1 ;
                cluster(number_of_cluster , A.data_matrix(i,j,k)+3 ) = size_counter(A,i,j,k,-1) ;
                A = delector(A,i,j,k,-1);
            end
            if A.data_matrix(i,j,k) == 1
                number_of_cluster = number_of_cluster +1 ;
                cluster(number_of_cluster , A.data_matrix(i,j,k)+3 ) = size_counter(A,i,j,k,1) ;
                A = delector(A,i,j,k,1);
            end
        end
    end
end