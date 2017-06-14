function delete_number = size_counter(A,i,j,k,type) 
    A.data_matrix(i,j,k) = 0 ;
    storage = particle.empty ; %只是要想要紀錄位置資訊，偷偷亂用一下XD  
    delete_number = 1;
    exacute_number = 0 ;
    while 1 
        for val_x = -1:1 
            for val_y = -1:1
                for val_z = -1:1
                    if k+val_z < 31 && k+val_z >0 
                        if A.data_matrix(mod(i+val_x-1,A.N_x)+1,mod(j+val_y-1,A.N_y)+1,k+val_z) == type 
                            A.data_matrix(mod(i+val_x-1,A.N_x)+1,mod(j+val_y-1,A.N_y)+1,k+val_z) = 0;
                            storage(delete_number) = particle( 1 , mod(i+val_x-1,A.N_x) + 1 , mod(j+val_y-1,A.N_y) + 1 ,  k+val_z );
                            delete_number = delete_number +1 ;
                        end
                    end
                end
            end              
        end
        
        exacute_number = exacute_number + 1;
        
        if  exacute_number >= delete_number 
            break ; 
        else
            i = storage(exacute_number).pos_x;
            j = storage(exacute_number).pos_y;
            k = storage(exacute_number).pos_z; 
        end
    end
end
