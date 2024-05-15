function fixeddof = supportDOFs(sup_all,sup_x,sup_y,sup_z)

fn1 = sup_all; fixeddof1 = reshape([fn1*3-2 fn1*3-1 fn1*3]'...
    ,1,3*length(fn1));                                                      %fixed dofs in all directions

fn2 = sup_x; fixeddof2 = (fn2*3-2)';                                        %fixed dofs in x direction                

fn3 = sup_y; fixeddof3 = (fn3*3-1)';                                        %fixed dofs in y direction 

fn4 = sup_z; fixeddof4 = (fn4*3)';                                          %fixed dofs in z direction 

fixeddof = union(fixeddof1,union(fixeddof2,union(fixeddof3,fixeddof4)));