%This function will produce cross product matrix for a vector A

%vector can be 3x1 or 1x3, does not matter

function[cross_mat] = CrossMatrix( A )

Ax = A(1);
Ay = A(2);
Az = A(3);

cross_mat  = [ 0 -Az Ay;
              Az  0  -Ax;
             -Ay  Ax  0 ];
         
end