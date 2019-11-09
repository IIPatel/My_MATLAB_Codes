clc; clear all; close all;

A=[3 7 13; 1 5 3; 12 3 -5];

i_row=1;
check=0;

while i_row<=3 
   ddc=sum(abs(A(i_row, :)))-abs(A(i_row,i_row));%diagonal dominance check
   
    if A(i_row,i_row)>=ddc
        check=check+1;
        fprintf('%d > %d; check %d; OK\n\n', A(i_row,i_row), ddc, check)
    end
    
    i_row=i_row+1;
end
if check==3
    disp('The matrix is diagonally dominant; Convergent')
else
    disp('The matrix is not diagonally dominant!')
end