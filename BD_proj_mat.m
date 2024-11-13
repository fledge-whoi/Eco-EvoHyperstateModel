% BD_proj_mat function
%
%
%Given
%   - a n x n x u_1 x ... x u_r multi-dimensional array B, 
%
%Returns
%   - the block-diagonal matrix A in which the blocks are given by the 
%     n x n matrix B(:,:,i_1,...,i_r) 
%
% 23/7/15

function A = BD_proj_mat(B)
siz=size(B);
siz=siz(2:end);
s=prod(siz);                    %size of the expected block-diagonal matrix
sk=s/siz(1);                    %number of block on the diagonal
siz=siz(2:end);                 %maximal value of each index i_1,...,i_r
A=zeros(s,s);

for i=1:sk 
 A=A+kron(Emat(i,i,sk),B(:,:,ind2sub(siz,i)));
end
end