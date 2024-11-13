%Emat function
%
%Given:
%   - size u
%   - entry indices i and j
%
%Returns:
%   - the matrix E of size u x u, with one in its i,j entry and zero elsewhere
%
%24/09/14

function E= Emat(i,j,u)
    E=zeros(u);
    E(i,j)=1;
end