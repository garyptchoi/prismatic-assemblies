% Construct a minimum rigidifying link pattern (MRP) for a 2x2x3 rectangular prismatic assembly
%
% Reference:
% G. P. T. Choi, S. Chen, L. Mahadevan, 
% ``Control of connectivity and rigidity in prismatic assemblies.''
% Preprint, arXiv:2005.05845, 2020.

L = 2; %length
M = 2; %width
N = 3; %height
ncube = L*M*N; %Number of cubes
nlink = 2*L*M*N-2 ;% theoretical lower bound for number of links
mat=zeros(ncube*18+nlink*3,ncube*24); 
% column: x1 y1 z1 x2 y2 z2 ...

% Order of the 8 vertices
%   8 - - 7
%  /     /|
% 5 - - 6 |
% | 4 - | 3
% |/    |/
% 1 - - 2

% length constraints:
% 12 quad boundary constraints, and 6 diagonal constraints
for i=0:ncube-1
    % edges
    %x1 x2
    mat(i*18+1,i*24+1)=-1;
    mat(i*18+1,i*24+4)=1;
    %x3 x4
    mat(i*18+2,i*24+7)=-1;
    mat(i*18+2,i*24+10)=1;
    %x5 x6
    mat(i*18+3,i*24+13)=-1;
    mat(i*18+3,i*24+16)=1;
    %x7 x8
    mat(i*18+4,i*24+19)=-1;
    mat(i*18+4,i*24+22)=1;
    
    %y1 y4
    mat(i*18+5,i*24+2)=-1;
    mat(i*18+5,i*24+11)=1;
    %y2 y3
    mat(i*18+6,i*24+5)=-1;
    mat(i*18+6,i*24+8)=1;
    %y5 y8
    mat(i*18+7,i*24+14)=-1;
    mat(i*18+7,i*24+23)=1;
    %y6 y7
    mat(i*18+8,i*24+17)=-1;
    mat(i*18+8,i*24+20)=1;
    
    %z1 z5
    mat(i*18+9,i*24+3)=-1;
    mat(i*18+9,i*24+15)=1;
    %z2 z6
    mat(i*18+10,i*24+6)=-1;
    mat(i*18+10,i*24+18)=1;
    %z3 z7
    mat(i*18+11,i*24+9)=-1;
    mat(i*18+11,i*24+21)=1;
    %z4 z8
    mat(i*18+12,i*24+12)=-1;
    mat(i*18+12,i*24+24)=1;
    
    % face diagonals
    %x1 y1 x3 y3
    mat(i*18+13,i*24+1)=-1;
    mat(i*18+13,i*24+2)=-1;
    mat(i*18+13,i*24+7)=1;
    mat(i*18+13,i*24+8)=1;
    %x5 y5 x7 y7
    mat(i*18+14,i*24+13)=-1;
    mat(i*18+14,i*24+14)=-1;
    mat(i*18+14,i*24+19)=1;
    mat(i*18+14,i*24+20)=1;
    %x1 z1 x6 z6
    mat(i*18+15,i*24+1)=-1;
    mat(i*18+15,i*24+3)=-1;
    mat(i*18+15,i*24+16)=1;
    mat(i*18+15,i*24+18)=1;
    %x4 z4 x7 z7
    mat(i*18+16,i*24+10)=-1;
    mat(i*18+16,i*24+12)=-1;
    mat(i*18+16,i*24+19)=1;
    mat(i*18+16,i*24+21)=1;
    %y1 z1 y8 z8
    mat(i*18+17,i*24+2)=-1;
    mat(i*18+17,i*24+3)=-1;
    mat(i*18+17,i*24+23)=1;
    mat(i*18+17,i*24+24)=1;
    %y2 z2 y7 z7
    mat(i*18+18,i*24+5)=-1;
    mat(i*18+18,i*24+6)=-1;
    mat(i*18+18,i*24+20)=1;
    mat(i*18+18,i*24+21)=1;
    
end

% Order of the cubes
%    [10]- [11]
%   /     / |
% [8] - [9] |
%  |     |  |
%  | [6] - [7]
%  |/    |/ |
% [4] - [5] |
%  |     |  |
%  | [2] - [3]
%  |/    |/
% [0] - [1]

rown=ncube*18+1;
linkpairs=[
    8*0+2, 8*1+1;
    8*0+4, 8*2+1;
    8*1+3, 8*3+2;
    8*2+3, 8*3+4;
    
    8*8+6, 8*9+5;
    8*8+8, 8*10+5;
    8*9+7, 8*11+6;
    8*10+7, 8*11+8;
    
    8*1+6, 8*5+2;
    8*2+8, 8*6+4;
    8*3+7, 8*7+3;
    
    8*4+5, 8*8+1;
    8*5+6, 8*9+2;
    8*6+8, 8*10+4;
    8*7+7, 8*11+3;
    
    8*0+6, 8*5+1;
    8*2+7, 8*7+4;

    8*0+7, 8*3+5;
    8*1+8, 8*6+2;
   
    8*4+3, 8*2+6;
   
    8*4+6, 8*9+1;
    8*5+8, 8*8+3;
    ];

for t=1:size(linkpairs,1)
    [mat,rown]=constrain(mat,rown,linkpairs(t,1),linkpairs(t,2));
end

disp(['rank = ',num2str(rank(mat))])
disp(['DOF = ',num2str(ncube*24-rank(mat))])

% generate plot
v = zeros(8*L*M*N,3);
f = [];
for i = 0:M-1 
    for j = 0:L-1
        for k = 0:N-1
            n = k*(M*L)+L*i + j + 1;
            v(8*n-7,:) = [2*j,2*i,2*k];
            v(8*n-6,:) = [2*j+1.3,2*i,2*k];
            v(8*n-5,:) = [2*j+1.3,2*i+1.3,2*k];
            v(8*n-4,:) = [2*j,2*i+1.3,2*k];
            v(8*n-3,:) = [2*j,2*i,2*k+1.3];
            v(8*n-2,:) = [2*j+1.3,2*i,2*k+1.3];
            v(8*n-1,:) = [2*j+1.3,2*i+1.3,2*k+1.3];
            v(8*n,:) = [2*j,2*i+1.3,2*k+1.3];
            f = [f; 
                8*n-7, 8*n-4, 8*n-5, 8*n-6;
                8*n-3, 8*n-2, 8*n-1, 8*n-0;
                8*n-7, 8*n-3, 8*n-0, 8*n-4;
                8*n-6, 8*n-5, 8*n-1, 8*n-2;
                8*n-7, 8*n-6, 8*n-2, 8*n-3;
                8*n-5, 8*n-4, 8*n-0, 8*n-1];
        end
    end
end

figure; hold on
%  plot the links
for i = 1:length(linkpairs)
    plot3(v(linkpairs(i,:),1), v(linkpairs(i,:),2), v(linkpairs(i,:),3),'Color',[255 51 51]/255,'LineWidth',3);
end
% plot the cubes
patch('Faces',f,'Vertices',v,'FaceColor',[255 229 204]/255,'LineWidth',3);
axis equal tight off
view([10 15])

%%
function [mat, rown] = constrain(mat,rown,i,j)
    mat(rown,i*3-2)=1;
    mat(rown,j*3-2)=-1;
    mat(rown+1,i*3-1)=1;
    mat(rown+1,j*3-1)=-1;
    mat(rown+2,i*3)=1;
    mat(rown+2,j*3)=-1;
    rown = rown+3;
end
