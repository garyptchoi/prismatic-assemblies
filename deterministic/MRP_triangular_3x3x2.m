% Construct a minimum rigidifying link pattern (MRP) for a 3x3x2 triangular prismatic assembly
%
% Reference:
% G. P. T. Choi, S. Chen, L. Mahadevan, 
% ``Control of connectivity and rigidity in prismatic assemblies.''
% Preprint, arXiv:2005.05845, 2020.

L = 3; %length
M = 3; %width
N = 2; %height
nprism = L*M*N; %Number of prisms
nlink = 2*L*M*N-2 ;% theoretical lower bound for number of links
mat=zeros(nprism*12+nlink*3,nprism*18); % 18 = 6 vertices x 3 coordinates

% column: x1 y1 z1 x2 y2 z2 ...

% Order of the 6 vertices 
%    6 
%   / \ 
% 4 - - 5
% |  3  |
% | / \ |
% 1 - - 2
% or
% 6 - - 5
% | \ / |
% |  4  |
% 3  |  2
%  \ | /
%    1


% length constraints:
% 9 edge length constraints, and 3 diagonal length constraints
for i=0:nprism-1
    if mod(i,9) == 0 || mod(i,9) == 2  || mod(i,9) == 4  || mod(i,9) == 6  || mod(i,9) == 8 % for 3x3 base only, change if needed
        % edges
        %x1 x2
        mat(i*12+1,i*18+1)=-1; % (x1-x2)
        mat(i*12+1,i*18+4)=1; % (x2-x1)
        %x2 y2 x3 y3
        mat(i*12+2,i*18+4)=1/2; %(x2-x3)
        mat(i*12+2,i*18+5)=-sqrt(3)/2; %(y2-y3)
        mat(i*12+2,i*18+7)=-1/2; % (x3-x2)
        mat(i*12+2,i*18+8)=sqrt(3)/2; %(y3-y2)
        %x1 y1 x3 y3
        mat(i*12+3,i*18+1)=-1/2; 
        mat(i*12+3,i*18+2)=-sqrt(3)/2;
        mat(i*12+3,i*18+7)=1/2;
        mat(i*12+3,i*18+8)=sqrt(3)/2;

        %x4 x5
        mat(i*12+4,i*18+10)=-1;
        mat(i*12+4,i*18+13)=1;
        %x5 y5 x6 y6
        mat(i*12+5,i*18+13)=1/2;
        mat(i*12+5,i*18+14)=-sqrt(3)/2;
        mat(i*12+5,i*18+16)=-1/2;
        mat(i*12+5,i*18+17)=sqrt(3)/2;
        %x4 y4 x6 y6
        mat(i*12+6,i*18+10)=-1/2; 
        mat(i*12+6,i*18+11)=-sqrt(3)/2;
        mat(i*12+6,i*18+16)=1/2; 
        mat(i*12+6,i*18+17)=sqrt(3)/2;

        %z1 z4
        mat(i*12+7,i*18+3)=-1;
        mat(i*12+7,i*18+12)=1;
        %z2 z5
        mat(i*12+8,i*18+6)=-1;
        mat(i*12+8,i*18+15)=1;
        %z3 z6
        mat(i*12+9,i*18+9)=-1;
        mat(i*12+9,i*18+18)=1;

        % face diagonals
        %x1 z1 x5 z5
        mat(i*12+10,i*18+1)=-1; %x1-x5
        mat(i*12+10,i*18+3)=-1; %z1-z5
        mat(i*12+10,i*18+13)=1; %x5-x1
        mat(i*12+10,i*18+15)=1; %z5-z1

        %x2 y2 z2 x6 y6 z6
        mat(i*12+11,i*18+4)=1/2;
        mat(i*12+11,i*18+5)=-sqrt(3)/2;
        mat(i*12+11,i*18+6)=-1;
        mat(i*12+11,i*18+16)=-1/2;
        mat(i*12+11,i*18+17)=sqrt(3)/2;
        mat(i*12+11,i*18+18)=1;

        %x3 y3 z3 x4 y4 z4
        mat(i*12+12,i*18+7)=1/2;
        mat(i*12+12,i*18+8)=sqrt(3)/2;
        mat(i*12+12,i*18+9)=-1;
        mat(i*12+12,i*18+10)=-1/2;
        mat(i*12+12,i*18+11)=-sqrt(3)/2;
        mat(i*12+12,i*18+12)=1;
    else
        % edges
        %x2 x3
        mat(i*12+1,i*18+4)=1;
        mat(i*12+1,i*18+7)=-1;
        %x2 y2 x1 y1
        mat(i*12+2,i*18+4)=1/2;
        mat(i*12+2,i*18+5)=sqrt(3)/2;
        mat(i*12+2,i*18+1)=-1/2;
        mat(i*12+2,i*18+2)=-sqrt(3)/2;
        %x1 y1 x3 y3
        mat(i*12+3,i*18+1)=1/2;
        mat(i*12+3,i*18+2)=-sqrt(3)/2;
        mat(i*12+3,i*18+7)=-1/2;
        mat(i*12+3,i*18+8)=sqrt(3)/2;

        %x5 x6
        mat(i*12+4,i*18+13)=1;
        mat(i*12+4,i*18+16)=-1;
        %x5 y5 x4 y4
        mat(i*12+5,i*18+13)=1/2;
        mat(i*12+5,i*18+14)=sqrt(3)/2;
        mat(i*12+5,i*18+10)=-1/2;
        mat(i*12+5,i*18+11)=-sqrt(3)/2;
        %x4 y4 x6 y6
        mat(i*12+6,i*18+10)=1/2;
        mat(i*12+6,i*18+11)=-sqrt(3)/2;
        mat(i*12+6,i*18+16)=-1/2;
        mat(i*12+6,i*18+17)=sqrt(3)/2;

        %z1 z4
        mat(i*12+7,i*18+3)=-1;
        mat(i*12+7,i*18+12)=1;
        %z2 z5
        mat(i*12+8,i*18+6)=-1;
        mat(i*12+8,i*18+15)=1;
        %z3 z6
        mat(i*12+9,i*18+9)=-1;
        mat(i*12+9,i*18+18)=1;

        % face diagonals
        %x1 y1 z1 x5 y5 z5
        mat(i*12+10,i*18+1)=-1/2;
        mat(i*12+10,i*18+2)=-sqrt(3)/2;
        mat(i*12+10,i*18+3)=-1;
        mat(i*12+10,i*18+13)=1/2;
        mat(i*12+10,i*18+14)=sqrt(3)/2;
        mat(i*12+10,i*18+15)=1;

        %x2 z2 x6 z6
        mat(i*12+11,i*18+4)=1;
        mat(i*12+11,i*18+6)=-1;
        mat(i*12+11,i*18+16)=-1;
        mat(i*12+11,i*18+18)=1;

        %x3 y3 z3 x4 y4 z4
        mat(i*12+12,i*18+7)=-1/2;
        mat(i*12+12,i*18+8)=sqrt(3)/2;
        mat(i*12+12,i*18+9)=-1;
        mat(i*12+12,i*18+10)=1/2;
        mat(i*12+12,i*18+11)=-sqrt(3)/2;
        mat(i*12+12,i*18+12)=1;
    end
end

% Order of the prisms
%       [15] - [16] - [17] 
%      /       /      /|
%    [12] - [13] - [14]|
%   /       /      /   |
% [9] - [10] - [11]    |
%  |                   |
%  |    [6] -  [7] -  [8]
%  |   /      /       /
%  | [3] -  [4] -  [5]
%  |/      /      /
% [0] -  [1] -  [2]

rown=nprism*12+1;
linkpairs=[
    6*0+2, 6*1+1;
    6*1+1, 6*2+1;
    6*2+3, 6*5+1;
    6*4+3, 6*5+3;
    6*3+2, 6*4+3;
    6*0+3, 6*3+1;
    6*9+5, 6*10+4;
    6*10+4, 6*11+4;
    6*11+6, 6*14+4;
    6*9+6, 6*12+4;
    6*0+4, 6*9+1;
    6*2+5, 6*11+2;
    6*3+6, 6*12+3;
    6*5+5, 6*14+2;
    6*1+4, 6*9+2;
    6*14+1, 6*10+2;
    6*2+6, 6*4+5;
    6*4+4, 6*12+1;
    6*1+6, 6*13+1;
    6*10+6, 6*13+4;
    6*3+3, 6*6+1;
    6*6+3, 6*7+3;
    6*7+2, 6*8+3;
    6*5+2, 6*8+2;
    6*12+6, 6*15+4;
    6*15+6, 6*16+6;
    6*16+5, 6*17+6;
    6*14+5, 6*17+5;
    6*6+6, 6*15+3;
    6*8+6, 6*17+3;
    6*12+2, 6*7+4;
    6*13+3, 6*17+1;
    6*16+4, 6*14+6;
    6*13+6, 6*15+5;
    ];


for t=1:size(linkpairs,1)
    [mat,rown]=constrain(mat,rown,linkpairs(t,1),linkpairs(t,2));
end

disp(['rank = ',num2str(rank(mat))])
disp(['DOF = ',num2str(nprism*18-rank(mat))])

% generate plot
v = zeros(6*L*M*N,3);
f = [];
edgelength = 2;
vs = 0.6;
height = 1.5;

for i = 0:M-1 
    for j = 0:L-1
        for k = 0:N-1
            n = k*(M*L)+L*i + j + 1;
            
            if mod(i,2) ==  mod(j,2)

                v(6*n-5,:) = [2*j-edgelength/2, 2*i+i*vs,2*k];
                v(6*n-4,:) = [2*j+edgelength/2, 2*i+i*vs,2*k];
                v(6*n-3,:)   = [2*j , 2*i+edgelength*sqrt(3)/2+i*vs,2*k];
                v(6*n-2,:) = [2*j-edgelength/2, 2*i+i*vs,2*k+height];
                v(6*n-1,:) = [2*j+edgelength/2, 2*i+i*vs,2*k+height];
                v(6*n,:)   = [2*j , 2*i+edgelength*sqrt(3)/2+i*vs,2*k+height];
            else
                v(6*n-5,:) = [2*j, 2*i+i*vs,2*k];
                v(6*n-4,:) = [2*j+edgelength/2, 2*i+edgelength*sqrt(3)/2+i*vs,2*k];
                v(6*n-3,:) = [2*j-edgelength/2, 2*i+edgelength*sqrt(3)/2+i*vs,2*k];
                v(6*n-2,:) = [2*j, 2*i+i*vs,2*k+height];
                v(6*n-1,:) = [2*j+edgelength/2, 2*i+edgelength*sqrt(3)/2+i*vs,2*k+height];
                v(6*n,:) = [2*j-edgelength/2, 2*i+edgelength*sqrt(3)/2+i*vs,2*k+height];
            end
            f = [f; 
                6*n-5 6*n-4 6*n-3 6*n-3;
                6*n-2 6*n-1 6*n 6*n;
                6*n-5 6*n-4 6*n-1 6*n-2;
                6*n-4 6*n-3 6*n 6*n-1;
                6*n-3 6*n-5 6*n-2 6*n;
                ];
        end
    end
end

fig = figure; hold on
%  plot the links
for i = 1:length(linkpairs)
    plot3(v(linkpairs(i,:),1), v(linkpairs(i,:),2), v(linkpairs(i,:),3),'Color',[255 51 51]/255,'LineWidth',3);
end
% plot the prisms
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
