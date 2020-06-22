% Construct a minimum connecting link pattern (MCP) for a 2x2x3 rectangular prismatic assembly
%
% Reference:
% G. P. T. Choi, S. Chen, L. Mahadevan, 
% ``Control of connectivity and rigidity in prismatic assemblies.''
% Preprint, arXiv:2005.05845, 2020.

L = 2; %length
M = 2; %width
N = 3; %height
ncube = L*M*N; %Number of cubes

% Order of the 8 vertices
%   8 - - 7
%  /     /|
% 5 - - 6 |
% | 4 - | 3
% |/    |/
% 1 - - 2

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
    8*0+3, 8*2+2;
    8*1+3, 8*3+2;
   
    8*0+5, 8*4+1;
    
    8*4+2, 8*5+1;
    8*4+3, 8*6+2;
    8*5+3, 8*7+2;
    
    8*4+5, 8*8+1;
    
    8*8+2, 8*9+1;
    8*8+3, 8*10+2;
    8*9+3, 8*11+2;
    ];

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

fig = figure; hold on
%  plot the links
for i = 1:length(linkpairs)
    plot3(v(linkpairs(i,:),1), v(linkpairs(i,:),2), v(linkpairs(i,:),3),'Color',[255 51 51]/255,'LineWidth',3);
end
% plot the cubes
patch('Faces',f,'Vertices',v,'FaceColor',[255 229 204]/255,'LineWidth',3);
axis equal tight off
view([10 15])
