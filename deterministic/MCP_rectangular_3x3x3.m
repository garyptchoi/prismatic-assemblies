% Construct a minimum connecting link pattern (MCP) for a 3x3x3 rectangular prismatic assembly
%
% Reference:
% G. P. T. Choi, S. Chen, L. Mahadevan, 
% ``Control of connectivity and rigidity in prismatic assemblies.''
% Proceedings of the Royal Society A, 476(2244), 20200485, 2020.

L = 3; %length
M = 3; %width
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
%       [24] - [25] - [26] 
%      /       /      /|
%    [21] - [22] - [23]|
%   /       /      /   |
% [18] - [19] - [20]   |
%  |                   |
%  |    [15] - [16] - [17]
%  |   /       /      /|
%  | [12] - [13] - [14]|
%  |/       /      /   |
% [9] - [10] - [11]    |
%  |                   |
%  |    [6] -  [7] -  [8]
%  |   /      /       /
%  | [3] -  [4] -  [5]
%  |/      /      /
% [0] -  [1] -  [2]


rown=ncube*18+1;
linkpairs=[
    8*0+2, 8*1+1;
    8*1+2, 8*2+1;
    8*0+3, 8*3+2;
    8*1+3, 8*4+2;
    8*2+3, 8*5+2;
    8*3+3, 8*6+2;
    8*4+3, 8*7+2;
    8*5+3, 8*8+2;
    
    8*0+5, 8*9+1;
    
    8*9+2, 8*10+1;
    8*10+2, 8*11+1;
    8*9+3, 8*12+2;
    8*10+3, 8*13+2;
    8*11+3, 8*14+2;
    8*12+3, 8*15+2;
    8*13+3, 8*16+2;
    8*14+3, 8*17+2;
    
    8*9+5, 8*18+1;
    
    8*18+2, 8*19+1;
    8*19+2, 8*20+1;
    8*18+3, 8*21+2;
    8*19+3, 8*22+2;
    8*20+3, 8*23+2;
    8*21+3, 8*24+2;
    8*22+3, 8*25+2;
    8*23+3, 8*26+2;
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
