function Vaff = affichage(V,perrow,Li,Co,bw)

% Display of NMF solutions, for image datasets. 
% More precisely, each columnof V is a vectorized gray level image, and 
% affichage displays these images of dimension LixCo in a grid with lig
% images per row. 
%
% Vaff = affichage(V,lig,Li,Co,bw)
%
% Input.
%   V           : (m x r) matrix whose colums contains vectorized images
%  perrow       : number of images per row in the display
%   (Li,Co)     : dimensions of images
%  bw           : if bw = 0: high intensities in black (default) 
%               :    bw = 1: high intensities in white 
%
% Output.
%   Vaff is the matrix which allows to diplay columns of matrix V as images 
%   in a grid with lig images per row.
% 
% Note: affichage in French means display 

if nargin == 4
    bw = 0;
end
V = max(V,0); 
[m,r] = size(V); 
% Normalize columns to have maximum entry equal to 1 
for i = 1 : r
    V(:,i) = V(:,i)/max(V(:,i));
end
% Construct the image Vaff to display 
Vaff = []; 
lin = ceil(r/perrow);
col = perrow;
lastindli = 1; 
lastindcol = 1; 
i = 1; 
for lini = 1 : lin
    for coli = 1 : col
        if coli == perrow
            Vaff(lastindli:lastindli+Li-1,lastindcol:lastindcol+Co-1) ... 
                = reshape(V(:,i),Li,Co); 
            Vaff = [Vaff; ones(1,size(Vaff,2))]; 
            lastindcol = 1; 
            lastindli = size(Vaff,1)+1; 
        else
            Vaff(lastindli:lastindli+Li-1,lastindcol:lastindcol+Co) ... 
                = [reshape(V(:,i),Li,Co) ones(Li,1)]; 
            lastindcol = lastindcol+Co+1; 
        end
        i = i+1; 
        if i > r
            break;
        end
    end
end

[m,n] = size(Vaff);
for i = 1 : m/Li-1
        Vaff2 = [Vaff(1:Li*i+i-1,:); ones(20,n); Vaff(Li*i+i:end,:)];
end
% Display image taking about 1/4 of screen at most
figure;  
[hi,wi] = size(Vaff); 
TwiThi = get(0,'ScreenSize'); 
maxs = ceil(0.5*[TwiThi(3), TwiThi(4)]);
disratio = maxs(1)/maxs(2); 
if wi >= hi*disratio
    widthi = maxs(1); 
    heighti = ceil(widthi*hi/wi)+40; % +40 for possible title 
else
    heighti = maxs(2); 
    widthi = ceil(heighti/hi*wi)+40; % +40 for possible title
end
set(gcf, 'Position',  [ceil(TwiThi(3)*0.1), ceil(TwiThi(4)*0.1), widthi, heighti])
if bw == 1
    imagesc(Vaff,[0 1]); 
else
    imagesc(1-Vaff,[0 1]);
end
colormap(gray); 
set(gca,'XTick',[],'YTick',[]); 
warning('on');