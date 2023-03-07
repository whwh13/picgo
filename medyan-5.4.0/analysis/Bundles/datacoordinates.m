% This code reads an input snapshot file. First it displays the whole
% system in the last frame, then it displays a biograph of the 
% filaments that are close (adjustable parameter). The final graph
% is of the bundles that are categorized by closeness and 
% alignment (adjustable parameter). 

% Papoian Lab, Univeristy of Maryland 7/9/2015 


% Script to read data points of each bead

fclose all;
clear all;close all;clf


% open the text file
f1=fopen('snapshottrial.traj','r');


% read each line of the text file in the cells of A
i=1;
while (~feof(f1))
    
    A{i}= fgetl(f1);
    i=i+1;
    
end


% Separate each frame whenever an empty line is read, or whenever a
% cell of A is empty
index=find(strcmp(A,''));

if zeros(1,6)==str2num(A{1})
    % if the starting point is initialized at zero concentrations
    frame{1}{1}=A{1};
    
    for j=1:length(index)-1
        i=1;
        for k=index(j)+1:index(j+1)-1
            frame{j+1}{i}=A{k};
            i=i+1;
        end
    end
    
else
    % if starting point is not initialized at zero
    for n=1:index(1)-1
        frame{1}{n}=A{n};
    end
    
    for j=1:length(index)-1
        i=1;
        for k=index(j)+1:index(j+1)-1
            frame{j+1}{i}=A{k};
            i=i+1;
        end
    end
    
end


% Assign each filament coordinates to the cells of d.
% For example d{1}{1}(1) is the x-coordinate of the first filament
% in the first frame.
for ii=1:length(frame)
    
    if zeros(1,6)==str2num(A{1})
        d{1}{1}=str2num(frame{1}{1});
    end
    
    m=length(frame{ii});
    a=1;
    for jj=1:m
        if frame{ii}{jj}(1)=='F'
            d{ii}{a}=str2num(frame{ii}{jj+1});
            a=a+1;
        end
    end
    
    
end


% Assign the x-, y-, and z- coordinates of each filament in the last
% frame to matrices xcoord, ycoord, zcoord.
for i=1:length(d{end})
    x=1;
    y=1;
    z=1;
    for k=1:3:length(d{end}{i})
        
        xcoord{i}{x}=d{end}{i}(k);
        x=x+1;
        
    end
    
    for k=2:3:length(d{end}{i})
        
        ycoord{i}{y}=d{end}{i}(k);
        y=y+1;
        
    end
    
    for k=3:3:length(d{end}{i})
        
        zcoord{i}{z}=d{end}{i}(k);
        z=z+1;
        
    end
    
end


% Now that the x-, y-, and z- coordinates are assigned to their
% respective cells, we can compute the best fit line for each
% filament of the last frame. The normal vectors, fitting errors,
% and a point on each line is output in N,Err, and P respectively.
hold on
figure(1)
for i = 1 : length(xcoord)
    clear xc yc zc
    
    for j = 1 : length(xcoord{i})
        
        xc(1,j)=xcoord{i}{j};
        yc(1,j)=ycoord{i}{j};
        zc(1,j)=zcoord{i}{j};
        
    end
    xc=xc';
    yc=yc';
    zc=zc';
    [Err(i),N(:,i),P(:,i)]=fit_3D_data(xc,yc,zc,'line','on','off');
    
    
    
end


% Check the shortest distance between two filaments. If the two
% filaments intersect or if the distance between the two filaments
% is less than a certain threshold, T, then we consider the
% two to be close. If two filaments are close then we place a 1
% in the matrix entry corresponding to the indicies of the filaments.
% If the two filaments are not close, then we place a 0 in that entry

T=5;

for i=1:length(xcoord)
    
    for j=1:i
        
        P1=[xcoord{i}{1},ycoord{i}{1},zcoord{i}{1}];
        P2=[xcoord{i}{end},ycoord{i}{end},zcoord{i}{end}];
        
        P3=[xcoord{j}{1},ycoord{j}{1},zcoord{j}{1}];
        P4=[xcoord{j}{end},ycoord{j}{end},zcoord{j}{end}];
        
        dist = DistBetween2Segment(P1, P2, P3, P4);
        
        if dist <= T
            B1(i,j)=1;
        else
            B1(i,j)=0;
        end
    end
end


% Now that we have a sparse matrix indicating the connectivity
% of each filament based on distance, we can find the filaments
% that we consider are in a bundle based on their closeness.
% First we remove the diagonals of B1 so that the graph will not
% display repitions. Then we display an image of the connections
% based on closeness. Lastly, we find the components of the
% filaments that we consider close. The list of bundles considered
% close are in the cell distbundles
I=eye(size(B1));C=B1-I;
h = view(biograph(C));
[nComponents,sizes,distbundles] = networkComponents(C);



% Now that we have the bundles that are close together, we can
% further characterize bundles by their alignmet. Using the criteria
% for the Nematic Director, found at:
% http://cmt.dur.ac.uk/sjc/thesis_dlc/node65.html, look at the
% alignment of the bundles with two or more filaments.
for i=1:length(distbundles)
    
    if length(distbundles{i})>1
        
        Nm=length(distbundles{i});
        
        for j=1:length(distbundles{i})
            
            for k=1:j
                
                if k==j
                    delta=1;
                else
                    delta=0;
                end
                
                Q(j,k)=(1/Nm)*(3/2)*...
                    abs(dot(N(:,distbundles{i}(j)),...
                    N(:,distbundles{i}(k))))-(1/2)*delta;
                
            end
            
        end
        
        QQ{i}=Q+transpose(Q)-triu(Q);
        clear Q
    end
    
end

% Next we find the indicies of the filaments in each bundle. The
% Parameter Ang can be adjusted depending on the angle we want to
% consider. If we want an alignment angle of 30 degrees or smaller,
% the value of Ang should be ~0.75. For angle of 15 degrees or
% smaller, Ang should be ~ 0.9330. As Ang approaches a value of 1
% the alignment becomes more parallel.

Ang = 0.4;

n=1;
for i=1:length(QQ)
    
    if max(eig(QQ{i}))>= Ang
        Bundles{n}=distbundles{i};
        n=n+1;
    end
    
end


% Now we display all the bundles using the best fit script.
hold on;
figure(2)
for i=1:length(Bundles)
    
    for j=1:length(Bundles{i})
        
        clear xb yb zb
        
        for k = 1 : length(xcoord{Bundles{i}(j)})
            
            xb(k,1)=xcoord{Bundles{i}(j)}{k};
            yb(k,1)=ycoord{Bundles{i}(j)}{k};
            zb(k,1)=zcoord{Bundles{i}(j)}{k};
            
        end
        
        [Errb(i),Nb(:,i),Pb(:,i)]=fit_3D_data(xb,yb,zb,'line','on','off');               
        
    end
end
title('System with Only Bundles')


