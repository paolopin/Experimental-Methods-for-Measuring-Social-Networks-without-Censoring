% Code programmed by Paolo Pin on August 26th, 2022

% Here is the explanation of what the code is doing, given by ChatGpt on May 2nd 2023:

% This MATLAB code is performing network analysis on different types of networks, including complete networks, truncated networks, random truncation networks, and random addition networks. The code starts by setting the random number generator seed to 1 and defining some parameters such as the number of steps, probability, and number of iterations.

% The code then proceeds to analyze complete networks. It first extracts the edges of each class of nodes and creates a sorted list of unique nodes. It initializes an adjacency matrix and populates it with the edges from the extracted class of nodes. It calculates the number of nodes, average degree, infection rate, average degree after truncation, assortativity, diameter, and the size of the giant component for each class of nodes. It repeats this process for truncated networks and records the corresponding network properties.

% The code then analyzes random truncation and random addition networks. For each network, it randomly selects a class of nodes and extracts the edges. It then constructs an adjacency matrix and calculates the number of nodes, infection rate, average degree after truncation, assortativity, diameter, and the size of the giant component. The code repeats this process for 1000 iterations and records the corresponding network properties.

% Finally, the code outputs the mean and standard deviation of the network properties for each type of network. 


% A=table2array(A);
% B=table2array(B);
% 
% classesA=table2array(unique(A(:,2)));
% classesB=table2array(unique(B(:,2)));

rng(1);
steps=4;
p=.5;
TT=100;

%% COMPLETE NETWORKS
ii=size(classesA);
for i=1:ii
    a=A(A(:,2)==classesA(i),:);
    nodes=sort(unique(a(:,3:4)));
    
    adj=zeros(numel(nodes));   % initialize adjacency matrix
    % across all edges
    for j=1:size(a,1); adj(find(nodes==a(j,3)),find(nodes==a(j,4)))=a(j,5); end
    
    n(i)=size(nodes,1);
    
    adj2=(adj>0);
    av_d(i)=mean(sum(adj2));

    for ii=1:TT
        seed=ceil(rand*n(i));
        infected=find_SImodel(adj2,seed,steps,p);
        perc_infT((i-1)*TT+ii)=infected/n(i);
    end
    %perc_infT(i)=mean(perc_inf);

    adj3=adj2+adj2';
    adj3=(adj3>0);

    av_du(i)=mean(sum(adj3));

    ass(i)=assortativity(adj3,0);
    
    C=giant_component(adj3);
    sizeC=size(C,1);

    diam(i)=diameter(C);
    perc_size(i)=sizeC/size(nodes,1);
    
    cc(i)=clust_coeff(adj3);
end
ClA=[ mean(n)  std(n) , mean(av_d)  std(av_d) , mean(perc_infT)  std(perc_infT) , mean(av_du)  std(av_du) , mean(ass)  std(ass) , mean(diam) std(diam) , mean(perc_size) std(perc_size) , mean(cc) std(cc) ];

n=[];
av_d=[];
ass=[];
av_du=[];
diam=[];
perc_size=[];
cc=[];


%% TRUNCATED NETWORKS

ii=size(classesB);
for i=1:ii
    a=B(B(:,2)==classesB(i),:);
    nodes=sort(unique(a(:,3:4)));
    
    adj=zeros(numel(nodes));   % initialize adjacency matrix
    % across all edges
    for j=1:size(a,1); adj(find(nodes==a(j,3)),find(nodes==a(j,4)))=a(j,5); end
    
    n(i)=size(nodes,1);
    
    adj2=(adj>0);
    av_d(i)=mean(sum(adj2));

    for ii=1:TT
        seed=ceil(rand*n(i));
        infected=find_SImodel(adj2,seed,steps,p);
        perc_infT((i-1)*TT+ii)=infected/n(i);
    end
    %perc_infT(i)=mean(perc_inf);
    
    adj3=adj2+adj2';
    adj3=(adj3>0);

    av_du(i)=mean(sum(adj3));

    ass(i)=assortativity(adj3,0);
    
    C=giant_component(adj3);
    sizeC=size(C,1);

    diam(i)=diameter(C);
    perc_size(i)=sizeC/size(nodes,1);
    
    cc(i)=clust_coeff(adj3);
end
ClB=[ mean(n)  std(n) , mean(av_d)  std(av_d) , mean(perc_infT)  std(perc_infT) , mean(av_du)  std(av_du) , mean(ass)  std(ass) , mean(diam) std(diam) , mean(perc_size) std(perc_size) , mean(cc) std(cc) ];

n=[];
av_d=[];
ass=[];
av_du=[];
diam=[];
perc_size=[];
cc=[];


%% Random truncation NETWORKS
for i=1:1000
    ii=size(classesA);
    c=(floor(rand*ii(1)))+1;
    a=A(A(:,2)==classesA(c),:);
    nodes=sort(unique(a(:,3:4)));
    
    adj=zeros(numel(nodes));   % initialize adjacency matrix
    % across all edges
    for j=1:size(a,1); adj(find(nodes==a(j,3)),find(nodes==a(j,4)))=a(j,5); end

    n(i)=size(nodes,1);
    
    adj2=(adj>0);

    for ii=1:TT
        seed=ceil(rand*n(i));
        infected=find_SImodel(adj2,seed,steps,p);
        perc_infT((i-1)*TT+ii)=infected/n(i);
    end
    %perc_infT(i)=mean(perc_inf);

    adjT=zeros(n(i),n(i));
    for jj=1:n(i)
        list=adj2(jj,:);
        if sum(list)<6
            adjT(jj,:)=list;
        else
            links=find(list>0);
            for kk=1:5
                pick=links(floor(rand*size(links,2))+1);
                adjT(jj,pick)=1;
                links=links(links~=pick);
            end
        end


    end

    av_d(i)=mean(sum(adjT));

    ass(i)=assortativity(adjT,1);
    
    adj3=adjT+adjT';
    adj3=(adj3>0);

    av_du(i)=mean(sum(adj3));
    
    C=giant_component(adj3);
    sizeC=size(C,1);

    diam(i)=diameter(C);
    perc_size(i)=sizeC/size(nodes,1);
    
    cc(i)=clust_coeff(adjT);

    if floor(i/100)==i/100
        [i]
    end
end
ClN=[ mean(n)  std(n) , mean(av_d)  std(av_d) , mean(perc_infT)  std(perc_infT) , mean(av_du)  std(av_du) , mean(ass)  std(ass) , mean(diam) std(diam) , mean(perc_size) std(perc_size) , mean(cc) std(cc) ];

n=[];
av_d=[];
ass=[];
av_du=[];
diam=[];
perc_size=[];
cc=[];

%% Random ADDITION NETWORKS

for i=1:1000
    ii=size(classesB);
    c=(floor(rand*ii(1)))+1;
    a=B(B(:,2)==classesB(c),:);
    nodes=sort(unique(a(:,3:4)));
    
    adj=zeros(numel(nodes));   % initialize adjacency matrix
    % across all edges
    for j=1:size(a,1); adj(find(nodes==a(j,3)),find(nodes==a(j,4)))=a(j,5); end

    n(i)=size(nodes,1);
    
    adj2=(adj>0);

    for ii=1:TT
        seed=ceil(rand*n(i));
        infected=find_SImodel(adj2,seed,steps,p);
        perc_infT((i-1)*TT+ii)=infected/n(i);
    end
    %perc_infT(i)=mean(perc_inf);

    adjT=zeros(n(i),n(i));
    for jj=1:n(i)
        list=adj2(jj,:);
        if sum(list)<6
            adjT(jj,:)=list;
        else
            links=find(list>0);
            for kk=1:5
                pick=links(floor(rand*size(links,2))+1);
                adjT(jj,pick)=1;
                links=links(links~=pick);
            end
        end


    end

    av_d(i)=mean(sum(adjT));

    ass(i)=assortativity(adjT,1);
    
    adj3=adjT+adjT';
    adj3=(adj3>0);

    av_du(i)=mean(sum(adj3));
    
    C=giant_component(adj3);
    sizeC=size(C,1);

    diam(i)=diameter(C);
    perc_size(i)=sizeC/size(nodes,1);
    
    cc(i)=clust_coeff(adjT);

    if floor(i/100)==i/100
        [i]
    end
end
ClN=[ mean(n)  std(n) , mean(av_d)  std(av_d) , mean(perc_infT)  std(perc_infT) , mean(av_du)  std(av_du) , mean(ass)  std(ass) , mean(diam) std(diam) , mean(perc_size) std(perc_size) , mean(cc) std(cc) ];

