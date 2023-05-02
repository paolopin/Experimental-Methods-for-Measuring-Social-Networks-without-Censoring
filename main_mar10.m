% Code programmed by Paolo Pin on March 10th, 2023

% Here is the explanation of what the code is doing, given by ChatGpt on May 2nd 2023:

% This is a script written in MATLAB that performs various network analysis tasks on different types of networks. Here is a brief overview of what the script does:

% The script reads in two tables A and B using the table2array function and stores them as arrays A and B.
The unique classes (networks) in A and B are identified and stored in arrays classesA and classesB using the unique function.
% For each network in classesA and classesB, the script computes various network metrics such as the number of nodes, average degree, assortativity, diameter, and clustering coefficient. It also runs a simulation of infection spreading using a SI model and records the percentage of infected nodes at each time step.
For the truncated networks in classesB, the script randomly adds edges to each node until the degree of the node is at least 6, and then computes the same network metrics as in step 3.
% For the complete networks in classesA, the script randomly removes edges from each node until the degree of the node is at most 5, and then computes the same network metrics as in step 3.
% The script stores the computed metrics for each type of network (complete, truncated, random addition, and random truncation) in a results array.
% Overall, the script performs network analysis tasks on different types of networks and compares their properties.

% A=table2array(A);
% B=table2array(B);
% 
% classesA=table2array(unique(A(:,2)));
% classesB=table2array(unique(B(:,2)));

% classesA are the complete networks, classesB are the truncated ones
rng(1);
steps=4;
p=.5;
TT=100;

n=[];
av_d=[];
ass=[];
av_du=[];
diam=[];
perc_size=[];
perc_infT=[];
cc=[];


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

    % simulated infections
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
perc_infT=[];
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

    % simulated infections
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
perc_infT=[];
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

%     % simulated infections
%     for ii=1:TT
%         seed=ceil(rand*n(i));
%         infected=find_SImodel(adj2,seed,steps,p);
%         perc_infT((i-1)*TT+ii)=infected/n(i);
%     end
%     %perc_infT(i)=mean(perc_inf);

    % random truncation
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

    % simulated infections
    for ii=1:TT
        seed=ceil(rand*n(i));
        infected=find_SImodel(adjT,seed,steps,p);
        perc_infT((i-1)*TT+ii)=infected/n(i);
    end
    %perc_infT(i)=mean(perc_inf);

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
perc_infT=[];
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
% 
%     % simulated infections
%     for ii=1:TT
%         seed=ceil(rand*n(i));
%         infected=find_SImodel(adj2,seed,steps,p);
%         perc_infT((i-1)*TT+ii)=infected/n(i);
%     end
%     %perc_infT(i)=mean(perc_inf);

    % random addition
    adjM=adj2;
    for jj=1:n(i)
        list=adj2(jj,:);
%         if sum(list)>5
%             adjM(jj,:)=list;
%         else
            links=find(list==0);
            for kk=1:14
                if size(links)>0
                    pick=links(floor(rand*size(links,2))+1);
                end
                adjM(jj,pick)=1;
                links=links(links~=pick);
            end
%        end


    end

    % simulated infections
    for ii=1:TT
        seed=ceil(rand*n(i));
        infected=find_SImodel(adjM,seed,steps,p);
        perc_infT((i-1)*TT+ii)=infected/n(i);
    end
    %perc_infT(i)=mean(perc_inf);

    av_d(i)=mean(sum(adjM));

    ass(i)=assortativity(adjM,1);
    
    adj3=adjM+adjM';
    adj3=(adj3>0);

    av_du(i)=mean(sum(adj3));
    
    C=giant_component(adj3);
    sizeC=size(C,1);

    diam(i)=diameter(C);
    perc_size(i)=sizeC/size(nodes,1);
    
    cc(i)=clust_coeff(adjM);

    if floor(i/100)==i/100
        [i]
    end
end
ClM=[ mean(n)  std(n) , mean(av_d)  std(av_d) , mean(perc_infT)  std(perc_infT) , mean(av_du)  std(av_du) , mean(ass)  std(ass) , mean(diam) std(diam) , mean(perc_size) std(perc_size) , mean(cc) std(cc) ];

results=[ ClA ; ClB ; ClN ; ClM];