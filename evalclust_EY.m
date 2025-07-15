function [sum_cen, r_sum_cen,shuff_comp_sum_cen]=evalclust_EY(mat,k);
%mat is the matrix of DFF calculated cells x time
% k is up to how many k's you wish to run the iteration
% the key plot to do using the out put is diff_cen, this will give us the distribution when adding new 'k's doesnt improve the distance to centroid calculations
% for Ewelina's paper, for all animals (mean +/- sem) we can plot a graph of sum_cen and r_sum_cen. This will give us the comparision of real data and the shuffleed data simulation
% for Ewelina's paper, we can also plot (for all animals, mean +/- sem),
% suff_comp)sum_cen... the peak of this graph will give us the optimal number clusters, both for HB and for TEL . 
% Keep in mind HB is not that time consuming but TEL takes awful amount of time

clust = zeros(size(mat,1),k);
sum_cent = zeros(1,k);
% clust = zeros(size(mat,1),k);
for i=1:k
[clust(:,i), C, sumd, D] = kmeans(mat,i,'Distance', 'correlation', 'MaxIter',100, 'replicates',50);
sum_cen(:,i)=sum(sumd)/mean(sum(D));
end

diff_cen=diff(sum_cen);

%for the elbow method see here: https://towardsdatascience.com/k-means-clustering-algorithm-applications-evaluation-methods-and-drawbacks-aa03e644b48a
% figure, subplot (1,2,1), plot (sum_cen), title('sumd/avergare of all D') , subplot (1,2,2), plot (diff(sum_cen)), title('slope of sumd/avergare of all D');

%% this section was the crappy evaluation method in MATLAB but I left it still 
% eva = evalclusters(mat,clust,'CalinskiHarabasz')
% eva = evalclusters(mat,clust,'silhouette') ;
% super slow
% eva = evalclusters(mat,'kmeans','gap','KList',[1:k]);

%%
% here we are shuffling the time course of each neurons seperately
r_mat= zeros(size(mat));
for i=1:size(mat,1)
r_mat(i,:)=mat(i,randperm(size(mat,2)));
end

% from here of classical clutering with this shuffled data
r_clust = zeros(size(r_mat,1),k);
r_sum_cent = zeros(1,k);

for i=1:k
[r_clust(:,i), r_C, r_sumd, r_D] = kmeans(r_mat,i,'Distance', 'correlation', 'MaxIter',100, 'replicates',50);
r_sum_cen(:,i)=sum(r_sumd)/mean(sum(r_D));
end

r_diff_cen=diff(r_sum_cen);

shuff_comp_sum_cen=r_sum_cen-sum_cen; % this is the subrtaction for the sumd
shuff_comp_diff=diff(r_sum_cen)-diff(sum_cen); % this is the subrtaction for the diff(sumd)

%% this section was the crappy evaluation method in MATLAB but I left it still 
% eva = evalclusters(r_mat,r_clust,'CalinskiHarabasz')
% eva = evalclusters(r_mat,r_clust,'silhouette') ;

%%
figure,
% subplot(1,2,1), plot(eva);

 
subplot (12,2,1), imagesc(mat(find(clust(:,6)==1),:)), colormap jet, colorbar, caxis([-10 30]), set(gca,'XTick',[]),set(gca,'YTick',[]),box off, colormap jet % set(gca,'XTick',[]) % title (num2str(1)) 
subplot (12,2,3), imagesc(mat(find(clust(:,6)==2),:)), colormap jet, colorbar, caxis([-10 30]),set(gca,'XTick',[]),set(gca,'YTick',[]),box off % set(gca,'XTick',[]) % title (num2str(2))
subplot (12,2,5), imagesc(mat(find(clust(:,6)==3),:)), colormap jet, colorbar, caxis([-10 30]), set(gca,'XTick',[]),set(gca,'YTick',[]),box off% set(gca,'XTick',[]) % title (num2str(3))
subplot (12,2,7), imagesc(mat(find(clust(:,6)==4),:)), colormap jet, colorbar, caxis([-10 30]), set(gca,'XTick',[]),set(gca,'YTick',[]),box off% set(gca,'XTick',[]) % title (num2str(4))
subplot (12,2,9), imagesc(mat(find(clust(:,6)==5),:)), colormap jet, colorbar, caxis([-10 30]), set(gca,'XTick',[]),set(gca,'YTick',[]),box off% set(gca,'XTick',[]) % title (num2str(5))
subplot (12,2,11), imagesc(mat(find(clust(:,6)==6),:)), colormap jet, colorbar, caxis([-10 30]),set(gca,'YTick',[]),box off % set(gca,'XTick',[]) % title (num2str(2)) 

subplot (12,2,2), imagesc(r_mat(find(r_clust(:,6)==1),:)), colormap jet, colorbar, caxis([-10 30]), set(gca,'XTick',[]),set(gca,'YTick',[]),box off, colormap jet % set(gca,'XTick',[]) % title (num2str(1)) 
subplot (12,2,4), imagesc(r_mat(find(r_clust(:,6)==2),:)), colormap jet, colorbar, caxis([-10 30]),set(gca,'XTick',[]),set(gca,'YTick',[]),box off % set(gca,'XTick',[]) % title (num2str(2))
subplot (12,2,6), imagesc(r_mat(find(r_clust(:,6)==3),:)), colormap jet, colorbar, caxis([-10 30]), set(gca,'XTick',[]),set(gca,'YTick',[]),box off% set(gca,'XTick',[]) % title (num2str(3))
subplot (12,2,8), imagesc(r_mat(find(r_clust(:,6)==4),:)), colormap jet, colorbar, caxis([-10 30]), set(gca,'XTick',[]),set(gca,'YTick',[]),box off% set(gca,'XTick',[]) % title (num2str(4))
subplot (12,2,10), imagesc(r_mat(find(r_clust(:,6)==5),:)), colormap jet, colorbar, caxis([-10 30]), set(gca,'XTick',[]),set(gca,'YTick',[]),box off% set(gca,'XTick',[]) % title (num2str(5))
subplot (12,2,12), imagesc(r_mat(find(r_clust(:,6)==6),:)), colormap jet, colorbar, caxis([-10 30]),set(gca,'YTick',[]),box off % set(gca,'XTick',[]) % title (num2str(2)) 


subplot (2,4,5), plot (sum_cen), hold on, plot (r_sum_cen,'r'), hold off,  title('sumd/avergare of all D'), hold on , 
subplot (2,4,7), plot (diff(sum_cen)),hold on, plot (diff(r_sum_cen)), hold off, title('slope of sumd/avergare of all D');
subplot (2,4,6), plot (r_sum_cen-sum_cen), title('difference to random sumd/avergare of all D'),
subplot (2,4,8), plot (diff(r_sum_cen)-diff(sum_cen)), title('slope of sumd/avergare of all D');