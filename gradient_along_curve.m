
%will fit a straight line on every possible 5 consecutives points and find the gradient of the fitted line.
%input tau in first column of 'var'. Do not include 0.
%input MSD/tau in second column of 'var'. Do not include 0.
%% Written by D.S.Jokhun on 05.08.2016

smoothening_pts=5   %num of points to use for finding tangent. Preferable to use odd number but can use even as well (default is 5, can use 9 if required)
linear_range=60:99  ;  %the range of x (in terms of row numbers) for which gradient of the straight line is required. Gradint for tau=15s to tau=25s for telomeres

log_var=[];
log_var(:,1)=log(var(:,1));
log_var(:,2)=log(var(:,2));

tangent=[];
gradient_series=[];
for fragment=1:size(log_var,1)-(smoothening_pts-1)   %we will take tangent from n points
    tangent(1:smoothening_pts,1)=log_var(fragment:fragment+(smoothening_pts-1),1);
    tangent(1:smoothening_pts,2)= polyval( polyfit(log_var(fragment:fragment+(smoothening_pts-1),1),log_var(fragment:fragment+(smoothening_pts-1),2),1) , log_var(fragment:fragment+(smoothening_pts-1),1) );
    dy_dx=diff(tangent(:,2))./diff(tangent(:,1));
    gradient_series(size(gradient_series,1)+1,1)=exp(log_var(fragment+floor(smoothening_pts/2),1));
    gradient_series(size(gradient_series,1),2)=mean(dy_dx);
end

gradient_series(:,2)=smooth(gradient_series(:,2),8,'moving');


%%% finding gradient of the line fitting the points within the range specified by linear_range
x=[];
y=[];
x(:,1)=log_var(linear_range,1);
y(:,1)=log_var(linear_range,2);
[coeff,~,~,~,~] = pca([x,y]);
gradient_of_PCA=(coeff(2,1)/coeff(1,1))
X=mean([var(min(linear_range)),var(max(linear_range))])
alpha_at_specific_point=gradient_of_PCA+1



% figure
% loglog(var(:,1),var(:,2))
% 
% figure
% semilogx(gradient_series(:,1),gradient_series(:,2));

alpha_series=[];
alpha_series(:,1)=gradient_series(:,1);
alpha_series(:,2)=gradient_series(:,2)+1;
% figure
% semilogx(alpha_series(:,1),alpha_series(:,2));