clc;
%% Part 3
% Loading the stimulus file
load('C:\Users\Asus\Desktop\Stimulus_Files\msq1D.mat');
% Reading the Neuron
neuron_struct = Func_ReadData ('000412.a01') ;
% Initializing some matrixes
Spike_Triggered = [] ;
MAX = [];
corrcont = [] ;
corr = [] ;
% Calculating " Spike triggered average "
for i = 1:length (neuron_struct) 
    events = neuron_struct(i).events ;
    Spike_Triggered_add = Func_StimuliExtraction( events , msq1D , 59.721395 ) ;
    Spike_Triggered = cat (3 ,Spike_Triggered , Spike_Triggered_add ) ;
    MAX = [MAX max(neuron_struct(i).events)] ;
end
STA = mean (Spike_Triggered , 3) ;
% Creating Control Matrix
for j = 1:length (Spike_Triggered)
    corr_add = sum(sum((STA).*(Spike_Triggered ( : , : , j ) ))) ;
    corr = [corr corr_add] ;
end
control_events = ceil(max(MAX)*rand(1,length(Spike_Triggered))) ;
control_matrix3D = Func_StimuliExtraction(control_events , msq1D , 59.721395 ) ;
for k = 1 : length (control_matrix3D)
   corrcont_add =  sum(sum(STA .* control_matrix3D(: , : , k))) ;
   corrcont = [corrcont corrcont_add] ;
end
% calculating p_values
p_value = zeros(16,16) ;
for l = 1:16
    for m = 1: 16
     [a p_value(l,m)] = ttest(Spike_Triggered(l,m,:)) ;
    end
end
% running hypothesis test for the two imaged vectors 
[p h] = ttest ( corr , corrcont) ;
h
%% fitting gaussians on datas
mean (corr);
mu1 = mean (corr);
mu2 = mean (corrcont);
sigma1 = var (corr);
sigma2 = var (corrcont);
dist1 =@(x) exp(-(x-mu1).^2 / (2*sigma1^2)) / sqrt(2*sigma1^2*pi);
dist2 =@(x) exp(-(x-mu2).^2 / (2*sigma2^2)) / sqrt(2*sigma2^2*pi);
threshold = fzero(@(x) dist1(x) - dist2(x), rand * (mu1 - mu2) + (mu1 + mu2))
STA_Success_Percentage = length(find((dist1(corr)-dist2(corr))>0))/length (corr)
%% plotting all the resualts 
figure
subplot (1,2,1)
imshow(STA, [-1 1],'InitialMagnification','fit') 
title ('Spike Triggered Average')
xlabel('Spatial')
ylabel('Temporal')
subplot (1,2,2)
imshow(1-p_value,'InitialMagnification','fit')
title ('1-p value')
xlabel('Spatial')
ylabel('Temporal')
figure
histogram(corr,'Normalization','probability','BinWidth',0.1)
hold on
histogram(corrcont,'Normalization','probability','BinWidth',0.1)
hold off
%% part 4
%% Creating the Correlation Matrix
correlation_matrix = zeros (256 , 256) ; % Initializing the correlation matrix
t = size (Spike_Triggered,3) ; % Number of all spikes
for u = 1 : 256 % uth element of the stimulus
    i = (ceil(u / 16)) ; % finding the row number of uth element
    j = mod(u,16) + 16 *(mod(u,16)==0) ; % finding the coloumn of the uth element
   for v = 1 : 256 % vth element of the stimulus
       k = (ceil(v / 16)) ; % finding the row number of vth element
       l = mod(v,16) + 16 *(mod(v,16)==0) ; % finding the coloumn of the vth element
       correlation_matrix(u,v) = sum(Spike_Triggered(j,i,:) .* Spike_Triggered (l,k,:))/ t ; % multiplying (i,j) and (k,l) through the third dimension and averaging 
   end
end
%% Showing the Image
[V , D] = eig(correlation_matrix) ; % Extravtiong Eigenvalues and Eigenvectors
z1 = reshape(V (:,256), 16,[]); % Extracting most significant eigenvector
z2 = reshape(V (:,255), 16,[]); % Extracting second most significant Eigenvector
z3 = reshape(V (:,254), 16,[]); % Extracting nth Eigenvector
z1 = imresize(z1, 20); % Magnifying the first image
figure
subplot (1,3,1) % Subplot of first image
a = imshowpair(z1,histeq(z1),'montage'); % Smoothing and increasing the contrast
a.CData =  a.CData (1:320,1:320) ; % Extracting the final image
xlabel ('spatial') 
ylabel ('Temporal')
z2 = imresize(z2, 20); % Magnifying the second image
subplot (1,3,2) % Subplot of second image
a = imshowpair(z2,histeq(z2),'montage'); % Smoothing and increasing the contrast
a.CData =  a.CData (1:320,1:320) ; % Extracting the final image
xlabel ('spatial')
ylabel ('Temporal')
z3 = imresize(z3, 20); % Magnifying the third image
subplot (1,3,3) % Subplot of third image
a = imshowpair(z3,histeq(z3),'montage'); % Smoothing and increasing the contrast
a.CData =  a.CData (1:320,1:320) ; % Extracting the final image
xlabel ('spatial')
ylabel ('Temporal')

%% Creating Control Matrix & Control Confidence Level
 Control_Matrix = zeros (256 , 256) ; % Initializing the control matrix
t = size (control_matrix3D,3) ;
for u = 1 : 256 % uth element of control matrix 
    i = (ceil(u / 16)) ; % finding the row number of uth element
    j = mod(u,16) + 16 *(mod(u,16)==0) ; % finding the coloumn number of uth element
   for v = 1 : 256 % vth element of control matrix 
       k = (ceil(v / 16)) ; % finding the row number of vth element
       l = mod(v,16) + 16 *(mod(v,16)==0) ; % finding the coloumn number of vth element
       Control_Matrix(u,v) = sum(control_matrix3D(j,i,:) .* control_matrix3D (l,k,:))/ t ; % multiplying (i,j) and (k,l) through the third dimension and averaging 
   end
end
rank = 1:1:30 ; % creating rank vector_horizontal axis of the eigenvalue plot
Control_Eigenvalue = eig(Control_Matrix) ; % Control matrix eigenvalues
Correlation_Eiegenvalue = eig(correlation_matrix) ; % Correlation matrix eigenvalues
figure
stem (fliplr(rank) , Control_Eigenvalue(227:256)) ; % Steming the Control Matrix eigenvalues
SD =  2 * sqrt(var(Control_Eigenvalue)) ; % Calculating the Standard deviation of Control Matrix eigenvalues
Confidence_Interval_up = Control_Eigenvalue(227:256) + SD ; % Specifying the upper limit of confidence zone
Confidence_Interval_down = Control_Eigenvalue(227:256) - SD ; %Specifying the lower limit of confidence zone
xlim ([1 30])
ylim ([0 3])
hold on 
stem (fliplr(rank) , Correlation_Eiegenvalue(227:256)) ;
plot (fliplr(rank) , Confidence_Interval_up , '--')
plot (fliplr(rank) , Confidence_Interval_down , '--')
legend('Control Matrix Eigenvalue','Correlation Matrix Eigenvalue','Confidence Inteval (Upper Limit)','Confidence Interval (Lower Limit)')
xlabel ('Rank')
ylabel ('Eigenvalue')
hold off
%% Histogram 3D
corr1 = [] ;
corr2 = [] ;
corrcont1 = [];
corrcont2 = [];
Eigen_Matrix1 = reshape(V (:,256), 16,[]) ; % Reshaping the first Eigenvector
Eigen_Matrix2 = reshape(V (:,255), 16,[]) ; % Reshaping the second Eigenvector
for j = 1:length (Spike_Triggered)
    corr_add1 = sum(sum((Eigen_Matrix1).*(Spike_Triggered ( : , : , j ) ))) ; % Projecting spike tiggered stimulus on the first eigenvector
    corr1 = [corr1 corr_add1] ;
end
for j = 1:length (Spike_Triggered)
    corr_add2 = sum(sum((Eigen_Matrix2).*(Spike_Triggered ( : , : , j ) ))) ; % Projecting spike tiggered stimulus on the first eigenvector
    corr2 = [corr2 corr_add2] ;
end
for k = 1 : length (control_matrix3D)
   corrcont_add1 =  sum(sum(Eigen_Matrix1 .* control_matrix3D(: , : , k))) ; % Projecting random stimulus on the first eigenvector
   corrcont1 = [corrcont1 corrcont_add1] ;
end
for k = 1 : length (control_matrix3D)
   corrcont_add2 =  sum(sum(Eigen_Matrix2 .* control_matrix3D(: , : , k))) ; % Projecting random stimulus on the second eigenvector
   corrcont2 = [corrcont2 corrcont_add2] ;
end
figure
histogram2(corr1,corr2) % Plotting Histogram of " projections of spike triggered stimulus on the first and second eigenvector "
hold on
histogram2(corrcont1,corrcont2) % Plotting Histogram of " projections of random stimulus on the first and second eigenvector "
legend('Spike' , 'Control')

%%
[x, y] = meshgrid(linspace(-5, 5));
sigma1 = sqrt(var (corr1)) ; % Calculating the standard deviation of SpikeTriggered projection on first eigenvalue
sigma2 = sqrt(var (corr2)) ; % Calculating the standard deviation of SpikeTriggered projection on second eigenvalue
mean1 = mean (corr1) ; % Calculating the mean of SpikeTriggered projection on first eigenvalue
mean2 = mean (corr2) ; % Calculating the mean of SpikeTriggered projection on second eigenvalue
rho = regression(corr1,corr2) ; % Calculating "RHO" parameter
constant1 = 1/(2*pi*sigma1*sigma2*sqrt (1-rho^2)) ;
constant2 = -1 / (2*(1-rho^2)) ;
z1 = constant1 * exp(constant2 * ((x-mean1).^2/sigma1.^2+(y-mean2).^2/sigma2.^2-(2*rho.*(x-mean1).*(y-mean2))./sigma1.*sigma2)) ; % Fitting the distribution on a gaussian
sigma1 = sqrt(var (corr1)) ; % Calculating the standard deviation of random projection on first eigenvalue
sigma2 = sqrt(var (corr2)) ; % Calculating the standard deviation of random projection on second eigenvalue
mean1 = mean (corrcont1) ; % Calculating the mean of random projection on first eigenvalue
mean2 = mean (corrcont2) ; % Calculating the mean of random projection on second eigenvalue
rho = regression(corrcont1,corrcont2) ; % Calculating "RHO" parameter
constant1 = 1/(2*pi*sigma1*sigma2*sqrt (1-rho^2)) ;
constant2 = -1 / (2*(1-rho^2)) ;
z2 = constant1 * exp(constant2 * ((x-mean1).^2/sigma1.^2+(y-mean2).^2/sigma2.^2-(2*rho.*(x-mean1).*(y-mean2))./sigma1.*sigma2)) ; % Fitting the distribution on a gaussian
zdiff = z1 - z2 ; % The difference of two distribution  
figure
contour(x,y,zdiff,'ShowText','on') % Finding the contour of the difference of the Spike triggered from random stimulus
%%
x = corr1 ; 
y = corr2 ;
sigma1 = sqrt(var (corr1)) ; % Calculating the standard deviation of SpikeTriggered projection on first eigenvalue
sigma2 = sqrt(var (corr2)) ; % Calculating the standard deviation of SpikeTriggered projection on second eigenvalue
mean1 = mean (corr1) ; % Calculating the mean of SpikeTriggered projection on first eigenvalue
mean2 = mean (corr2) ; % Calculating the mean of SpikeTriggered projection on second eigenvalue
rho = regression(corr1,corr2) ; % Calculating "RHO" parameter
constant1 = 1/(2*pi*sigma1*sigma2*sqrt (1-rho^2)) ;
constant2 = -1 / (2*(1-rho^2)) ;
z1 = constant1 * exp(constant2 * ((x-mean1).^2/sigma1.^2+(y-mean2).^2/sigma2.^2-(2*rho.*(x-mean1).*(y-mean2))./sigma1.*sigma2)) ; % Fitting the distribution on a gaussian
sigma1 = sqrt(var (corr1)) ; % Calculating the standard deviation of random projection on first eigenvalue
sigma2 = sqrt(var (corr2)) ; % Calculating the standard deviation of random projection on second eigenvalue
mean1 = mean (corrcont1) ; % Calculating the mean of random projection on first eigenvalue
mean2 = mean (corrcont2) ; % Calculating the mean of random projection on second eigenvalue
rho = regression(corrcont1,corrcont2) ; % Calculating "RHO" parameter
constant1 = 1/(2*pi*sigma1*sigma2*sqrt (1-rho^2)) ;
constant2 = -1 / (2*(1-rho^2)) ;
z2 = constant1 * exp(constant2 * ((x-mean1).^2/sigma1.^2+(y-mean2).^2/sigma2.^2-(2*rho.*(x-mean1).*(y-mean2))./sigma1.*sigma2)) ; % Fitting the distribution on a gaussian
zdiff = z1 - z2 ; % The difference of two distribution 
PCA_Success_Percentage = length (find(zdiff>0)) / length (corr1)  % Percentage of Spikes