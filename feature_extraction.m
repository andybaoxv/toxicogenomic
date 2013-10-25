clear all,close all,clc;

%% Step 1: Load data
tmp = load('data.mat','data');
data = tmp.data;
Y = data(16,:)';

%% Step 2: Plot the sample ACF and PACF
% figure
% subplot(2,1,1),autocorr(Y)
% subplot(2,1,2),parcorr(Y)

%% Step 3: Fit ARMA(p,q) models
LOGL = zeros(4,4);
PQ = zeros(4,4);
for p = 1:4
    for q = 1:4
        mod = arima(p,0,q);
        [fit,~,logL] = estimate(mod,Y,'print',false);
        LOGL(p,q) = logL;
        PQ(p,q) = p+q;
    end
end

%% Step 4: Calculate the BIC
LOGL = reshape(LOGL,16,1);
PQ = reshape(PQ,16,1);
[aic,bic] = aicbic(LOGL,PQ+1,25);
aic = reshape(aic,4,4);
[row_aic,col_aic] = find(aic==min(min(aic)));
bic = reshape(bic,4,4);
[row_bic,col_bic] = find(bic==min(min(bic)))
