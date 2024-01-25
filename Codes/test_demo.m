clear;clc;
warning off 
bits = [16];  

%%   OfficeHome    
db_name   =   'OfficeHome_PtoA';
load ./demo_data/OfficeHome_PtoA_source
load ./demo_data/OfficeHome_PtoA_target 

all_methods = [ "ATH" ];
for met = 1:length(all_methods)                                                                                                                            
    current_method = all_methods(met)  
    for ii = 1:length(bits)
        bits(ii)
        mAP_cross = zeros(10,1);
        mAP_single = zeros(10,1);
        for times = 1:10  
            switch(current_method) 
            case 'ATH' 
                method = 'ATH';
                param.r = bits(ii);                                 
                param.alpha1 = 10^-3;        param.alpha2 = 10^-1; 
                param.beta1 = 10^-2;        param.beta2 = 10^-2;
                param.lambda2 = 10^-2;
                param.T = 10;
                [As, At] = ATH(Xs, Ys, Xt_train, Yt_train, param);
                Bs = (Xs * As > 0);
                Bt_train = (Xt_train * At > 0);
                Bt_test  = (Xt_test * At > 0); 
            otherwise 
                fprintf('Method Not Found!' ); 
            end

           %% cross-domain retrieval
            Bs = compactbit(Bs);
            Bt_test = compactbit(Bt_test);	       
            WtrueTestTraining = compute_S(Ys, Yt_test) ;
            Dhamm = hammingDist(Bt_test, Bs);
            [recall, precision, ~] = recall_precision(WtrueTestTraining', Dhamm);
            mAP_cross(times) = area_RP(recall, precision);
            clear  mAP   Dhamm

           %% single-domain retrieval
            Bt_train = compactbit(Bt_train);      
            WtrueTestTraining = compute_S(Yt_train, Yt_test) ;
            Dhamm = hammingDist(Bt_test, Bt_train);
            [recall, precision, ~] = recall_precision(WtrueTestTraining', Dhamm);
            mAP_single(times) = area_RP(recall, precision);
            clear  mAP  Dhamm
        end
        mAP_cross = sum(mAP_cross) / times;
        mAP_single = sum(mAP_single) / times;
        save(['demo_result/', method, '_', db_name, '_', num2str(bits(ii)), '_bits.mat'], 'mAP_cross', 'mAP_single');
    end
end

