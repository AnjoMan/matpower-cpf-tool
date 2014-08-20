function test_cpf()% close all; clear all; clc
    close all;
    
    base = loadcase('case118_mod');
%     base = loadcase('case30_mod');
%     faults = defineFaults(base,2);
%     mFault = faults{20509};


    mFault = Fault('dual', {[8],[],[],[]});
    myCases = mFault.applyto(base, false, false);


    myCase = myCases{1};
%     myCase = myCases{2};


    out = cpf(myCase,-1,1,0.025,true, true, true)
    
%     out = cpf(base, 7, 5, 0.025, true, true)
%     [max_lambda, predicted, corrected, combined] = cpf(myCase, -1, 5,0.025,true);
    
    
    

    %%time test:
%     times = [];
%     fprintf('__________|\n');
%     indices = randi([1,length(faults)], 10,1);
%     for i = indices(:)'
%        mFault = faults{i};
%        
%        myCase = mFault.applyto(base);
%        myCase = myCase{1};
%        
%        tic;
%        out = cpf(myCase,-1,5,0.025);
%        time = toc;
%        
%        times = [times out.time];
%        fprintf('.');
%     end
%     
%     fprintf('\nAvg computation time: %f\n', mean(times));
%     for i = 1:length(indices),
%        fprintf('Fault %d: %f seconds\n', indices(i), times(i)); 
%     end
end