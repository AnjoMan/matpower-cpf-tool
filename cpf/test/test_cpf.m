function test_cpf()% close all; clear all; clc
    close all;
    
    base = loadcase('case30_mod');
    faults = defineFaults(base,1);
    mFault = faults{5};
    myCase = mFault.applyto(base);


    myCase = myCase{1};



    out = cpf(myCase,-1,5,0.025,true,true)

%     out = cpf(base, 7, 5, 0.025, true, true)
%     [max_lambda, predicted, corrected, combined] = cpf(myCase, -1, 5,0.025,true);
    
    
    
end