SVM_THETA_HOME=pwd; 
setenv('SVM_THETA_HOME', SVM_THETA_HOME); 


%% Compiling libsvm 
disp('Compiling libsvm'); 
cd([SVM_THETA_HOME '/lib/libsvm-3.1']); 
[status, result] = system('make') ;
if status == 0
    disp('compiled libsvm successfully'); 
else    
    error('Could not compile libsvm'); 
end
setenv('LIBSVM_ROOT', [SVM_THETA_HOME '/lib/libsvm-3.1']); 
cd(SVM_THETA_HOME);

%% Compiling vjlibsvmapi 
disp('Compiling vjlibsvmapi'); 
cd([SVM_THETA_HOME '/src']); 
setenv('MATLAB_HOME', matlabroot); 
[status, result] = system('make');
if status == 0
    disp('compiled vjlibsvmapi successfully'); 
else    
    error('Could not compile vjlibsvmapi'); 
end
cd(SVM_THETA_HOME);

%% Add the paths
disp('Adding the paths'); 
addpath(genpath([SVM_THETA_HOME '/lib'])); 
addpath( [SVM_THETA_HOME '/src'] ); 
warning('on', 'lmb:verbose'); 