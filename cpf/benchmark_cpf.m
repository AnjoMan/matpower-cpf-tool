base = loadcase('case118mod.mat');


repeat = 20;



times = zeros(1, repeat);

hWaitbar = waitbar(0, 'Benchmarking CPF');

for iter = 1:repeat,
	waitbar(iter/repeat, hWaitbar);
	
	tic;
	
	result = cpf(base, -1,10, 0.025);
	
	times(iter) = toc;
	
end

pause(0.2); delete(hWaitbar);

fprintf('Average completion time: %f\n', mean(times));