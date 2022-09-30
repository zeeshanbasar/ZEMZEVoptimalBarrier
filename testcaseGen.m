%%% TEST CASE GEN %%%

n = 30; % num. of samples

% set positions
r0 = [normrnd(-800,300,[n/3,1]) normrnd(-800,300,[n/3,1]) normrnd(2500,200,[n/3,1]);...
normrnd(0,350,[n/3,1]) normrnd(0,350,[n/3,1]) normrnd(3000,100,[n/3,1]);...
normrnd(800,300,[n/3,1]) normrnd(800,300,[n/3,1]) normrnd(2000,400,[n/3,1])];


% set velocities
v0 = [normrnd(-25,3,[n/3,1]) normrnd(-25,3,[n/3,1]) normrnd(-80,5,[n/3,1]);...
normrnd(0,10,[n/3,1]) normrnd(0,10,[n/3,1]) normrnd(-80,5,[n/3,1]);...
normrnd(25,3,[n/3,1]) normrnd(25,3,[n/3,1]) normrnd(-80,5,[n/3,1])];


% set mass
mass = [normrnd(1900, 50, [n/2,1]); normrnd(2100, 50, [n/2,1])];


% concat
IC = [r0 v0 mass];


% save
save('init-test_small.mat', 'IC') %% change file name here