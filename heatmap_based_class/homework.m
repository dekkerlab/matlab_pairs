

% the idea behind the homework is to
% use some of discussed MatLab functions
% in order to "simulate" some Hi-C data:
% simluated heatmap should have the following
% features:
%  - checkerboard pattern
%      (aka interactions of hetero and/or euchromatin)
%  - distance decay of interaction density
%      (consequence of tight chromatin packing in 3D)


% hints:
% we'll operate in a binned heatmap level, so
% assume some number of bins: nbins - as many as you like

% of course, index of a bin - is proportional to genomic
% distance

% generate some vector (1D array) with nbins values,
% that decay with index (any decay if fine!)
% try to use toeplitz function to create a 2D matrix
% with generated decay (see the workshop script and "google" for details)


% here is a hint on how to generate some checkerboard in an effecient way:
alternating_pattern = [1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1]

% use outer product of vectors to generate 2D pattern
heatmap(transpose(alternating_pattern) * alternating_pattern)
% read the beginning of https://en.wikipedia.org/wiki/Outer_product or "google"
% think about what's going on

% you've got 2 patterns as a separate matrices/heatmaps/whatever
% how would you combine the 2 to get some "simulated" Hi-C data ?
% try whatever seems appropriate


% show the result using heatmap/imagesc .


