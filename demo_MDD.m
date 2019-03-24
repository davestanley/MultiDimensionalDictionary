
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % MDD Demo % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


%%

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % MDD Setup % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%% Set up paths and formatting

% Format
format compact
format short g

% Check if in MDD folder
if ~exist(fullfile('.','sample_data.mat'), 'file')
    error('Should be in MDD folder to run this code.')
end


% Add MDD toolbox to Matlab path if needed
if ~exist('MDD','class')
  addpath(genpath(pwd));
end

plot_on = false;


%% Load a sample dataset

% Load some sample simulated data
load('sample_data.mat');
load('sample_data_meta2.mat');
%%
% The file _sample_data.mat_ contains 3 variables: dat, axis_vals, and
% axis_names. In particular, _dat_ is a 4-dimensional cell array:

whos dat

%%
% Each cell in _dat_ cell contains time series output from a 
% neural network simulation. For example, here is the contents of
% _dat(1,1,2,1):
dat(1,1,2,1)

%%
% This is data from 20 simulated interneurons. Here is a plot of this data:

if plot_on
    % Plot one of the cells in dat
    figure; plot(dat{1,1,2,1}); title('Inhibitory (I) cell voltage');
    ylabel('Membrane voltage (mV)'); xlabel('Time (ms)');
end

%%%
% The other variables, _axis_vals_, and _axis_names_ contain metadata about
% what's in _dat_. We'll explore these more later.

%% Import the data into an MDD object
% Let's import this data into an MDD object. There are several ways to
% do this, but we'll do the simplest, which is to call the class constructor.

% Import data into MDD
mdd = MDD(dat,axis_vals,axis_names);


%% What does an MDD object do and why use one?
% When working with high-dimensional data, like that cell array _dat_, it's
% often difficult to keep track of what each axis and entry represents.
% *MDD objects provide a way of organizing and manipulating this 
% information.*

%%% 
% Let's print a summary of the data contained within _mdd_:
mdd.printAxisInfo

%%
% This tells us that _mdd_ is a 4-dimensional object, just like _dat_. The
% size of the axes are [3,3,2,8]. The following lines (Axis1-4) describe each of the
% dimensions of the data. This particular dataset is from a
% series of neural network simulations in which two parameters are varied.
% The values of these parameters (param1 and param2) are given by the first
% two axes. The remaining two axes store data about a particular
% simulation, specifically the different cell types and different variables
% associated with each cell type. 


% xp.meta stores meta data for use by the user as they see fit.
% Here we will add some custom info to xp.meta. This can be whatever
% you want. Here, I will use this to provide information about what is
% stored in each of the matrices in the xp.data cell array. (Alternatively,
% we could also make each of these matrices an MDD object!)
meta = struct;
meta.datainfo(1:2) = MDDAxis;
meta.datainfo(1).name = 'time(ms)';
meta.datainfo(1).values = time;
meta.datainfo(2).name = 'cells';
meta.datainfo(2).values = [];
mdd.meta = meta;
clear meta

%%
% Using some of MDD's functions, we can visualize this data

% Tip: don't try to understand what recursiveFunc is doing - instead, try
% putting break points in the various function handles to see how this
% command works.
close all;

% Pull out a 2D subset of the data
clc
mdd4 = mdd(:,:,'E','v');
mdd4.printAxisInfo

% Set up plotting arguments         (NOTE: subplot_handle = @(xp) xp_subplot_grid(xp,op);)
function_handles = {@xp_subplot_grid,@xp_matrix_basicplot};   % Specifies the handles of the plotting functions
dimensions = {{'param1','param2'},{'data'}};                % Specifies which axes of xp each function handle
                                                                % will operate on. Note that dimension 'data' refers to the 
                                                                % the contents of each element in xp.data (e.g. the matrix of
                                                                % time series data). It must come last.
function_arguments = {{},{}};	% This allows you to supply input arguments to each of the 
                                % functions in function handles. For
                                % now we'll leave this empty.
                                                                
% Run the plot. Note the "+" icons next to each plot allow zooming. 
figl; recursiveFunc(mdd4,function_handles,dimensions,function_arguments);

%%
%%%
% Looking in our _mdd_ object, it contains 3 properties:
mdd

%%
% The property _data_ stores our original 4-D cell array, exactly
% the same as _dat_. The property _meta_ stores metadata optional
% additional metadata about _mdd_ and is currently empty.
% The _axis_ property contain information about each
% axis in our 4D array. For example, in _dat_, the 3rd axis represents the
% type of neural population being simulated.

mdd.axis(3)

%%
% In this particular dataset, there are different neuron populations:
% excitatory (E) cells and inhibitory (I) cells.
% The property _name_ indicates
mdd.axis(3).values{1}
mdd.axis(3).values{2}

% Let's get
% some information about this object.

%%%
% This provides some information about the size of xp and what each of the
% 4 axes represents. 


%% %%%%%%%%%%%%%%%%%%%%%%%%%

%% Load a sample dataset and create an MDD object

% Load some sample simulated data
load('sample_data.mat');

%%%
% The file _sample_data.mat_ contains 3 variables: dat, axis_vals, and
% axis_names. I'll get more into what these contain later but, for now,
% we'll just use these to create an MDD object.

%%% 
% There are several ways to import data, but for now we'll call the class
% constructor, which is the simplest. 

% Import data into MDD
mdd = MDD(dat,axis_vals,axis_names);


%% Explore this object
% Now we have an MDD object called _mdd_. Let's see what it contains.

% Print number of dimensions
ndims(mdd)

%%

% Print size
size(mdd)

%%

% Print size
fprintf('size=%s',mat2str(size(mdd)))

%%
% Therefore, _mdd_ is a 4-dimensional 





