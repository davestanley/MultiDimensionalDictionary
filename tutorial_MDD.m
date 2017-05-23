% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % MDD Tutorial % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%% % % % % % % % % % % % % % % MDD Setup % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%% Set up paths 
% Get ready...

% Format
format compact

% Check if in MDD folder
if ~exist(fullfile('.','sample_data.mat'), 'file')
    error('Should be in MDD folder to run this code.')
end


% Add MDD toolbox to Matlab path if needed
if ~exist('MDD','class')
  addpath(genpath(pwd));
end


%% Load some sample data

% Load some sample simulated data
load('sample_data.mat');

% We loaded 3 variables: dat, axis_values, and axis_names. dat contains a
% 4D cell array containing time series data. It represents several
% simulations of a network of coupled excitatory (E) and
% inhibitory (I) cells.

% Axis_names gives the names of what is varied along each axis of dat.
% These are:
% 1) Applied current to E cells (E_Iapp);
% 2) Inhibitory synapse decay constant Tau (I_E_tauD);
% 3) Population name (excitatory or inhibitory cells - E or I);
% 4) Name of state variable
% These are stored in axis_names:
fprintf('axis_names = ')
disp(axis_names)

% The possile values that each of these axes can take on are listed in
% axis_vals. For example, the population axis, axis 3, can be either E or
%  I for excitatory or inhibitory cells respectively.
fprintf('axis_vals{3} = ')
disp(axis_vals{3}')

% Note that the number of entries in axis_vals must 1:1 match up with the 
% size of dat.
fprintf('axis_vals = ')
disp(axis_vals);
fprintf('size(dat) = ')
disp(size(dat));


% Thus, dat is of the form (E_Iapp, I_E_tauD, population, variable).
% For example:
figure; plot(dat{1,1,2,1}); title('Inhibitory (I) cell voltage'); ylabel('Vm'); xlabel('Time (ms)');
fprintf('E_Iapp parameter value = ')
disp(axis_vals{1}(1))       % E_Iapp parameter value

fprintf('I_E_tauD parameter value = ')
disp(axis_vals{2}(1))       % I_E_tauD parameter value

fprintf('Inhibitory cells label = ')
disp(axis_vals{3}{2})       % Inhibitory cells label

fprintf('Voltage label = ')
disp(axis_vals{4}{1})       % Voltage label


%% Import into MDD object

% All of this information can be imported into an MDD object.

% Create MDD object
xp = MDD;

% Import the data
xp = xp.importData(dat,axis_vals,axis_names);


% MDD objects are essentially cell arrays (or matricies), but with the
% option to index using strings instead of just integers. 
% Thus, they are analogous to dictionaries in Python.
% (This core functionality is implemented by the multidimensional
% dictionaries (MDD), which MDD inherits, and to which MDD adds
% plotting functionality.)
disp(xp);


% At its core, MDD has 3 fields. xp.data stores the actual data (either a 
% matrix or a cell array). Note: it is okay that there are some empties.
disp(xp.data);
size(xp.data)

% Next, xp.axis stores axis labels associated with each dimension in
% xp.data.
disp(xp.axis(1));

% Axis.values stores the actual axis labels. These can be numeric...
disp(xp.axis(1).values);

% ...or string type. As we shall see below, these axis labels can be
% referenced via index or regular expression.
disp(xp.axis(4).values);

% Axis.name field stores the name of the dimension. In this case, it is
% named after the parameter in the model that was varied.
disp(xp.axis(1).name)

% Axis.axismeta is for internal use and is currently empty.
xp.axis(1).axismeta

% All of this information can be obtained in summary form by running
% printAxisInfo
xp.printAxisInfo

% xp.meta stores meta data for use by the user as they see fit.
% Here we will add some custom info to xp.metadata. This can be whatever
% you want. Here, I will use this to provide information about what is
% stored in each of the matrices in the xp.data cell array. (Alternatively,
% we could also make each of these matrices an MDD object!)
meta = struct;
meta.datainfo(1:2) = MDDAxis;
meta.datainfo(1).name = 'time(ms)';
meta.datainfo(1).values = time;
meta.datainfo(2).name = 'cells';
meta.datainfo(2).values = [];
xp.meta = meta;
clear meta


%% % % % % % % % % % % % % % % MDD BASICS % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%% Importing data

% MDD objects are just matrices or cell arrays with extra functionality to
% keep track of what each dimension is storing. 

% Above, we imported data into an MDD object along with axis values and names. 
% But it is also possible to import a high dimensional matrix alone.

% For example let's say we have a random cell array.
mydata = dat;

% We can import the data 3 different ways:
%   1) Calling the public method, importData
xp2 = MDD;
xp2 = xp2.importData(mydata);

%   2) Using a class/static method without having to make the object first by using
%      an uppercase method name.
xp2_alt = MDD.ImportData(mydata);

%   3) Calling the class constructor directly
xp2_alt2 = MDD(mydata);

% The resulting objects are equivalent
isequal(xp2, xp2_alt, xp2_alt2)

%% Importing data with axis values
% We didn't supply any axis names/values, so default values were assigned
xp2.printAxisInfo;

% We can instead import axis values along with the data
xp2 = MDD;
xp2 = xp2.importData(mydata, axis_vals);
xp2.printAxisInfo

% This double argument interface works for the other two methods as well:
xp2_alt = MDD.ImportData(mydata, axis_vals);
xp2_alt2 = MDD(mydata, axis_vals);

% The resulting objects are again equivalent
isequal(xp2, xp2_alt, xp2_alt2)

% These axis values can be acquired from an object using the exportAxisVals method
ax_vals = xp2.exportAxisVals;

%% Importing data with axis names

% Axis names can be assigned in this way as well
xp2 = xp2.importData(mydata, axis_vals, axis_names);
xp2.printAxisInfo

% This triple argument interface also works with the other 2 approaches:
xp2_alt = MDD.ImportData(mydata, axis_vals, axis_names);
xp2_alt2 = MDD(mydata, axis_vals, axis_names);

% The resulting objects are yet again equivalent
isequal(xp2, xp2_alt, xp2_alt2)

% The axis names are accessible from an object via the exportAxisNames method
ax_names = xp2.exportAxisNames;

clear xp2_alt xp2_alt2

%% Exporting Data to 2D Table

% Multi-dimensional data is often represented in 2D table form, with 1 column
% representing the data, and other columns representing parameters associated
% with the data. The multidimensional data from an MDD object can be exported
% to a 2D table as below.

[data_column, axis_val_columns, axis_names] = xp.exportDataTable();

% The table can be previewed in the command window by calling:
xp.exportDataTable(true);

% The number of rows printed to screen can be changed from the default of 10 by
% passing a second argument.
xp.exportDataTable(true, 5); % printing 5 rows to screen


%% Importing Data from 2D Table

% As mentioned before, multi-dimensional data is often represented in 2D table form, 
% with 1 column representing the data, and other columns representing parameters 
% associated with the data. MDD can import this 2D data, as we generated
% previously.

% First let us inspect the tabular data sizes:

fprintf('Data column size: %s\n', num2str(size(data_column)))
fprintf('Axis values columns size: %s\n', num2str(size(axis_val_columns)))
fprintf('Axis names size: %s\n', num2str(size(axis_names)))

% As with importData, there are 2 interfaces for importing:
xp3 = MDD;
xp3 = xp3.importDataTable(data_column, axis_val_columns, axis_names); % lowercase object method
xp3.printAxisInfo

%  or

xp3 = MDD.ImportDataTable(data_column, axis_val_columns, axis_names); % uppercase class method
xp3.printAxisInfo


%% MDD Indexing

% Indexing works just like with normal matrices and cell arrays and axis
% labels are updated appropriately.
clc
xp4 = xp(:,:,1,8);                  % ## Update - This does the same thing as xp.subset([],[],1,8), which was the old way
                                    % of pulling data subsets. Note that [] and : are equivalent.
xp4.printAxisInfo

% Similarly, can index axis values using substring matching (via strfind internally)
% Pull out sodium mechs only
xp5 = xp(:,:,1,'iNa');
xp5.printAxisInfo

% Pull out synaptic state variables
xp5 = xp(:,:,1,'_s');
xp5.printAxisInfo

% Same as before, but using regular expression syntax:
%   '/regularExpressionString/' will get passed to regexp as 'regularExpressionString' (ie without the enclosing forward slashes)
xp5 = xp(:,:,1,'/_s$/');
xp5.printAxisInfo

% Can also reference a given axis based on its index number or based on its
% name
disp(xp.axis(4))
disp(xp.axis('populations'))

% Lastly, you can reference xp.data with the following shorthand
% (This is the same as xp.data(:,:,1,8). Regular expressions dont work in this mode)
mydata = xp{:,:,1,8};
mydata2 = xp.data(:,:,1,8);
disp(isequal(mydata,mydata2));

clear mydata mydata2 xp4 xp5


%% % % % % % % % % % % % % % % PLOTTING EXAMPLES % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%% Plot 2D data
% Tip: don't try to understand what recursiveFunc is doing - instead, try
% putting break points in the various function handles to see how this
% command works.
close all;

% Pull out a 2D subset of the data
clc
xp4 = xp(:,:,'E','v');
xp4.printAxisInfo

% Set up plotting arguments
function_handles = {@xp_subplot_grid,@xp_matrix_basicplot};   % Specifies the handles of the plotting functions
dimensions = {{'E_Iapp','I_E_tauD'},{'data'}};                % Specifies which axes of xp each function handle
                                                                % will operate on. Note that dimension 'data' refers to the 
                                                                % the contents of each element in xp.data (e.g. the matrix of
                                                                % time series data). It must come last.
function_arguments = {{},{}};	% This allows you to supply input arguments to each of the 
                                % functions in function handles. For
                                % now we'll leave this empty.
                                                                
% Run the plot. Note the "+" icons next to each plot allow zooming. 
figl; recursiveFunc(xp4,function_handles,dimensions,function_arguments);

%% Plot 3D data 

close all;

% Pull out a 3D subset of data (parameter sweeps and the 2 cell
% types)
clc
xp4 = xp(:,1:2,:,'v');
xp4.printAxisInfo

% This will plot E cells and I cells (axis 3) each in separate figures and
% the parameter sweeps (axes 1 and 2) as subplots.
dimensions = {{'populations'},{'I_E_tauD','E_Iapp'},{'data'}};
recursiveFunc(xp4,{@xp_handles_newfig,@xp_subplot_grid,@xp_matrix_imagesc},dimensions);

% Note that here we produced rastergrams instead of time series by
% submitting a different function to operate on dimension zero.

%% Plot 3D data re-ordered

% Alternatively, we can put E and I cells in the same figure. This
% essentially swaps the population and tauD axes.
dimensions = {{'I_E_tauD'},{'populations','E_Iapp'},'data'};
recursiveFunc(xp4,{@xp_handles_newfig,@xp_subplot_grid,@xp_matrix_imagesc},dimensions);

%% Plot 4D data

close all;

% Pull out sodium channel state variables for E and I cells.
clc
xp4 = xp(1:2,1:2,:,6:7);
xp4.printAxisInfo

dimensions = {'populations',{'E_Iapp','I_E_tauD'},'variables',0};       % Note - we can also use a mixture of strings and index locations to specify dimensions. Dimension "0" corresponds to data.

% Note that here we will supply a function argument. This tells the second
% subplot command to write its output to the axis as an RGB image, rather than
% as subplots. This "hack" enables nested subplots.
xp_subplot_grid_options.display_mode = 1;
function_arguments = {{},{},{xp_subplot_grid_options},{}};

if verLessThan('matlab','8.4'); error('This will not work on earlier versions of MATLAB'); end
recursiveFunc(xp4,{@xp_handles_newfig,@xp_subplot_grid,@xp_subplot_grid,@xp_matrix_basicplot},dimensions,function_arguments);

%% Plot multiple dimensions adaptively.

close all;

% Another option is to use @xp_subplot_grid_adaptive, which will plot the data using axes in
% descending order of the size of the axis values, and plot remaining
% combinations of axis values across figures.

recursiveFunc(xp4,{@xp_subplot_grid_adaptive,@xp_matrix_basicplot},{1:4,0});

%% Combine and Plot two MDD objects
close all;
clc
xp3 = xp(2,:,'E','v');
xp3.printAxisInfo

xp4 = xp(:,3,'E','v');
xp4.printAxisInfo

% Notice that xp3 and xp4 are overlapping at
% Axis 1: E_Iapp (numeric) -> 10
% Axis 2: I_E_tauD (numeric) -> 15

% Attempt to merge them
xp5 = merge(xp3,xp4); % or xp5 = xp3.merge(xp4);

% This throws a warning that there is an overlap, and sets xp5 = xp3
% We will disregard the message by setting the third argument to true, allowing 
% xp4 to overwrite xp3.
xp5 = merge(xp3,xp4, true); % or xp5 = xp3.merge(xp4, true);

dimensions = {[1,2],0};
figl; recursiveFunc(xp5,{@xp_subplot_grid,@xp_matrix_imagesc},dimensions);


%% % % % % % % % % % % % % % % ADVANCED MDD / MDD USAGE % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%% Modifying MDD.data directly

% While it is encouraged to use importData, MDD.data can also be written
% to directly. For example, we can take the average across all cells by
% doing the following.
xp2 = xp;
for i = 1:numel(xp2.data)
    xp2.data{i} = mean(xp2.data{i},2);
end
disp(xp2.data);         % (Note the size is 10001x1 instead of by 10001x20 or 10001x80)

% However, you cannot make modifications that would destroy the 1:1 matchup
% between data and axis (commenting out error producing commands below).
xp2 = xp;
xp2.data{5} = randn(100);       % Ok
% xp2.axis = xp2.axis(1:2);       % Not ok
% mydata = reshape(xp2.data,[3,3,16]);
% xp2.data = mydata;              % Not ok.

%% Method packDim
% Analogous to cell2mat.
clear xp2 xp3 xp4 xp5

% Start by taking a smaller subset of the original xp object.
% xp2 = xp.subset(2,2,[],[1,3,5:8]);      % Selection based on index locations
xp2 = xp.subset(2,2,:,'/(v|^i||ISYN$)/');  % Same thing as above using regular expression. Selects everything except the _s terms. "^" - beginning with; "$" - ending with
xp2 = xp2.squeeze;
xp2.printAxisInfo;

% Note that xp2 is sparse (there are some empty cells)
disp(xp2.data);         % (E cells don't receive AMPA synapses, and I cells don't receive GABAA synapses)

% Now pack dimension two (columns) of xp2 into xp2.data.
src = 2;                    % Take 2nd dimension in xp2
dest = 3;                   % Pack into 3rd dimension in xp2.data matrix
xp3 = xp2.packDim(src,dest);

% Check dimensionality of xp3.data
disp(xp3.data)             % The dimension "variables", which was dimension 2
                           % in xp2, is now dimension 3 in xp3.data.
                           % Now xp3.data is time x cells x variables

% View axis of xp3
xp3.printAxisInfo;            % The dimension "variables" is now missing

% Alternatively, you can use a regular expression to select the dimension
% you want to pack; if the destination dimension is left empty, the first
% dimension which is not occupied in any of the matrix entries of xp.data
% will be used.
xp3 = xp2.packDim('var');
xp3.printAxisInfo;

% Note some of this data is sparse! We can see this sparseness by plotting
% as follows (note the NaNs)
temp1 = squeeze(xp3.data{1}(100,:,:));  % Pick out a random time point
temp2 = squeeze(xp3.data{2}(100,:,:));  % Pick out a random time point
figure; 
subplot(211); imagesc(temp1);
ylabel('Cells');
xlabel(xp2.axis(2).name); 
set(gca,'XTick',1:length(xp2.axis(2).values)); set(gca,'XTickLabel',strrep(xp2.axis(2).values,'_',' '));

subplot(212); imagesc(temp2);
ylabel('Cells');
xlabel(xp2.axis(2).name); 
set(gca,'XTick',1:length(xp2.axis(2).values)); set(gca,'XTickLabel',strrep(xp2.axis(2).values,'_',' '));

%% Method unpackDim (undoing packDims)
% When packDim is applied to an MDD object, say to pack dimension 3, the
% information from the packed axis is stored in the MDDAxis
% matrix_dim_3, a field of xp3.meta.
xp3.meta.matrix_dim_3.printAxisInfo

% If dimension 3 of each cell in xp3.data is unpacked using unpackDim,
% xp3.meta.matrix_dim_3 will be used to provide axis info for the new
% MDD object.
xp4 = xp3.unpackDim(dest, src);
xp4.printAxisInfo;

% Unless new axis info is provided, that is.
xp4 = xp3.unpackDim(dest, src, 'New_Axis_Names'); % The values can also be left empty, as in xp4 = xp3.unpackDim(dest, src, 'New_Axis_Names', []);
xp4.printAxisInfo;

xp4 = xp3.unpackDim(dest, src, 'New_Axis_Names', {'One','Two','Three','Four','Five','Six'});
xp4.printAxisInfo;

%% Use packDim to average across cells

xp2 = xp;
xp2 = xp(:,:,:,'v');  % Same thing as above using regular expression. Selects everything except the _s terms. "^" - beginning with; "$" - ending with
xp2 = xp2.squeeze;
%
% Average across all cells
xp2.data = cellfun(@(x) mean(x,2), xp2.data,'UniformOutput',0);

% % Convert xp2.data from a matrix into an MDD object as well. This is
% % useful for keeping track of axis names. 
% mat_ax_names = {'Time','Cell Number'};
% mat_ax_values = {1:10001, []};
% 
% % xp2.data = Cell_2_MDD(xp2.data,mat_ax_names,mat_ax_values);

% Pack E and I cells together
src=3;
dest=2;
xp3 = xp2.packDim(src,dest);


% Plot 
figl; recursiveFunc(xp3,{@xp_subplot_grid,@xp_matrix_basicplot},{[1,2],[]},{{},{}});

%% Use packDim to average over synaptic currents
% Analogous to cell2mat
% See also plotting material by Hadley Wickham

% First, pull out synaptic current variables
xp2 = xp(:,:,:,'/(ISYN$)/');  % Same thing as above using regular expression. Selects everything except the _s terms. "^" - beginning with; "$" - ending with
xp2.printAxisInfo;

% Second, put this into matrix form, so we can average over them
xp3 = xp2.packDim(4,3);
disp(xp3.data)              % xp3.data is now 3D, with the 3rd dim denoting synaptic current
xp3 = xp3.squeeze;
xp3.printAxisInfo;

% Average across membrane currents
xp3.data = cellfun(@(x) nanmean(x,3), xp3.data,'UniformOutput',0);

% Plot 
recursiveFunc(xp3,{@xp_handles_newfig,@xp_subplot_grid,@xp_matrix_basicplot},{[3],[1,2],[0]});

%% Test mergeDims
% Analogous to Reshape.

% This command combines two (or more) dimensions into a single dimension.
xp2 = xp.mergeDims([3,4]);
xp2.printAxisInfo;


%% Advanced testing
clear xp2 xp3 xp4 xp5 xp6
% Test squeezeRegexp
xp2 = xp(:,1,:,end); xp2.printAxisInfo

xp2b = xp2.squeezeRegexp('var'); xp2b.printAxisInfo
xp2b = xp2.squeezeRegexp('I_E_tauD'); xp2b.printAxisInfo
xp2b = xp2.squeezeRegexp('populations'); xp2b.printAxisInfo

%% To implement
% 
% Implement the following:
% + Make DynaSimPlotExtract more general
% + Starting work on dsPlot2 - any new requests?
