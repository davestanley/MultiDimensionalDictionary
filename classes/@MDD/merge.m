function obj_out = merge(obj1, obj2)
% merge - multidimensional merge of 2 MDD objects without checking for overlap
%
% Usage: obj_out = qmerge(obj1,obj2)
%        obj_out = qmerge(obj1,obj2, forceMergeBool)
%
% Inputs:
%   obj1/2: MDD objects

% Get object properties
Nd1 = ndims(obj1);
axis_vals1 = obj1.exportAxisVals;
axis_names1 = obj1.exportAxisNames;

Nd2 = ndims(obj2);
axis_vals2 = obj2.exportAxisVals;
axis_names2 = obj2.exportAxisNames;

% Find overlapping axis names
[~, indInt1, indInt2] = intersect(axis_names1, axis_names2, 'stable');
[uniqueName1, indUnique1] = setdiff(axis_names1, axis_names2, 'stable');
[uniqueName2, indUnique2] = setdiff(axis_names2, axis_names1, 'stable');

% Find merged number of axes
nAx = length(indInt1) + length(indUnique1) + length(indUnique2);

% Confirm that unique axis names only have 1 value
if ~isempty(indUnique1)
    for k = indUnique1
        if length(axis_vals1{k}) ~= 1
            error('Non-unique axis has more than 1 value')
        end
    end
end
if ~isempty(indUnique2)
    for k = indUnique2
        if length(axis_vals2{k}) ~= 1
            error('Non-unique axis has more than 1 value')
        end
    end
end

% Add axes to front, up to number merged axes
% then shift dim so that front added are on back
if nAx > Nd1
    obj1 = obj1.shiftdim(Nd1-nAx);
    obj1 = obj1.shiftdim(nAx-Nd1);
    obj1.axis(Nd1+1:nAx).name = deal('uniqueAxObj2');
end
if nAx > Nd2
    obj2 = obj2.shiftdim(Nd2-nAx);
    obj2 = obj2.shiftdim(nAx-Nd2);
    obj2.axis(Nd2+1:nAx).name = deal('uniqueAxObj1');
end

% Permute so axes match
permInd = zeros(1, nAx);
permInd(indInt1) = indInt2;
if ~isempty(indUnique2)
    permInd( find(permInd == 0, length(indUnique2), 'last') ) = indUnique2;
end
permInd(permInd == 0) = Nd2+1:nAx;
obj2 = obj2.permute(permInd);

% Add unique axis names and corresponding single values
if nAx > Nd1
    uniqueAxNameInd1 = strcmp(obj1.exportAxisNames, 'uniqueAxObj2');
    [obj1.axis(uniqueAxNameInd1).name] = deal(uniqueName2{:});
    
    uniqueAxNameInd1 = find(uniqueAxNameInd1);
    for k = 1:length(indUnique2)
        obj1.axis(uniqueAxNameInd1(k)).values = obj2.axis(indUnique2(k) == permInd).values;
    end
end
if nAx > Nd2
    uniqueAxNameInd2 = strcmp(obj2.exportAxisNames, 'uniqueAxObj1');
    [obj2.axis(uniqueAxNameInd2).name] = deal(uniqueName1{:});
    
    uniqueAxNameInd2 = find(uniqueAxNameInd2);
    for k = 1:length(indUnique1)
        obj2.axis(uniqueAxNameInd2(k)).values = obj1.axis(indUnique1(k)).values;
    end
end

% Fill in obj_out
obj_out = obj1.reset; % make blank obj
obj_out.axis(1:nAx) = deal(obj1.axisClass); % add empty axes
dataIndObj1 = cell(1,nAx);
dataIndObj2 = cell(1,nAx);
objOutDataSize = zeros(1,nAx);
for iAx = 1:nAx
    % 1) Add names to obj_out
    obj_out.axis_pr(iAx).name = obj1.axis(iAx).name; % add name
    
    % 2) Expand obj_out axes to take on unique values in obj1 and obj2
    if isnumeric(obj1.axis(iAx).values)
        obj_out.axis_pr(iAx).values = unique([obj1.axis(iAx).values(:); obj2.axis(iAx).values(:)], 'sorted'); % add unique values in sorted order
    elseif iscellstr(obj1.axis(iAx).values)
        obj_out.axis_pr(iAx).values = unique([obj1.axis(iAx).values(:); obj2.axis(iAx).values(:)], 'stable'); % add unique values in merged order
    end
    
    % get indicies for data
    [~, dataIndObj1{iAx}] = intersect(obj_out.axis_pr(iAx).values, obj1.axis(iAx).values);
    dataIndObj1{iAx} = dataIndObj1{iAx}';
    [~, dataIndObj2{iAx}] = intersect(obj_out.axis_pr(iAx).values, obj2.axis(iAx).values);
    dataIndObj2{iAx} = dataIndObj2{iAx}';
    
    objOutDataSize(iAx) = length(obj_out.axis_pr(iAx).values);
end

% 3) Fill in data in obj_out from obj1 and obj2
obj_out.data_pr = cell(objOutDataSize);
if isnumeric(obj1.data)
    obj_out.data_pr(dataIndObj1{:}) = num2cell(obj1.data);
else
    obj_out.data_pr(dataIndObj1{:}) = obj1.data;
end
if isnumeric(obj2.data)
    obj_out.data_pr(dataIndObj2{:}) = num2cell(obj2.data);
else
    obj_out.data_pr(dataIndObj2{:}) = obj2.data;
end

% Turn cellnum data into numeric
emptyCells = cellfun(@isempty, obj_out.data_pr);
if ~any(emptyCells(:)) && iscellnum(obj_out.data_pr)
    obj_out.data_pr = cell2mat(obj_out.data_pr);
end

% Combine meta
obj_out = obj_out.importMeta(catstruct(obj1.meta, obj2.meta));

end