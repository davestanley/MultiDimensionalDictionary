# MDD
MATLAB Multidimensional Dictionary class. MDD allows the implementation of multidimensional Python-esque dictionaries in MATLAB.

## In a nutshell
MDD can be interpreted in several different ways. 
- A map/dictionary that associates multiple keys with a single value
- An N-dimensional table (a table is equivalent to an MDD object in 2-dimensions)
- A matrix or cell array that can be indexed by using strings and regular expressions

## Overview
MDD is a MATLAB tool for managing high-dimensional data that often arises in scientific data analysis. It extends the core functionality of Matlab cells and matrices to advanced labeling and indexing options. 

## Examples
MDD allows you to do things like this:
- Initialize an MDD object with axis labels:
```
>> mdd = MDD(stocks_data,axis_labels);
```
```
mdd = 
                        2017         2018         2019  
--------------------------------------------------------
'Tech_Apple'        |    120.1       175.3       191.3
'Tech_Microsoft'    |     62.6        90.3       117.2
'Retail_CVS'        |     82.2        78.4        56.4
'Retail_Walmart'    |     68.3       100.1        98.2

```
- Index using numerical queries
```
>> mdd2 = mdd(:,'2017 < x < 2019');
```
```
mdd2 = 
                        2018  
--------------------------------------------------------
'Tech_Apple'        |    175.3
'Tech_Microsoft'    |     90.3
'Retail_CVS'        |     78.4
'Retail_Walmart'    |    100.1
```
- Index using regular expressions
```
>> mdd3 = mdd('Tech_',:);
```
```
mdd3 = 
                        2017         2018         2019  
--------------------------------------------------------
'Tech_Apple'        |    120.1       175.3       191.3
'Tech_Microsoft'    |     62.6        90.3       117.2

```

## Details
An MDD object can be thoughout of as a MATLAB cell array (or matrix) with some additional functionality. At its core, it extends the way cells and matrices are indexed by allowing string labels to be assigned to each dimension of the cell array, similar to how row and column names are assigned to a table (e.g., in Pandas or SQL). MDD objects can then be indexed, sorted, merged, and manipulated according to these labels. (See section: MDD subscripts and indexing)

Additionally, MDD includes methods for performing operations on high dimensional data. The goal is to modularize the process of working with high dimensional data. Within MDD, functions designed to work on low dimensional data (1 or 2 dimension) can each be assigned to each operate on different dimensions of a higher dimensional object. Chaining several of these functions together can allow the entire high dimensional object to be processed. The advantage of this modular approach is that functions can be easily assigned other dimensions or swapped out entirely, without necessitating substantial code re-writes. (See section: Running functions on MDD objects)

## Demo
http://htmlpreview.github.com?https://github.com/davestanley/MultiDimensionalDictionary/blob/master/html/demo_MDD.html

**WARNING**: This is a work in progress. The core code is implemented, but documentation is missing. I will be updating it more over the coming months.

This code is similar to the multidimensional map function implemented by David Young. However, this implementation does not use MATLAB maps and instead is adds functionality to traditional MATLAB matrices and cell arrays.



**Related commands**: 
- MATLAB Map Containers (https://www.mathworks.com/help/matlab/map-containers.html)
- Multidimensional implementation of this by David Young (https://www.mathworks.com/matlabcentral/fileexchange/33068-a-multidimensional-map-class)
