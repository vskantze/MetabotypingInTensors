function cell_equal_solutions = EqualSolutions(solution_cell)
%EqualSolutions - checks which cluster solutions are the same.
%
% Syntax:  cell_equal_solutions = EqualSolutions(solution_cell)
%
% Inputs:
%    solution_cell - a cell with stacked cluster solutions
%
% Outputs:
%    cell_equal_solutions - nested cell with the indeces of the equal solutions              
%
% Example: 
%    cell_equal_solutions = EqualSolutions(solution_cell)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required:
% CSV-files required: 
%
% See also: 
% Author: Viktor Skantze
% email: viktor.skantze@fcc.chalmers.se
% August 2020
%------------- BEGIN CODE --------------                                                         
num_solutions_to_compare = length(solution_cell);                          % Number of solutions to compare
list_solutions_left_to_compare = 1:num_solutions_to_compare;               % Prepare to loop and compare over other solutions
cell_equal_solutions = cell(num_solutions_to_compare,1);
indeces_to_remove = [];
for i_solution = 1:num_solutions_to_compare                                % Loop over all solutions to be compared
    solution_current = solution_cell{i_solution};
    solution_current = SortCellArrays(solution_current);
    
    if ~isempty(indeces_to_remove)
        for i_r = indeces_to_remove                                        % remove already taken indeces
            list_solutions_left_to_compare(list_solutions_left_to_compare...
                ==i_r) = [];
        end
    end
    list_solutions_left_to_compare(list_solutions_left_to_compare...
        ==i_solution) = [];                                                % Remove current solution from list so that you don't compare on itself
    indeces_to_remove = [];
    for i_comparing_sol = list_solutions_left_to_compare                   % Loop over the other solutions to compare with current
        solution_compare = solution_cell{i_comparing_sol};
        if length(solution_compare) == length(solution_current)
            solution_compare = SortCellArrays(solution_compare);
            score = cellfun(@isequal,solution_compare,solution_current);
            if sum(score)== length(solution_compare)                       % if they are the same, then save the indeces in the group with the same solutions
                cell_equal_solutions{i_solution} = unique([cell_equal_solutions{i_solution} i_solution i_comparing_sol ]);
                indeces_to_remove = [indeces_to_remove i_comparing_sol];
            end
        else
            continue
        end 
    end
end

cell_equal_solutions = cell_equal_solutions(~cellfun('isempty',cell_equal_solutions)); % output the nested cell with equal solutions


end

function out = SortCellArrays(solution_current)                            % Sorting the nested cell arrays on length
    lengths = cellfun(@length, solution_current);
    [~, idx] = sort(lengths);
    solution_current = solution_current(idx);
    
    new_indeces = 1:length(solution_current);
    previous = lengths(1);
    for i = 2:length(solution_current)
        current = lengths(i);
        if current == previous
           if sum(solution_current{i}.^2) < sum(solution_current{i-1}.^2)
               new_indeces(i-1) = i;
               new_indeces(i) = i-1;
           end
        end
        previous = idx(i);
    end

    out = solution_current(new_indeces);
end
%------------- END CODE -------------- 