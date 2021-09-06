function constraint_matrix = build_constraints(input_constraints)
%BUILD_CONSTRAINTS Automated constrain matrix construction
%  Currently supports simple upper and lower contraints

no_rows     = 2 * size(input_constraints, 1);
no_columns  = size(input_constraints, 1) + 1;

constraint_matrix = zeros(no_rows, no_columns);

for i = 1:size(input_constraints, 1)
    constraint = input_constraints(i, :);
    
    if constraint(1) == 0 && constraint(2) == 0
        mat_1 = [0; 0];
        mat_2 = [-1; -1];
    else
        mat_1 = [-1; 1];
        mat_2 = [constraint(1); -constraint(2)];
    end
    
    constraint_matrix(2*(i-1)+1:2*i, i)     = mat_1;
    constraint_matrix(2*(i-1)+1:2*i, end)   = mat_2;
    
end

end
