
"""
    evaluate_expression(patch::Dict, expr::string; verbose::Int=0)

Attempt to parse and evaluate a given expression `expr`. To include a value from a variable,
interpolate it into the expression string using the interpolation (dollar-sign) construct. 

# Example
```
julia> gamma = snap["gamma"] 
julia> snap["patches"][1]["var"]("gamma*π/4*sqrt(bx^2 + by^2 + bz^2"))
```
"""
function evaluate_expression(patch, expr; all=false, verbose = 0)

    expr = replace(expr, "sqrt" => "√")

    if occursin("log10", expr)
        expr = replace(expr, "log10" => "lg")
    end

    if occursin("log2", expr)
        expr = replace(expr, "log2" => "lgg")
    end

    if occursin("log", expr)
        expr = replace(expr, "log" => "ln")
    end
    ops = "+-*/^"

    methods = ["cos", "sin", "tan", "√", "abs", "atan", "ln", "lgg", "lg"]
    bools = ["true", "false"]

    """ 
    if there are any mathematical expressions inside the expression, 
    locate their indices 
    """
    idxs = Int[]
    if any(occursin.(methods, expr))
        where_methods = [findall(method, expr) for method in methods]
        if !isempty(where_methods)
            for ids in where_methods
                isempty(ids) && continue
                for id in ids
                    append!(idxs, id)
                end
            end
        end
    end


    # find the quantities in the expression that are being used 
    all_vars = patch["idx"]["dict"] |> keys
    for var in all_vars
        m = match(Regex("[$var]{$(length(var))}"), expr)
        if !isnothing(m) && !(m.offset in idxs)  
            expr = replace(expr, "$var" => """patch["var"]("$var")""")
        end
    end

    for var in all_vars
        expr = replace(expr, """patch["var"]("$var")""" => """patch["var"]("$var", all=$all)""")
    end


    where_bools = [findall(bool, expr) for bool in bools]
    if !isempty(where_bools)
        for ids in where_bools
            isempty(ids) && continue
            for id in ids
                append!(idxs, id)
            end
        end
    end

    # vectorize the operations and mathemetical expressions
    for op in ops
        if occursin(op, expr)
            expr = replace(expr, op => " .$op ")
        end
    end


        # vectorize the operations and mathemetical expressions
    for method in methods
        if occursin(method, expr)
            expr = replace(expr, method => " $method.")
        end
    end

    expr = replace(expr, "ln" => "log")
    expr = replace(expr, "lg" => "log10")
    expr = replace(expr, "lgg" => "log2")

    expr = replace(expr, "√" => "sqrt") # sqrt works and √ doesn't for some reason

    # parse and evaluate
    verbose >= 1 && println("evaluating expression $expr")
    parsed = Meta.parse(expr)
    try
        # eval() looks for variables in global scope, so we have to force the
        # patch argument to be in its local scope by using the `let` construct
        return @eval begin
            let patch = $patch 
                $parsed
            end
        end
    catch e
        println("Unable to parse expression")
        throw(e)
    end

end
