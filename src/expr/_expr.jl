macro eval_expression(patch, expr)
    return :(expr)
end

function evaluate_expression(patch, expr; verbose = 0)

    verbose > 0 && println("Parsing expression $expr")

    expr = replace(expr, "sqrt" => "√")
    ops = "+-*/^"
    methods = ["cos", "sin", "tan", "exp", "√"]
    where_methods = [findall(method, expr) for method in methods]
    if !isempty(where_methods)
        idxs = Int[]
        for ids in where_methods
            isempty(ids) && continue
            for id in ids
                append!(idxs, id)
            end
        end
    end

    all_vars = patch["idx"]["dict"] |> keys |> collect
    for var in all_vars
        m = match(Regex("[$var]{$(length(var))}"), expr)
        if !isnothing(m) && !(m.offset in idxs)  
            expr = replace(expr, "$var" => """localpatch["var"]("$var")""")
        end
    end

    for (method, op) in zip(methods, ops)
        if occursin(op, expr)
            expr = replace(expr, op => " .$op ")
        end
        if occursin(method, expr)
            expr = replace(expr, method => " $method.")
        end
    end
    expr = replace(expr, "√" => "sqrt")
    verbose > 1 && println("evaluating expression $expr")
    parsed = Meta.parse(expr)
    return @eval begin
        let localpatch = $patch
            $parsed
        end
    end
    # try
    #     return eval(parsed)
    # catch e
    #     throw(e)
    # end
    # #@eval_expression(patch, expr)
    # test = @eval function(localpatch)
    #     # localpatch = patch
    #     $parsed
    # end
    # test(patch)
end
