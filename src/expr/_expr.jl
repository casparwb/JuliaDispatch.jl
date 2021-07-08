"""
18.06.2020: not finished. to-do: finish evaluate() and expression_parser()
what is locals()?
"""


mutable struct Expressions
    verbose::Int
    expressions::Array{String}
    ops::String
    known::Array{String}
    pvars::Array{String}
    patch::Dict
end

function init(patch = nothing, verbose = 0)
    """
    Initialize expression
    """
    expressions = []
    ops = " +-*/()"
    known = []
    pvars = []

    if patch != nothing
        for (k, v) in patch["keys"]
            k != "numbers" ? append!(pvars, v) : nothing
        end
        verbose ? println("patch_variables: $pvars") : nothing
    end
    verbose > 1 println("EXPRESSIONS: ") :nothing ## WHAT?

    Expression(verbose, expressions, ops, known, pvars, patch)

end

function split(expressions::Expressions, rhs)
    """ Split a RHS into words, based on operators in expressions.ops """

    words = Array{String}[]

    add(w) = (w != "" && !occursin(w, words)) ? push!(words, w) : nothing
    word = ""
    for c in rhs
        if c in expressions.ops
            add(word)
            word = ""
        else
            word += c
        end
    end
    add(word)
    return words
end

function register(obj::Expressions, expression)
    """ Register expressions, checking each RHS for words that are either
        known patch variables, previsouly defined, or scalars
    """

    verbose > 1 ? println("$expression") : nothing

    lhs, rhs = split(expression, "=")
    if word in obj.known || word in obj.pvars
        continue
    else
        if word in globals()
            println("WARINING: $word in globals()")
        end
        try
            eval(word)
            push!(obj.known, word)
        catch
            println("ERROR: could not evaluate $word")
            return false
        end
    end

    push!(obj.known, lhs)
    return true
end

function evaluate(obj::Expressions)
    """ After all expressions registered, evaluate """
    verbose ? println("EVALUATE: ") : nothing

    ll = locals()


end
