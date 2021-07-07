using Unitful, Latexify, UnitfulLatexify, UnitfulRecipes

function get_unit(iv)

    if startswith(iv, "u") || startswith(iv, "v")
        unit = u"m/s"
end