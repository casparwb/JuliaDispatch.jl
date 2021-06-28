


function nml_to_dict(filepath)

    namelist = Dict{String, Any}()

    open(filepath, "r") do file
        lines = readlines(file)

        idx = 1
        while true
            if idx == length(lines) break end
            line = lines[idx]
            if startswith(line, "&") && lines[idx+1] != "/"
                println(line)
                words = split(line)
                key = replace(
                      replace(words[1][2:end], "_params" => "")
                                             , "_nml" => "")


                namelist[key] = Dict{String, Any}()
                if length(words) > 1
                    for word in words
                        if word == "/" break end
                        params = strip.(split(words, "="))
                        namelist[key][params[1]] = params[2]
                    end
                    idx += 1
                else
                    i = 1
                    while true
                        newline = lines[idx+i]
                        if newline != "/"
                            params = strip.(split(words, "="))

                            if params[2] == "t" || params[2] == "f"
                                val = ("true", "false")[params[2] .== ["t", "f"]][1]
                                namelist[key][params[1]] = val
                            else
                                val = params[2]
                                try
                                    val = parse(Int, params[2])
                                catch
                                    val = parse(Float32, params[2])
                                finally
                                    namelist[key][params[1]] = val
                                end
                            end
                            i += 1
                        else
                            idx += i
                            break
                        end
                    end
                end
            else
                idx += 1
            end
        end
    end

    return namelist


end
