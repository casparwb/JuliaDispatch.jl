@everywhere using SharedArrays

""" Read snapshots with  """
function read_snapshots_p(N;run="", data="", verbose=0)

    nsnaps = get_n_snapshots(run=run, data=data)
    # snaps = SharedArray{Dict}(nsnaps)

    @distributed for i = 0:nsnaps-1
        snap = snapshot(i, run=run, data=data)
    end

    # return snaps
end

function read_snapshots(;run="", data="", verbose=0)

    nsnaps = get_n_snapshots(run=run, data=data)
    snaps = Array{Dict, 1}(undef, nsnaps)

    for i = 0:nsnaps-1
        snaps[i+1] = snapshot(i, run=run, data=data)
    end

    return snaps
end
