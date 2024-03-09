include("src/APXTBC.jl")


datasets = [
    "23_wiki_talk.txt"


]

path_opt = ["pfm"]
for p in path_opt
    get_max_temporal_bc(p,datasets)
end
