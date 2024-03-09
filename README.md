# MANTRA
Repository associated to the paper "MANTRA: Temporal Betweenness Centrality Approximation through Sampling".

# Network file format
The input file must be a temporal edge list, one line for each temporal edge (a temporal edge is a triple `(u,v,t)`, where `u` and `v` are two vertices and `t` is one time in which the edge from `u` to `v` is available). All the graph are considered directed (if a graph is undirected it must have both temporal edges `(u,v,t)` and `(v,u,t)`). Nodes can be identified by any string (not necessarily a number), while time instants must be integer valued. The three elements of a temporal edge can be separated by any string, which can be specified while using the `load_temporal_graph` function (see below). In the following, we assume that this string is just a white-space.

The file can include duplicate lines and self-loops, however the tool will ignore them. In case your file contains duplicate lines or self-loops lines and you want to manually remove them, you can first execute the following `bash` command.

```
awk '($1!=$2)' <in_file_name> | sort | uniq > <out_file_name>.txt
```

You can then use the newly generated file as input to the tool.

# Note on multi-threading
In order to properly run the multi-threading implementations of our approaches you need to set the number of threads that you desire to use with: 
```
 julia --threads 16
```
where the number indicates the number of thread to use, for more information we refer to the official documentation: [Starting Julia with multiple threads](https://docs.julialang.org/en/v1/manual/multi-threading/)

# How to reproduce the experiments in the paper
To reproduce all the experiment you need to: 

(i) Compute the exact temporal betweenness values running the command `julia --threads <nthreads> compute_exact_temporal_bc.jl`

(ii) Compute MANTRA's approximation, running `julia --threads <nthreads> compute_mantra.jl`

(iii) Compute ONBRAS's approximation, running `julia --threads <nthreads> compute_onbra.jl`

where `<nthreads>` is the number of assigned threads.
All the results will be automatically saved in the `scores` and `times` folders.

Finally, to summarize the results run ```julia analyze_results.jl```, this command will create a new folder called `analysis` with all the results showed in the paper.

# How to use MANTRA to analyze temporal graphs

Here we describe how to run MANTRA on a general input temporal graph.

## MANTRA parameters
The algorithm has the following input parameters:

(-) Temporal Graph (`temporal_graph`)

(-) The absolute approximation precision epsilon (`Float64`)

(-) The confidence value delta (`Float64`)

(-) Whether the BigInt data type must be used (`Bool`)

(-) The estimator that must be used (`String` - options: 1- `ob`:ONBRA estimator, 2-`trk` (default) TRK estimator, 3-`rtb` random temporal betweenness estimator)

(-) The temporal path optimality that must be considered while computing the temporal betweenness (`String` - options: 1- `sh` (default) Shortest, 2- `sfm` Shortest-foremost, 3- `pfm` Prefix-foremost)

(-) Whether to use the VC-Upperbound or the Variance Upperbound (`Bool` , `False` default, i.e., Variance Upperbound is used as default)

(-) The geometric sampler value (`Float64`, default `1.2`)

(-) Diameter of the analyzed input graph (`Int64`, `-1` as defalut value. If the value is set to `-1` the algorithm approximates the diameter using the fast-sampling method described in our paper)

## Example on how to run MANTRA

Open Julia Repl in the main folder (where all the .jl scripts are) using the command 

```julia --threads <number of threads>```

Next, include MANTRA

```include("src/MANTRA.jl")```


Load the temporal graph that must be analyzed

```tg = load_temporal_graph("graphs/00_workplace.txt"," ")```

notice: the temporal graph loader has two input parameters: (i) the temporal graph file and (ii) the file separator

Run MANTRA,

``` resut = progressive_mantra(tg,0.07,0.1,false,"ob","pfm")```

the above command will compute a (0.07,0.1)-Approximation using the ONBRA (ob) estimator for the Prefix Foremost Temporal Betweenness

Save the output to file

```save_results_progressive_sampling("workplace","mantra_ob_pfm",result[1],result[2],result[4],1,0.07)```
