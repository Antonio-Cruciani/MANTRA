# APXTBC
Repository associated to the paper "On approximating the temporal betweennes centrality through sampling"

# Network file format
The input file must be a temporal edge list, one line for each temporal edge (a temporal edge is a triple `(u,v,t)`, where `u` and `v` are two vertices and `t` is one time in which the edge from `u` to `v` is available). All the graph are considered directed (if a graph is undirected it must have both temporal edges `(u,v,t)` and `(v,u,t)`). Nodes can be identified by any string (not necessarily a number), while time instants must be integer valued. The three elements of a temporal edge can be separated by any string, which can be specified while using the `load_temporal_graph` function (see below). In the following, we assume that this string is just a white-space.

The file can include duplicate lines and self-loops, however the tool will ignore them. In case your file contains duplicate lines or self-loops lines and you want to manually remove them, you can first execute the following `bash` command.

```
awk '($1!=$2)' <in_file_name> | sort | uniq > <out_file_name>.txt
```

You can then use the newly generated file as input to the tool.

# Note on multithreading
In order to properly run the multi-threading implementations of our approaches you need to set the number of threads that you desire to use with: 
```
 julia --threads 16
```
where the number indicates the number of thread to use, for more information we refer to the official documentation: [Starting Julia with multiple threads](https://docs.julialang.org/en/v1/manual/multi-threading/)

# How to use the temporal graph analysis tool
We assume that the `julia` (or `julia --threads 16`) REPL has been initialized within the directory `APXTBC` and the following command has been already executed (after having installed all he required packages)
```
include("src/APXTBC.jl")
```

# Reading the temporal network
Assuming that the file `workplace.txt` is contained in the directory `graphs/test/` included in the directory `APXTBC`, the corresponding temporal network can be loaded by executing the following command (note that, in this case, the separator between the three elements of a temporal edge is the space character).
```
tg = load_temporal_graph("graphs/workplace.txt", " ")
```
We can print some basic statistics about the temporal network by executing the following command
```
print_stats(tg, graph_name="Workplace network")
```
The result in this case should be the following output.
```
====================================================
Temporal network: Workplace network
====================================================
Number of nodes 92
Number temporal of edges 19654
Number of unique time stamps 7104
====================================================
```
# Running the multi-thread implementation
To run the multi-thread implementation of our algorithms you need to add the `threaded_` prefix to the algorithm name you wish to run, for example:
```
threaded_temporal_shortest_betweenness(tg,50,false)
```
runs the multi-thread implementation of the exact shortest-temporal betweenness algorithm `shortest_betweenness(tg,50,false)`
Here we show ho to run the multi-threading versions, in order to run the single-threaded remove the prefix `threaded_` from the functions' name.

# Computing the exact temporal betweenness scores

The values of the (*)-temporal betweenness of the graph `tg` can be computed by executing the following command.

#### shortest-temporal betweenness

```
stb,t_stb = threaded_temporal_shortest_betweenness(tg,50,false)
```

#### shortest-foremost-temporal betweenness
```
sftb,t_sftb = threaded_temporal_shortest_foremost_betweenness(tg,50,false)
```
#### prefix-foremost-temporal betweenness

```
pftb,t_pftb = threaded_prefix_foremost_betweenness(tg,50,false)
```

The second parameter specifies after how many processed nodes a message has to be printed on the console (in order to verify the status of the computation). If this parameter is `0`, then no output is produced. The third parameter specifies whether the big integer implementation has to be used (in case there too many shortest paths). The execution of the above command should require less than a minute. The values returned are the array of the (*)-temporal betweenness values and the execution time.

The values of the (*)-temporal betweenness and the execution time can be saved in the scores and the times directory, respectively, as follows.

Shortest-temporal betweenness
```
save_results("workplace", "sh", stb, t_stb)
```

Shortest-foremost-temporal betweenness
```
save_results("workplace", "sfm", sftb, t_sftb)
```

Prefix-foremost-temporal betweenness
```
save_results("workplace", "pfm", pftb, t_pftb);
```

# Computing the approximations 
Similarly to the computation of the (*)-temporal betweenness, we can compute its approximations using the sampling-based algorithms described in the paper: Random Temporal Betweenness (rtb), ONBRA (ob), Temporal Riondato and Kornaropoloulos (trk). Depending on your needs, you can execute the fixed sample size, or progressive sampling version of these algorithms.

## Fixed sample size
Every `fixed sample size` approximation algorithm takes as input: (1) temporal graph; (2) sample size; (3) verbose step; and, if not prefx-moremost also (4) boolean value for bigInt data-structure.

#### Random Temporal Betweenness (rtb)
##### shortest-temporal betweenness

```
apx_tbc,t_est = threaded_rtb(tg,20,10,false)
```

```
save_results_sampling("workplace", "rtb_sh", apx_tbc, 20,t_est)
```
##### shortest-foremost-temporal betweenness

```
apx_bc,t_est = threaded_rtb_shortest_foremost(tg,20,10,false)
```

```
save_results_sampling("workplace", "rtb_sfm", apx_tbc, 20,t_est)
```

##### prefix-foremost-temporal betweenness

```
apx_bc,t_est = threaded_rtb_prefix_foremost(tg,10,20)
```

```
save_results_sampling("workplace", "rtb_pfm", apx_tbc,20, t_est)
```

#### ONBRA (ob)
##### shortest-temporal betweenness

```
apx_tbc,t_est = threaded_onbra(tg,20,10,false)
```

```
save_results_sampling("workplace", "onbra_sh", apx_tbc,20, t_est)
```
##### shortest-foremost-temporal betweenness

```
apx_bc,t_est = threaded_onbra_shortest_foremost(tg,20,10,false)
```

```
save_results_sampling("workplace", "onbra_sfm", apx_tbc, 20,t_est)
```

##### prefix-foremost-temporal betweenness
```
apx_bc,t_est = threaded_onbra_prefix_foremost(tg,20,10)
```

```
save_results_sampling("workplace", "onbra_pfm", apx_tbc,20, t_est)
```

#### Temporal Riondato and Kornaropoloulos (trk)
##### shortest-temporal betweenness
```
apx_tbc,t_est = threaded_trk(tg,20,10,false)
```

```
save_results_sampling("workplace", "trk_sh", apx_tbc, 20,t_est)
```
##### shortest-foremost-temporal betweenness
```
apx_bc,t_est = threaded_trk_shortest_foremost(tg,20,false)
```

```
save_results_sampling("workplace", "trk_sfm", apx_tbc,20, t_est)
```

##### prefix-foremost-temporal betweenness
```
apx_bc,t_est = threaded_trk_prefix_foremost(tg,20)
```

```
save_results_sampling("workplace", "trk_pfm", apx_tbc,20, t_est)
```

## Progressive sampling algorithms
Here we illustrate how to run the `progressive sampling` approximation algorithms.

#### Progressive-Random Temporal Betweenness (p-rtb)

##### shortest-temporal betweenness
```
apx_bc,n_s,t_est  = progressive_rtb(tg,2,10,false)
```

##### shortest-foremost-temporal betweenness
```
apx_bc,n_s,t_est = progressive_shortest_foremost_rtb(tg,2,100,false)
```
##### prefix-foremost-temporal betweenness

```
apx_bc,n_s,t_est  = progressive_prefix_foremost_rtb(tg,2,10)
```
The second parameter is the threshold value that the algorithm uses to establish when to stop, the third how many processed nodes a message has to be printed on the console (in order to verify the status of the computation). If this parameter is `0`, then no output is disabled, the last one (if present) indicates if the algorithm has to use the BigINT data-structure.
To save the output we can call the function

```
save_results_progressive_sampling("workplace","p_rtb",apx_bc,n_s,t_est,2)
```
Where the third element is the tolerance parameter, the fourth is the estimation time, the fifth is the overall number of samples used by the algorithm to converge, and the last one is the geometric value for the progressive sampling (this algorithm does not use this technique, so you can put an arbitrary value).

#### Progressive-ONBRA (p-ob)

##### shortest-temporal betweenness
```
apx_bc,n_s,xi,t_est  = threaded_progressive_onbra(tg,350,0.1,0.1,1.5,50,false)
```

##### shortest-foremost-temporal betweenness
```
apx_bc,n_s,xi,t_est  = threaded_progressive_onbra_shortest_foremost(tg,350,0.1,0.1,1.5,50,false)
```
##### prefix-foremost-temporal betweenness

```
apx_bc,n_s,xi,t_est  = threaded_progressive_onbra_prefix_foremost(tg,350,0.1,0.1,1.5,50)
```

The second element is the starting sample size, the third is the accuracy $\varepsilon\in (0,1)$, the fourth is the acceptable failure probability $\delta\in (0,1)$, the fifth is the geometric schedule parameter, the sixth how many processed nodes a message has to be printed on the console (in order to verify the status of the computation). If this parameter is `0`, then no output is disabled, the last one (if present) indicates if the algorithm has to use the BigINT data-structure.
To save the results we can use
```
save_results_progressive_sampling("workplace","p_onbra",apx_bc,n_s,t_est,350,xi)
```

#### Progressive-ONBRA (p-trk)

##### shortest-temporal betweenness
```
apx_bc,lb,ub,n_s,ss,t_est = threaded_progressive_trk(tg,0.1,0.1,50,false,-1,100,10,true)
```

##### shortest-foremost-temporal betweenness
```
apx_bc,lb,ub,n_s,ss,t_est = threaded_progressive_trk_shortest_foremost(tg,0.1,0.1,50,false,-1,100,10,true)
```

##### shortest-foremost-temporal betweenness

```
apx_bc,lb,ub,n_s,ss,t_est = threaded_progressive_trk_prefix_foremost(tg,0.01,0.1,100,-1,100,10,true)
```


The second element is the accuracy $\varepsilon\in (0,1)$, the third is the acceptable failure probability $\delta\in (0,1)$, the fourth how many processed nodes a message has to be printed on the console (in order to verify the status of the computation). If this parameter is `0`, then no output is disabled, the fifth one (if present) indicates if the algorithm has to use the BigINT data-structure.

From the sixth to the ninth:

-	Diameter of the graph, default value -1
- 	Start factor, default 100 (do not change)
- Sample step, default 10, is the number of consecutive samples at each iteration (every 10 steps, the algorithm checks the convergence)
- Hoeffeding's inequality, default false:
	- if false, estimates the shortest-temporal diameter then uses the bound on the VC Dimension to establish the sample size
	- if ture, uses the Hoeffeding's inequality to establish the sample size

To save the results we can run
```
save_results_progressive_sampling("workplace","p_rtk",apx_bc,n_s,t_est,ss,-1.0)
```
