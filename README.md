# E2Pilot.jl

This software code is to support our manuscript "Maximizing Heavy-Duty E-Truck Decarbonization by Carbon-Optimized Timely Transportation".

## Installation

### Step 0 Download the code
Due to the storage restriction of Github, we may not upload the binary data (e.g., network data, carbon intensity data, etc.) via Git. 
Therefore, we recommend the user to download the bundle of the code and data (named `cfo-bundle.zip`)  from the [release page](https://github.com/researcher542/cfo/releases/) of this repository.


### Step 1 Install Julia 
This software code is written in the [Julia](https://julialang.org/) Language.
The users are kindly directed to follow the [installation guide](https://julialang.org/downloads/) to install the Julia.

### Step 2 Install the dependency
Once Julia is installed, we can install the dependency via the following command:
```[Julia]
cd path_to_project_root              # Go to the root folder of the project. 
cd E2PilotCFO
julia --project=. -t auto            # Enter the Julia REPL 
import Pkg
Pkg.instantiate()
```
This step might take 5-10 minutes depending on the hardware performance and the network condition.

### Troubleshooting

#### Newer version of Julia
Note that is mainly tested on Julia v1.11. However, it is generally compatible with newer Julia version. We can use the following command to update the dependency to the newer versions
```
Pkg.update()
```


### Tested Platforms
This software has been tested on the following platforms:
- Ubuntu 2404 virtual machine with 8 cores and 16GB of RAM
- Windows 11 with Intel i5-13400F CPU and 64GB of RAM
- Macbook Pro 2024 with a M4 Pro chip and 24 GB of RAM.

For all the platforms, the Julia v1.11 is used.

## Demo 

### If you are in the Julia REPL. That is, you can see a green `julia>` in the terminal.
We can simple run the following command to execute the demo.
```
include("paper/get_one_result.jl")  # Run the demo
```

### If you are not in the Julia REPL
You can run the command `julia --project=. -t auto` to enter the Julia REPL. Or alternatively, we can run the following command from the bash terminal to execute the demo.
```[Julia]
julia --project=. -t auto paper/get_one_result.jl
```
Here `--project=.` specifies the project file. `-t auto` specifies the multi-threading mode. In the `get_one_result.jl` script, we give an example of how to use the high-level API to get the result of a single instance.



### Expected Results
The run time to finish the above simulation is about 3 minutes on an M4 Pro chip with 8 threads.

The computed carbon-optimized solution is stored in a variable called `ds_res`. The fastest path solution is stored in the variable `fast_ds_res`.

At the end of the program, it will print a summary of the fastest path and the carbon optimized solution.
```
****************************************
Results for the carbon-optmized solution.
Carbon footprint: 1949.1335454807295 kg
Time: 25.25003914865546 hours
Energy: 2441.781585081727 kWh
****************************************
Results for the fastest path.
Carbon footprint: 4087.667898828323 kg
Time: 21.058834627041275 hours
Energy: 2629.960641730901 kWh
```



### How to get all the results

We can run the following command to get the required results for the first pair of origin and destination, including instances with different objectives, different deadlines, and different alphas, etc. 

```
julia --project=. -t auto paper/get_one_src_des.jl --idx 1 --region oldmap
```
Here `1` stands for the first pair. We can change it to values between 1 and 400. 
The `region` can take the values of `oldmap` and `eu' for US network and EU network, respectively.

**Note**: The above command requires us to exit the Julia REPL (with function `exit()` inside Julia REPL) first.

**Note**: The above script might take several hours and consume about 16GB RAM to complete for one pair of origin and destination.


## Useful Scripts
All the experiments scriptsare placed in the E2PilotCFO/paper/ folder.
- `get_one_result.jl`: The script to get the result of a single instance. 
- `get_one_src_des.jl`: The script to get the results of all the instances for one origin and destination.
- `method_cmp.jl`: The script to run the experiments to compare the performance of our approach with other two alternatives: MIP formulation with Gurobi solver and the time-expanded and battery-expanded graph approach.
- `get_cs_net.jl`: The scripts to get the charging station network used to store some pre-computed energy-efficient paths between charging stations. Used to make the algorithm faster
- `parse_cambium_data_sep.jl`: The script to parse the data of the Cambium dataset which store the carbon intensity data for each state and have projection to the future.


   

