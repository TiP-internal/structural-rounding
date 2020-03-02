
# Structural Rounding
~~See example.py if you'd like to write code that makes use of our structural rounding algorithms~~

See ```example.yaml``` if you'd like to create a config file for use with structural rounding.


## Compiling
- **Python** Simply run ```make python``` to compile.
- **C++** Simply run ```make cpp``` to compile.


## Generating Synthetic Graphs
Compile the generator by running ```make generator```.
Then, running ```./create_graphs.sh small``` will create 5 graphs per parameter setting (3600 graphs in total) each with 4 million edges in expectation.
To create smaller graphs, use ```./create_graphs.sh test``` which creates graphs with 100 thousand edges in expectation.
To create larger graphs, use ```./create_graphs.sh medium``` (40 million edges) or ```./create_graphs.sh large``` (200 million edges).

Note that the largest graphs can exceed 600MB in size even with a space efficient encoding.

You can also create different sizes of graphs using ```./create_graphs.sh custom <edges> <directory>```.


## Running Experiments
- **Python** Once compiled, run ```python3 main.py ...optional arguments...```

    - Arguments set in the command line have precedent over those set in the config file.
    - Putting -s alone as an argument will look for *config.yaml* by default.
    - Putting -r alone as an argument will write to *results.csv* by default.


- **C++** Once compiled, run ```./main ...optional arguments...```

    - Currently, no optional help command.
    - Currently, no config file parsing.

    Further information on our optional arguments:

    ```-h, --help``` show this message and exit

    ```-p, --problem``` **required** the problem to solve: [vertex_cover]

    ```-c, --class``` the graph class to edit to: [bipartite]

    ```-e, --edit``` the editing algorithm [remove_octset]

    ```-l, --lift``` the lifting algorithm [greedy, nieve]

    ```-a, --approx``` accessory approximation algorithms for the problem [dfs, heuristic, std]

    ```-g, --graph``` **required** the path to the graph file/dir

    ```-s, --spec``` the (optional) config file (.yaml) where these arguments can also be assigned

    ```-r, --results``` **required** the results file (.csv) where results are printed to


~~Once you have created synthetic graphs, you can reproduce our experimental results by running ```make small_data```.
Use ```make medium_data``` or ```make large_data``` if appropriate.
You can run our experiments on different sizes of graphs using ```python main.py <directory>```.
Note that the ```make``` commands additionally disable Python's random hashing feature so that results are consistent between runs.~~

~~The variance experiments can also be reproduced using ```make test_data```.
Do the large amount of repetition, it is not recommended to use graphs with more than 100 thousand edges in the variance tests.
To run the variance experiments on graphs of a different size, use ```python variance.py <directory>```.~~
