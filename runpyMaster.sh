mpirun -np 17 -hostfile mpihosts -mca plm_rsh_no_tree_spawn 1 --map-by node python MatrixMultErrorDetection.py
mpirun -np 33 -hostfile mpihosts -mca plm_rsh_no_tree_spawn 1 --map-by node python MatrixMultErrorDetection.py
mpirun -np 65 -hostfile mpihosts -mca plm_rsh_no_tree_spawn 1 --map-by node python MatrixMultErrorDetection.py
mpirun -np 129 -hostfile mpihosts -mca plm_rsh_no_tree_spawn 1 --map-by node python MatrixMultErrorDetection.py
