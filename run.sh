pyhton3 allen_cahn.py dumbel --eps 0.2
pyhton3 allen_cahn.py dumbel --eps 0.01

python3 plot.py data/dumbel/eps_0.01000000.vtu
python3 plot.py data/dumbel/eps_0.01000015.vtu
python3 plot.py data/dumbel/eps_0.01000101.vtu

python3 plot.py data/dumbel/eps_0.2000000.vtu
python3 plot.py data/dumbel/eps_0.2000015.vtu
python3 plot.py data/dumbel/eps_0.2000101.vtu


python3 allen_cahn.py bump --eps 0.2
python3 allen_cahn.py bump --eps 0.01

python3 plot.py data/bump/eps_0.01000000.vtu
python3 plot.py data/bump/eps_0.01000015.vtu
python3 plot.py data/bump/eps_0.01000101.vtu

python3 plot.py data/bump/eps_0.2000000.vtu
python3 plot.py data/bump/eps_0.2000015.vtu
python3 plot.py data/bump/eps_0.2000101.vtu

