import subprocess as sub

for r in range(100, 500):
    sub.call("sbatch submit.sh {0} 4 3 20 50".format(r), shell=True)

exit()
