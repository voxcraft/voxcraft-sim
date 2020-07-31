import subprocess as sub


sub.call("cp ../../build/voxcraft-sim .", shell=True)

sub.call("cp ../../build/vx3_node_worker .", shell=True)

sub.call("rm replicate.hist", shell=True)

print("saving history as replicate.hist")

sub.call("./voxcraft-sim -i data > replicate.hist", shell=True)
