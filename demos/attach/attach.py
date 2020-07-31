import subprocess as sub


sub.call("cp ../../build/voxcraft-sim .", shell=True)

sub.call("cp ../../build/vx3_node_worker .", shell=True)

sub.call("rm attach*.hist", shell=True)

sub.call("cp tests/test1.vxa data/base.vxa", shell=True)

print("Running Test 1 as attach_test1.hist")

sub.call("./voxcraft-sim -i data > attach_test1.hist", shell=True)

sub.call("cp tests/test2.vxa data/base.vxa", shell=True)

print("Running Test 2 as attach_test2.hist")

sub.call("./voxcraft-sim -i data > attach_test2.hist", shell=True)

sub.call("cp tests/test3.vxa data/base.vxa", shell=True)

print("Running Test 3 as attach_test3.hist")

sub.call("./voxcraft-sim -i data > attach_test3.hist", shell=True)
