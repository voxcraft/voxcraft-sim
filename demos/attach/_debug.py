import subprocess as sub


sub.call("rm a_debug.hist", shell=True)

print("saving history as a_debug.hist")

sub.call("./voxcraft-sim -i _debug > a_debug.hist", shell=True)
