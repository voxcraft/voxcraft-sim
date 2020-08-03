# Reporting voxcraft-sim issues on github

Thank you for using voxcraft-sim. 
We encourage you to open issues. 
We rely on this software for research, and peer review is the best way to find and eliminate errors.

However, at the time of writing, our team consists of just three people. 
We are students, scientists, and parents. We are not tech support. 

So: please try debugging before you open an issue. 
If you fix a specific bug or make a general improvement in your fork, please submit a pull request. 

If your issue persists, please think about the best way to communicate and demonstrate it. 

If we cannot replicate your demo, we cannot verify that your issue exists.

### Create a demo with as few voxels as possible. The demo should include:

**1.** &nbsp;&nbsp; A base.vxa file.

**2.** &nbsp;&nbsp; A robot.vxd file. 

**3.** &nbsp;&nbsp; Your cmake log.

**4.** &nbsp;&nbsp; The specific GPU you are using.
Please use ```nvidia-smi -L``` to determine this. 
If you are using google colab, this might change so please get the current gpu. 

**5.** &nbsp;&nbsp; Does your issue also occur on the master branch? 
Please confirm all source files are being tracked by git and are committed and pushed (```git status``` to confirm). 
Then please post the git commit that you are using (```git log -n 1``` to display)

