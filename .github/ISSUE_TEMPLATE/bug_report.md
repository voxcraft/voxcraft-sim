---
name: Bug report
about: Create a bug report to help us improve voxcraft-sim
title: "[Bug]"
labels: "? - Needs Triage, bug"
assignees: ''

---

**Summary of the Bug**
A clear and concise description of what the bug is.

# Steps to Reproduce

Thank you for using voxcraft-sim. 
We encourage you to open issues. 
We rely on this software for research, and peer review is the best way to find and eliminate errors.

However, at the time of writing, our team consists of just three people. 
We are students, scientists, and parents. We are not tech support. 

So -- please try debugging before you open an issue. 
If you fix a specific bug or make a general improvement in your fork, please submit a pull request. 

If your issue persists, please think about the best way to communicate and demonstrate it. 

If we cannot replicate your demo, we cannot verify that your issue exists.

Follow this guide http://matthewrocklin.com/blog/work/2018/02/28/minimal-bug-reports to craft a minimal bug report. This helps us reproduce the issue you're having and resolve the issue more quickly.


**Please complete the following information:**
    - Your base.vxa file
    - Your robot.vxd file
    - Your cmake log

    - What GPU are you using?
        - GPU Name:
        - Please use  ```nvidia-smi -L``` to determine this.
    - What branch are you on?
        - Branch Name:
        - Are there any uncommitted files in your git workspace?
        - Yes / No: 
        - Please use ```git status``` and to confirm.
    - If you have uncommitted files we will not be able to reproduce your bug.
        - What is your git commit hash?
        - Commit hash:
        - Please use ```git log -n 1``` to display it.
 

**Additional details**
Add any other details or context about the problem here.