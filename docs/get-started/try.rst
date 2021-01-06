Try
===

Google Colab provides a free online GPU environment.

.. code-block:: python

    # First, go to menu->Runtime->Change runtime type, select GPU.
    # Then, run the script:
    !git clone https://github.com/liusida/gpuVoxels.git; cd gpuVoxels/; git checkout dev; git pull; 
    !mkdir build; cd build; rm * -rf; cmake ../gpuVoxels/Voxelyze3; make -j 10; 
    !cd build; ./Voxelyze3 -i ../gpuVoxels/zoo/basic/ -o output.xml -f > ../a.history
    print("Done! Start downloading the generated a.history. Please use your local VoxCAD to open the file.")
    from google.colab import files
    files.download('a.history')

Readers can also check this readonly example:

https://colab.research.google.com/drive/1yiqw7Uq3W3CgYCinXq4t808M2l7uuLv1?usp=sharing