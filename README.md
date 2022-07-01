# KSW2 - Genomic Alignment Acceleration


Team number: xohw22-088 <br />
Project name: KSW2 - Genomic Alignment Acceleration <br />

### Structure of the repository
src: contains source files <br />
bitstream: contains .xclbin files <br />

### Building the bitstream

You are required to have Vitis 2020.2 and Xilinx Runtime installed together with the necessary files to build for the Xilinx Alveo U280.

To build the project simply type:
```
$> make all TARGET=hw
```

### Running our implementation on FPGA

Firstly, run the following command to build the host: 
```
$> make host

```

Then, run the following commands to run everything on FPGA: <br />
```
$> ./host bitstream/ksw2.xclbin

```
This will run 100k alignments of the FPGA with randomly generated sequences 256 chars long.
