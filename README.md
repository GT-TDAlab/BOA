# BOA

BOA (stands for "Bucket-Order-Assemble") is a parallel de novo genome assembly framework that utilizes
hypergraph and graph partitioning. This library is designed to improve
assembly quality as well as expose a high degree of parallelism for standalone assemblers.
Our approach models the assembly problem as one of either graph or hypergraph partitioning.
We support two implementations:
HyperBOA (which uses hypergraph partitioning) and
GraphBOA (which uses graph partitioning).
Input is a set of DNA short reads.
The user also can specify a standalone (sequential) assembler of choice, to use in the last step of our framework.
The output is a set of BOA-assembled contigs.
(Based on results, we recommend using HyperBOA.)


**Source Code:** [https://github.com/GT-TDAlab/BOA]

## License

BOA is distributed under BSD License. For details, see [`LICENSE`](LICENSE).

## Contributors

- [Priyanka Ghosh](https://www.linkedin.com/in/priyankaghosh/)
- [Xiaojing An](https://xiaojingan.com/)
- [Patrick Keppler](https://www.linkedin.com/in/patrick-keppler-231138140/)
- [Sureyya Emre Kurt](https://www.linkedin.com/in/semrekurt/)
- [Umit V. Catalyurek](http://cc.gatech.edu/~umit)
- [Sriram Krishnamoorthy](https://hpc.pnl.gov/people/sriram/)
- [P. Sadayappan](https://www.cs.utah.edu/~saday/)
- [Aravind Sukumaran Rajam](https://eecs.wsu.edu/~asukumar/)
- [Ananth Kalyanaraman](https://eecs.wsu.edu/~ananth/)

## Contact

For questions or support [open an issue](https://github.com/GT-TDAlab/BOA/issues/new) or contact contributors via <priyanka.ghosh15@gmail.com>, <anxiaojing@gatech.edu>.

## Citation
Citation for the BOA (BibTeX):

```bibtex
    @inproceedings{Xiaojing22BOA,
        author =  {Xiaojing An, Priyanka Ghosh, Patrick Keppler, Sureyya Emre Kurt, \"Umit V. \c{C}ataly\"urek, Sriram Krishnamoorthy, P. Sadayappan, Aravind Sukumaran Rajam, Ananth Kalyanaraman},
        title = {BOA: A Partitioned View of Genome Assembly},
        booktitle = {Recomb-Seq 2022-12th RECOMB Satellite Workshop on Massively Parallel Sequencing},
        year   = {2022},
        KEYWORDS = {genome assembly, graph partitioning, hypergraph partitioning},
    }
```

## How to build

Directly run `make`

    make -j8


### Dependencies

1) MPI library (preferably MPI-3 compatible)
2) C++14 (or greater) compliant compiler
3) Zoltan 3.83 and ParMetis 4.0.3

#### Installing Zoltan, ParMetis:

>You can download ParMetis through:
>```shell
>    wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz
>```
>Then unzip, and follow the readme.
>Note that, please change in `./metis/include/metis.h`: from `#define IDXTYPEWIDTH 32` to `#define IDXTYPEWIDTH 64`
>
>You can download Zoltan through:
>```shell
>    wget https://github.com/sandialabs/Zoltan/archive/v3.83.tar.gz
>```
>See [Zoltan install instruction](https://cs.sandia.gov/Zoltan/ug_html/ug_usage.html)
>
>Please configure to install with `--with-id-type=ullong` option and with the third party library ParMetis `--with-parmetis
>--with-parmetis-incdir=$ParMetis_include_dir --with-parmetis-libdir=$ParMetis_lib_dir`
>
>Then update [`config.txt`](./config.txt):
>```shell
>    PARMETIS_INSTALL=where ParMetis is installed
>    ZOLTAN_INSTALL=where Zoltan is installed
>    ZOLTAN_HOME=where is the source files of Zoltan
>```

Or

>You can call the script below (without any parameters, ParMetis and Zoltan will be installed at `./thirdparty`):
>
>```shell
>    ./scripts/get_partitioner.sh
>```
>The usage of it:
>```shell
>Usage: ./scripts/get_thirdparties.sh 'whether install Zoltan with ParMetis(0/1)' 'where the thirdparties should be download at' 'where config.txt is'
>    'whether Zoltan install with ParMetis', 0 (No) or 1 (Yes), is default as 1
>    'where thirdparties (Zoltan, ParMetis) should be download at' is default as './thirdparty'
>    'where config.txt is at' is default as './config.txt'
>```

## How to run

At this time, we support the input files to be in the FASTA format and utilize paired-end reads.

```shell
mpiexec -np $procs $BINARY -f $READS -m $READ_MATES -b $MAX_BATCH_SIZE -p $NUM_OF_PARTITIONS -t $PARTITIONER_TYPE -d $SEED -o $OUTPUT_DIR_NAME
```

Mandatory input arguments:

- -f : the input reads file in .fasta format
- -m : the input read mates file in .fasta format for paired-end information

Optional arguments:

- -b : number of k-mers to include in a batch for collective communication.
  Default value set to 100M (100,000,000). If memory is short, consider
  reducing to 50M.
- -p : the number pf partitions. In default, it is the number of MPI processes
- -t : the partitioner used: Zoltan (in default) or ParMetis
- -d : seed for Zoltan partitioning. Default value: 0
- -o : the path for output storing the partitioned reads. Default value: `./output`

For example:

```shell
mpiexec -np 4 ./pbucketing -f ./example_reads/N_1.fa -m ./example_reads/N_2.fa -b 100000000 -p 4 -t Zoltan -o ./output
```
Outputs:
```shell
Start Reading the input dataset
Total number of reads: 11
Completed Reading the input dataset
Start Reading the input dataset
Total number of reads: 11
Completed Reading the input dataset
rank:0, reads_data_size: 36, local_nlines: 6, reads_disp: 0
Starting sliding window
rank:1, reads_data_size: 36, local_nlines: 6, reads_disp: 3
rank:2, reads_data_size: 38, local_nlines: 6, reads_disp: 6
rank:3, reads_data_size: 28, local_nlines: 4, reads_disp: 9
Completed sliding window
Starting sliding window
Completed sliding window
Total distinct k-mer bucket entries across all proc's: 111
Number of outliers (buckets discarded): 828
Number of buckets: 939
percentage of buckets discarded: 88.178914
Completed kmer counting
Average time for performing k-mer bucketing across all procs (secs): 0.000900
Maximum time for performing k-mer bucketing across all procs (secs): 0.000906
Begin Zoltan partitioning
K-mer size: 31, L-mer size: 8, Number of Processes: 4, MAX_KMER_COUNT: 100000000, Number of Partitions: 4,Average time for MPI_FIle_Open (secs): 0.006305
Average time for MPI parallel file reading (secs): 0.005442
Average time for MPI_FIle_Close (secs): 0.001294
Average time for reading and storing the Reads in each proc's memory (secs): 0.013092
Total number of reads: 11
Average time for MPI_FIle_Open (secs): 0.005800
Average time for MPI parallel file reading (secs): 0.002349
Average time for MPI_FIle_Close (secs): 0.001211
Average time for reading and storing the Reads in each proc's memory (secs): 0.009600
Total number of reads: 11
rank: 0, reads_size: 36
transfer kmers last round
rank: 0, reads_size: 36
transfer kmers last round
Generating input for Zoltan
Total distinct (k+1)-mer entries across all proc's: 944
Number of batch iterations: 2
Total distinct k-mer bucket entries across all proc's: 111
Number of outliers (buckets discarded): 828
Number of buckets: 939
percentage of buckets discarded: 88.178914
Average time for Alltoall across all procs (secs): 0.000028
Average time for pack_sbuf_time across all procs (secs): 0.000000
Average time for AlltoallV across all procs (secs): 0.000022
Average time for unpack_rbuf_time across all procs (secs): 0.000199
Average time for unpack_rbuf_time:sort across all procs (secs): 0.000022
Average time for unpack_rbuf_time:insert across all procs (secs): 0.000020
Average time for unpack_rbuf_time:acc across all procs (secs): 0.000090
Average time for partial Sliding Window across all procs (secs): 0.000000
Average time for Sliding Window across all procs (secs): 0.000752
Average time for lmer Sliding Window across all procs (secs): 0.000009
Average time for vector insert across all procs (secs): 0.000025
Average time for temp_map insert across all procs (secs): 0.000075
Average number of re-calculations of min l-mer across all procs: 628.750000
Average time for performing k-mer bucketing across all procs (secs): 0.000900
Maximum time for performing k-mer bucketing across all procs (secs): 0.000906
Total number of pins: 329
Number of parts given to Zoltan: 4
Average time for partitioning across all procs (secs):   0.001530
Maximum time for partitioning across all procs (secs):   0.001530
Average time for data migration across all procs (secs): 0.000036
Maximum time for data migration across all procs (secs): 0.000036
Total number of reads partitioned across 4 procs: 11
Total number of parts across 4 procs: 4
```

## Primary tested compilers and architectures

* GCC 8.3.0, cray-mpich 7.7.10, x86, SUSE Linux Enterprise Server 15

## Acknowledgements

* Supported by U.S. National Science Foundation (awards: CCF 1946752, 1919122, 1919021)
