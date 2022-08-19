/*
 * @HEADER
 *
 * Modified from Zoltan's example code.
 * The function is created to run with the hypergraph already partitioned and distributed
 * among the MPI procs.
 *
 * @HEADER
 */

#ifndef ZOLTAN_PARTITION_H
#define ZOLTAN_PARTITION_H

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include "zoltan.h"
#include "defs.h"

#ifdef DEBUG
    template<typename S>
    void checkPtrs(std::vector<S> & ptrs, size_t expected_elements, S expected_min, S expected_max) {
        if (ptrs.size() != expected_elements) {
            std::cerr << "ptrs.size() != expected_elements : " << ptrs.size() << " " << expected_elements << std::endl;
            MPI_Abort(MPI_COMM_WORLD, -99);
        }
        if (ptrs[0] != expected_min) {
            std::cerr << "ptrs[0] != " << expected_min << ": " << ptrs[0] << std::endl;
            MPI_Abort(MPI_COMM_WORLD, -99);
        }
        for (size_t i=1; i<expected_elements; i++)
            if (ptrs[i] < ptrs[i-1]) {
                std::cerr << "(ptrs[i] < ptrs[i-1]) (i=" << i << ") : " << ptrs[i] << "  " << ptrs[i-1] << std::endl;
                MPI_Abort(MPI_COMM_WORLD, -99);
            }
        if (ptrs[expected_elements-1] != expected_max) {
            std::cerr << "ptrs[expected_elements-1] != expected_max : " << ptrs[expected_elements-1] << " " << expected_max << std::endl;
            MPI_Abort(MPI_COMM_WORLD, -99);
        }
    }
#endif
class HGRAPH_DATA {
    public:
        int numMyVertices;  /* number of vertices that I own initially */
        ZOLTAN_ID_TYPE startingGID;        /* starting global ID of these vertices, range [startingGID,  startingGID + numMyVertice) is the all vertices GID*/
        std::vector<size_t> & vtxDataPtrs; /* The vertex data lengths. It contains numMyVertices+1 elements, with vtxDataPtrs[i] and vtxDataPtrs[i+1] representing the start and ending positions of data for the ith vertex*/
        std::vector<read_t> & vtxData; /* The vertex data */

        std::vector<ZOLTAN_ID_TYPE> & edgeGID;       /* global ID of each of my hyperedges */
        std::vector<int> & nborIndex;     /* index into nborGID array of edge's vertices */
        std::vector<ZOLTAN_ID_TYPE> & nborGID;  /* Vertices of edge edgeGID[i] begin at nborGID[nborIndex[i]] */

        // Output for partitioning results
        std::vector<size_t> partPtr; /* indicate the starting position of data for each part in vtxData */
        std::vector<int> myParts; /* indicate what parts (global part ID) this proc has, and it's useful when the number of process is more than the number of parts*/

        HGRAPH_DATA (int _numMyVertices, ZOLTAN_ID_TYPE _startingGID, std::vector<size_t> & _vtxDataPtrs, std::vector<read_t> & _vtxData,
                std::vector<ZOLTAN_ID_TYPE> & _edgeGID, std::vector<int> & _nborIndex, std::vector<ZOLTAN_ID_TYPE> & _nborGID)
            : numMyVertices(_numMyVertices), startingGID(_startingGID), vtxDataPtrs(_vtxDataPtrs), vtxData(_vtxData),
            edgeGID(_edgeGID), nborIndex(_nborIndex), nborGID(_nborGID)
    {
#ifdef DEBUG
        checkPtrs(vtxDataPtrs, (size_t)numMyVertices+1, (size_t)0, (size_t)vtxData.size());
        checkPtrs(nborIndex, (size_t)edgeGID.size()+1, (int)0, (int)nborGID.size());
#endif
        partPtr.resize(2);
        partPtr[0] = 0;
        partPtr[1] = vtxData.size();
    }

        ~HGRAPH_DATA() { partPtr.clear(); myParts.clear(); }
        HGRAPH_DATA(const HGRAPH_DATA & ) = delete;
        HGRAPH_DATA(HGRAPH_DATA && ) = delete;
        HGRAPH_DATA& operator=(const HGRAPH_DATA &) = delete;
        HGRAPH_DATA& operator=(HGRAPH_DATA &&) = delete;

        size_t getNumMyParts() { return partPtr.size() - 1; }
        void getPartPos( size_t partID, size_t & start, size_t & end )
        {
#ifdef DEBUG
            if (partID >= getNumMyParts()) {
                std::cerr << "partID input is bigger than the number of parts (" << getNumMyParts() << "): " << partID << std::endl;
                MPI_Abort(MPI_COMM_WORLD, -99);
            }
#endif
            start = partPtr[partID];
            end = partPtr[partID + 1];
        }

        /* number of my hyperedges */
        int getNumMyHEdges() { return edgeGID.size(); }
        /* number of vertices in my hyperedges */
        int getNumAllNbors() { return nborGID.size(); }
        /* Starting position for the vertex i*/
        size_t getVtxDataStart(int i)
        {
            if (i >= numMyVertices) {
                std::cerr << "vertex ID is higher than the number of vertices (" << numMyVertices << "): " << i << std::endl;
                MPI_Abort(MPI_COMM_WORLD, -99);
            }
            return vtxDataPtrs[i];
        }
        /* Ending position for the vertex i*/
        size_t getVtxDataEnd(int i)
        {
            if (i >= numMyVertices) {
                std::cerr << "vertex ID is higher than the number of vertices (" << numMyVertices << "): " << i << std::endl;
                MPI_Abort(MPI_COMM_WORLD, -99);
            }
            return vtxDataPtrs[i+1];
        }
        /* The data size for the vertex i*/
        int getVtxDataSize(int i)
        {
            if (i >= numMyVertices) {
                std::cerr << "vertex ID is higher than the number of vertices (" << numMyVertices << "): " << i << std::endl;
                MPI_Abort(MPI_COMM_WORLD, -99);
            }
            return int(vtxDataPtrs[i+1] - vtxDataPtrs[i]);
        }
        int getVtxDataSizeByte(int i)  { return getVtxDataSize(i)*sizeof(read_t); }
};

class GRAPH_DATA {
    public:
        int numMyVertices;  /* number of vertices that I own initially */
        ZOLTAN_ID_TYPE startingGID;        /* starting global ID of these vertices, range [startingGID,  startingGID + numMyVertice) is the all vertices GID*/
        std::vector<size_t> & vtxDataPtrs; /* The vertex data lengths. It contains numMyVertices+1 elements, with vtxDataPtrs[i] and vtxDataPtrs[i+1] representing the start and ending positions of data for the ith vertex*/
        std::vector<read_t> & vtxData; /* The vertex data */

        std::vector<int> & nborIndex;     /* sizeof (getNumMyVertices + 1), index into nborGID array of vertices's neighbor list*/
        std::vector<ZOLTAN_ID_TYPE> & nborGID;  /* Neighbor vertices of vertex (i + startingGID) begins at nborGID[nborIndex[i]] ends at nborGID[nborIndex[i+1]] */

        // Output for partitioning results
        std::vector<size_t> partPtr; /* indicate the starting position of data for each part in vtxData */
        std::vector<int> myParts; /* indicate what parts this proc has, and it's useful when the number of process is more than the number of parts*/

        GRAPH_DATA (int _numMyVertices, ZOLTAN_ID_TYPE _startingGID, std::vector<size_t> & _vtxDataPtrs, std::vector<read_t> & _vtxData,
                std::vector<int> & _nborIndex, std::vector<ZOLTAN_ID_TYPE> & _nborGID)
            : numMyVertices(_numMyVertices), startingGID(_startingGID), vtxDataPtrs(_vtxDataPtrs), vtxData(_vtxData),
            nborIndex(_nborIndex), nborGID(_nborGID)
    {
#ifdef DEBUG
        checkPtrs(vtxDataPtrs, (size_t)numMyVertices+1, (size_t)0, (size_t)vtxData.size());
        checkPtrs(nborIndex, (size_t)numMyVertices+1, (int)0, (int)nborGID.size());
#endif
        partPtr.resize(2);
        partPtr[0] = 0;
        partPtr[1] = vtxData.size();
    }

        ~GRAPH_DATA() { partPtr.clear(); myParts.clear(); }
        GRAPH_DATA(const GRAPH_DATA & ) = delete;
        GRAPH_DATA(GRAPH_DATA && ) = delete;
        GRAPH_DATA& operator=(const GRAPH_DATA &) = delete;
        GRAPH_DATA& operator=(GRAPH_DATA &&) = delete;

        size_t getNumMyParts() { return partPtr.size() - 1; }
        void getPartPos( size_t partID, size_t & start, size_t & end )
        {
#ifdef DEBUG
            if (partID >= getNumMyParts()) {
                std::cerr << "partID input is bigger than the number of parts (" << getNumMyParts() << "): " << partID << std::endl;
                MPI_Abort(MPI_COMM_WORLD, -99);
            }
#endif
            start = partPtr[partID];
            end = partPtr[partID + 1];
        }

        /* number of edges */
        int getNumMyEdges() { return nborGID.size(); }
        /* Starting position for the vertex i*/
        size_t getVtxDataStart(int i)
        {
            if (i >= numMyVertices) {
                std::cerr << "vertex ID is higher than the number of vertices (" << numMyVertices << "): " << i << std::endl;
                MPI_Abort(MPI_COMM_WORLD, -99);
            }
            return vtxDataPtrs[i];
        }
        /* Ending position for the vertex i*/
        size_t getVtxDataEnd(int i)
        {
            if (i >= numMyVertices) {
                std::cerr << "vertex ID is higher than the number of vertices (" << numMyVertices << "): " << i << std::endl;
                MPI_Abort(MPI_COMM_WORLD, -99);
            }
            return vtxDataPtrs[i+1];
        }
        /* proc for a vertex*/
        int getProcID(ZOLTAN_ID_TYPE vgid) { return retrieve_read_proc_owner(vgid); }
        /* The data size for the vertex i*/
        int getVtxDataSize(int i)
        {
            if (i >= numMyVertices) {
                std::cerr << "vertex ID is higher than the number of vertices (" << numMyVertices << "): " << i << std::endl;
                MPI_Abort(MPI_COMM_WORLD, -99);
            }
            return int(vtxDataPtrs[i+1] - vtxDataPtrs[i]);
        }
        int getVtxDataSizeByte(int i)  { return getVtxDataSize(i)*sizeof(read_t); }
};

class ZDATA {
    public:
        HGRAPH_DATA & hg;

        int myRank;
        std::vector<ZOLTAN_ID_TYPE> vtxGID;  /* global ID of these vertices */

        ZDATA(HGRAPH_DATA & _hg, int _myRank) : hg(_hg), myRank(_myRank) {
            vtxGID.resize(hg.numMyVertices);
            for (size_t i=0; i<hg.numMyVertices; i++) vtxGID[i] = hg.startingGID + i;
            hg.myParts.resize(hg.numMyVertices); // It's used to store the part id for each vertex for now
            for (size_t i=0; i<hg.numMyVertices; i++) hg.myParts[i] = myRank;
        }

        ~ZDATA() { vtxGID.clear(); }
        ZDATA(const ZDATA & ) = delete;
        ZDATA(ZDATA && ) = delete;
        ZDATA& operator=(const ZDATA &) = delete;
        ZDATA& operator=(ZDATA &&) = delete;
};

class ZGRAPHDATA {
    public:
        GRAPH_DATA & g;

        int myRank;
        std::vector<ZOLTAN_ID_TYPE> vtxGID;  /* global ID of these vertices */

        ZGRAPHDATA(GRAPH_DATA & _g, int _myRank) : g(_g), myRank(_myRank) {
            vtxGID.resize(g.numMyVertices);
            for (size_t i=0; i<g.numMyVertices; i++) vtxGID[i] = g.startingGID + i;
            g.myParts.resize(g.numMyVertices); // It's used to store the part id for each vertex for now
            for (size_t i=0; i<g.numMyVertices; i++) g.myParts[i] = myRank;
        }

        ~ZGRAPHDATA() { vtxGID.clear(); }
        ZGRAPHDATA(const ZGRAPHDATA & ) = delete;
        ZGRAPHDATA(ZGRAPHDATA && ) = delete;
        ZGRAPHDATA& operator=(const ZGRAPHDATA &) = delete;
        ZGRAPHDATA& operator=(ZGRAPHDATA &&) = delete;
};

class ZOLTAN_VERTICES_EXCHANGE {
    public:
        int num;
        ZOLTAN_ID_PTR globalGids;
        ZOLTAN_ID_PTR localGids;
        int *procs, *toPart;

        ZOLTAN_VERTICES_EXCHANGE() {}
        ~ZOLTAN_VERTICES_EXCHANGE() {
            num=0;
            Zoltan_LB_Free_Part(&globalGids, &localGids,
                    &procs, &toPart);
        }
        ZOLTAN_VERTICES_EXCHANGE(const ZOLTAN_VERTICES_EXCHANGE & ) = delete;
        ZOLTAN_VERTICES_EXCHANGE(ZOLTAN_VERTICES_EXCHANGE && ) = delete;
        ZOLTAN_VERTICES_EXCHANGE& operator=(const ZOLTAN_VERTICES_EXCHANGE &) = delete;
        ZOLTAN_VERTICES_EXCHANGE& operator=(ZOLTAN_VERTICES_EXCHANGE &&) = delete;
};

/*
 * Output the hgraph vtxGID after partitioning
 * */
void print_hgraph_read_gid (ZDATA & zdata);

/*
 * Output the graph vtxGID after partitioning
 * */
void print_graph_read_gid (ZGRAPHDATA & zgraphdata);

/*
 * calling zoltan to partition the vertices
 * Inputs:
 *     argc, argv: inputs for the executable to initalize Zoltan
 *     hg: structure stores the data for the hypergraph located in current
 *     mpi proc, including data for the vertices
 *     its vertices and the
 *     parts: the number of partitions
 */
int perform_zoltan_partition_hg (int argc, char **argv,
        HGRAPH_DATA & hg, const int parts);

/*
 * calling zoltan with ParMetis to partition the vertices
 * Inputs:
 *     argc, argv: inputs for the executable to initalize Zoltan
 *     g: structure stores the data for the graph located in current
 *     mpi proc, including data for the vertices
 *     its vertices and the
 *     parts: the number of partitions
 */
int perform_zoltan_partition_g (int argc, char **argv,
        GRAPH_DATA & g, const int parts);
#endif
