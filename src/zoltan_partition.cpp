/*
 * @HEADER
 *
 * Modified from Zoltan's example code.
 * The function is created to run with the hypergraph and graph already partitioned and distributed
 * among the MPI procs.
 *
 * @HEADER
 */

#include "zoltan_partition.h"
#include <map>

extern int rank, size;
extern std::string outputDirName;
extern unsigned int partitioning_seed;

/* Application defined query functions */

// For hypergraph partitioning
static int get_number_of_vertices(void *data, int *ierr)
{
    HGRAPH_DATA *hg = (HGRAPH_DATA *)data;
    *ierr = ZOLTAN_OK;
    return hg->numMyVertices;
}

static void get_vertex_list(void *data, int sizeGID, int sizeLID,
        ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
        int wgt_dim, float *obj_wgts, int *ierr)
{
    HGRAPH_DATA *hg= (HGRAPH_DATA *)data;
    *ierr = ZOLTAN_OK;

    /* In this example, return the IDs of our vertices, but no weights.
     * Zoltan will assume equally weighted vertices.
     */

    for (int i=0; i<hg->numMyVertices; i++){
        globalID[i] = hg->startingGID + i;
        localID[i] = i;
    }
}

#ifdef EVAL_ZOLTAN
static void get_vertex_list_post_migrate(void *data, int sizeGID, int sizeLID,
        ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
        int wgt_dim, float *obj_wgts, int *ierr)
{
    ZDATA *zdata = (ZDATA *)data;
    *ierr = ZOLTAN_OK;

    /* In this example, return the IDs of our vertices, but no weights.
     * Zoltan will assume equally weighted vertices.
     */

    for (int i=0; i<zdata->hg.numMyVertices; i++){
        globalID[i] = zdata->vtxGID[i];
        localID[i] = i;
    }
}
#endif

static void get_hypergraph_size(void *data, int *num_lists, int *num_nonzeroes,
        int *format, int *ierr)
{
    HGRAPH_DATA *hg = (HGRAPH_DATA *)data;
    *ierr = ZOLTAN_OK;

    *num_lists = hg->getNumMyHEdges();
    *num_nonzeroes = hg->getNumAllNbors();

    /* We will provide compressed hyperedge (row) format.  The alternative is
     * is compressed vertex (column) format: ZOLTAN_COMPRESSED_VERTEX.
     */

    *format = ZOLTAN_COMPRESSED_EDGE;

    return;
}

static void get_hypergraph(void *data, int sizeGID, int num_edges, int num_nonzeroes,
        int format, ZOLTAN_ID_PTR edgeGID, int *vtxPtr,
        ZOLTAN_ID_PTR vtxGID, int *ierr)
{
    int i;

    HGRAPH_DATA *hg = (HGRAPH_DATA *)data;
    *ierr = ZOLTAN_OK;

    if ( (num_edges != hg->getNumMyHEdges()) || (num_nonzeroes != hg->getNumAllNbors()) ||
            (format != ZOLTAN_COMPRESSED_EDGE)) {
        *ierr = ZOLTAN_FATAL;
        return;
    }

    for (i=0; i < num_edges; i++){
        edgeGID[i] = hg->edgeGID[i];
        vtxPtr[i] = hg->nborIndex[i];
    }

    for (i=0; i < num_nonzeroes; i++){
        vtxGID[i] = hg->nborGID[i];
    }

    return;
}

// For data migrate
static int get_datasize_of_vertices(
  void *data,
  int num_gid_entries,
  int num_lid_entries,
  ZOLTAN_ID_PTR global_id,
  ZOLTAN_ID_PTR local_id,
  int *ierr
)
{
    HGRAPH_DATA *hg= (HGRAPH_DATA *)data;
    *ierr = ZOLTAN_OK;
    ZOLTAN_ID_TYPE the_local_id = *local_id;
    if (the_local_id >= hg->numMyVertices) {
        printf("(the_local_id >= hg->numMyVertices)\n");
        *ierr = ZOLTAN_FATAL;
    }
    return hg->getVtxDataSizeByte(the_local_id);
}

static void get_data_of_vertices (
  void *data,
  int num_gid_entries,
  int num_lid_entries,
  ZOLTAN_ID_PTR global_id,
  ZOLTAN_ID_PTR local_id,
  int dest_proc,
  int size,
  char *buf,
  int *ierr
)
{
    HGRAPH_DATA *hg = (HGRAPH_DATA *)data;

    ZOLTAN_ID_TYPE the_local_id = *local_id;
    if (the_local_id >= hg->numMyVertices) {
        printf("(the_local_id >= hg->numMyVertices)\n");
        *ierr = ZOLTAN_FATAL;
    }
    if (size >= hg->getVtxDataSizeByte(the_local_id))
        *ierr = ZOLTAN_OK;
    else
        *ierr = ZOLTAN_FATAL;

    size_t start_pos = hg->getVtxDataStart(the_local_id);
    int datalength=hg->getVtxDataSize(the_local_id);
    for (int i=0; i<datalength; i++)
        ((read_t *)buf)[i] = hg->vtxData[i + start_pos];
}

// Step 1,2,3
// Move forward the vertices that stay in the current proc and disregard
// ones that are going be transferred to other proc
//
// Step 4, 5
// preallocate space for vertices that are coming to this proc, except
// vtxData array, since we don't know how many space by only looking at
// the IDs.
static void reorganize_vdataparts (
  void *data,
  int num_gid_entries,
  int num_lid_entries,
  int num_import,
  ZOLTAN_ID_PTR import_global_ids,
  ZOLTAN_ID_PTR import_local_ids,
  int *import_procs,
  int *import_to_part,
  int num_export,
  ZOLTAN_ID_PTR export_global_ids,
  ZOLTAN_ID_PTR export_local_ids,
  int *export_procs,
  int *export_to_part,
  int *ierr
)
{
    ZDATA *zdata = (ZDATA *)data;
    // only support num id entries equal to 1
    if (num_gid_entries == 1 || num_lid_entries == 1)
        *ierr = ZOLTAN_OK;
    else
        *ierr = ZOLTAN_FATAL;

    // 1. label ones to export and initial parts
    // init my parts as myRank

    bool export_sign[zdata->hg.numMyVertices] = {false};
    for (size_t i=0; i<num_export; i++) {
        // update signs for exporting to other procs
#ifdef DEBUG
        if (export_local_ids[i] >= zdata->hg.numMyVertices)
            *ierr = ZOLTAN_FATAL;
#endif
        if (export_procs[i] != zdata->myRank)
            export_sign[export_local_ids[i]] = true;
        else // update myparts that is within the procs
            zdata->hg.myParts[export_local_ids[i]] = export_to_part[i];
    }

    // 2. remove export ones
    // TODO: optimize for different cases
    size_t count=0;
    size_t pos=0;
    for (size_t i=0; i<zdata->hg.numMyVertices; i++) {
        if (!export_sign[i]) {
            zdata->vtxGID[count] = zdata->vtxGID[i];
            zdata->hg.myParts[count] = zdata->hg.myParts[i];
            size_t vtxidatastart=zdata->hg.getVtxDataStart(i);
            size_t vtxidataend=zdata->hg.getVtxDataEnd(i);
            // set the starting position of the vertex moved
            zdata->hg.vtxDataPtrs[count]=pos;
            for (size_t k=vtxidatastart; k<vtxidataend; pos++, k++)
                zdata->hg.vtxData[pos] = zdata->hg.vtxData[k];
            count++;
        }
    }
    // set the ending position of the last vertex remains
    zdata->hg.vtxDataPtrs[count]=pos;
    zdata->hg.vtxData.resize(pos); // shrink vtxData to remove unused memory

    // 3. update num of vertices after removing exported ones
    zdata->hg.numMyVertices = count;

    // 4. count the number of vertices to be imported into myRank only for
    // preallocation
    for (size_t i=0; i<num_import; i++) {
        if (import_procs[i] != zdata->myRank)
            count++;
    }

    // 5. Preallocate vertex data ptrs, vtxGID, and myParts for receiving
    // imported data
    // for vtxData, we have to resize for each vertex, since we don't know
    // the size of the data right now
    zdata->hg.vtxDataPtrs.resize(count+1);
    zdata->vtxGID.resize(count);
    zdata->hg.myParts.resize(count);
}

// Append the data for a vertex imported from other procs
static void append_data_imported(
  void *data,
  int num_gid_entries,
  ZOLTAN_ID_PTR global_id,
  int size,
  char *buf,
  int *ierr
)
{
    ZDATA *zdata = (ZDATA *)data;

    size_t localvid = zdata->hg.numMyVertices;
    zdata->hg.numMyVertices++;

    zdata->vtxGID[localvid] = *global_id;
    zdata->hg.myParts[localvid] = -1; // as the part id is not sent (I guess the reason is to avoid unnecessary communication)
    size_t start_pos = zdata->hg.vtxData.size();
    // the starting position of the current vertex should be set already
    if (start_pos != zdata->hg.vtxDataPtrs[localvid]) {
        printf("start_pos != zdata->hg.vtxDataPtrs[localvid]\n");
        *ierr = ZOLTAN_FATAL;
    }
    zdata->hg.vtxData.resize(start_pos + size/sizeof(read_t));
    for (int i=0; i<size/sizeof(read_t); i++)
        zdata->hg.vtxData[start_pos + i] = ((read_t *)buf)[i];
    // set the ending position of current vertex (and the start position
    // of the next vertex)
    zdata->hg.vtxDataPtrs[localvid+1] = start_pos + size/sizeof(read_t);
}

// Step 1,2
// Assigne part id for vertices imported and appended in the vertex array,
// but with part id as -1 (as the part ids are unknown when recieving the
// vertex global id and vertex data)
// Step 3,4,5
// Sort vertices by part id, involving updating vtxDataPtrs, vtxGID, vtxData
static void sort_byparts(
  void *data,
  int num_gid_entries,
  int num_lid_entries,
  int num_import,
  ZOLTAN_ID_PTR import_global_ids,
  ZOLTAN_ID_PTR import_local_ids,
  int *import_procs,
  int *import_to_part,
  int num_export,
  ZOLTAN_ID_PTR export_global_ids,
  ZOLTAN_ID_PTR export_local_ids,
  int *export_procs,
  int *export_to_part,
  int *ierr
)
{
    ZDATA *zdata = (ZDATA *)data;
    *ierr = ZOLTAN_OK;

    // 1. build map between global vertex id to the part assign
    size_t num_v_imported = 0;
    std::map<ZOLTAN_ID_TYPE, int> m;
    for (size_t i=0; i<num_import; i++) {
        if (import_procs[i] != zdata->myRank) {
            num_v_imported++;
            m.insert(std::pair<ZOLTAN_ID_TYPE, int>(import_global_ids[i], import_to_part[i]));
        }
    }
    if (num_v_imported != m.size()) {
        printf("num_v_imported != m.size()\n");
        *ierr = ZOLTAN_FATAL;
    }

    // 2. assign part ids for vertices imported and without part id
    for(size_t pos = zdata->hg.numMyVertices - num_v_imported; pos < zdata->hg.numMyVertices; pos++) {
        if (zdata->hg.myParts[pos] != -1) {
            printf("zdata->hg.myParts[pos] != -1\n");
            *ierr = ZOLTAN_FATAL;
        }
        if (m.find(zdata->vtxGID[pos]) != m.end())
            zdata->hg.myParts[pos] = m[zdata->vtxGID[pos]];
        else {
            printf("(m.find(zdata->vtxGID[pos]) != m.end())\n");
            *ierr = ZOLTAN_FATAL;
        }
    }
    m.clear();

    // 3. generate maps between `part id` and its parts_sizes: `number of vertices` and `total vertices data size`
    std::map<int, std::pair<size_t, size_t>> parts_sizes;
    for(size_t i=0; i<zdata->hg.numMyVertices; i++) {
        if (parts_sizes.find(zdata->hg.myParts[i]) == parts_sizes.end()) {
            parts_sizes.insert(
                    std::pair<int, std::pair<size_t, size_t>>(
                        zdata->hg.myParts[i], std::pair<size_t, size_t>(1, zdata->hg.getVtxDataSize(i))));
        } else {
            parts_sizes[zdata->hg.myParts[i]].first += 1;
            parts_sizes[zdata->hg.myParts[i]].second += zdata->hg.getVtxDataSize(i);
        }
    }

    // 4. update sizes to ptr by prefix-summing them
    std::pair<size_t, size_t> p(0, 0);
    for (const auto &[key, value] : parts_sizes) {
        p.first += value.first;
        p.second += value.second;
        parts_sizes[key].first = p.first - parts_sizes[key].first;
        parts_sizes[key].second = p.second - parts_sizes[key].second;
    }

    // 5. sort vtxDataPtrs, vtxGID, vtxData by part ids
    //    - first back up the three arrays
    //    - then read from the back ups, and write to the zdata based on
    //    the parts_sizes
    std::vector<ZOLTAN_ID_TYPE> vtxGID_copy(zdata->vtxGID);
    std::vector<size_t> vtxDataPtrs_copy(zdata->hg.vtxDataPtrs);
    std::vector<read_t> vtxData_copy(zdata->hg.vtxData);
#ifdef DEBUG
    std::vector<int> myParts_copy(zdata->hg.myParts);
#endif

    for (size_t i=0; i<zdata->hg.numMyVertices; i++) {
#ifdef DEBUG
        int mypart=myParts_copy[i];
#else
        int mypart=zdata->hg.myParts[i];
#endif
        size_t newLocalID=parts_sizes[mypart].first;
        size_t newDataStart=parts_sizes[mypart].second;
        parts_sizes[mypart].first += 1;
        parts_sizes[mypart].second += vtxDataPtrs_copy[i+1]-vtxDataPtrs_copy[i];

        zdata->vtxGID[newLocalID] = vtxGID_copy[i];
        // the last element of vtxDataPtrs should stay the same and does not need to be assigned
        zdata->hg.vtxDataPtrs[newLocalID] = newDataStart;
#ifdef DEBUG
        zdata->hg.myParts[newLocalID] = myParts_copy[i];
#endif
        if (newDataStart != vtxDataPtrs_copy[i]) {
            for (size_t t=newDataStart, f=vtxDataPtrs_copy[i]; f<vtxDataPtrs_copy[i+1]; t++, f++)
                zdata->hg.vtxData[t] = vtxData_copy[f];
        }
    }
    vtxGID_copy.clear();
    vtxDataPtrs_copy.clear();
    vtxData_copy.clear();

#ifdef DEBUG
    myParts_copy.clear();
    int previous_part_id = -1;
    std::vector<size_t> partPtr_new;
    std::vector<int> myParts_new;
    for (size_t i=0; i<zdata->hg.numMyVertices; i++) {
        if (zdata->hg.myParts[i] != previous_part_id) {
            myParts_new.push_back(zdata->hg.myParts[i]);
            partPtr_new.push_back(i);
        }
        previous_part_id = zdata->hg.myParts[i];
    }
    partPtr_new.push_back(zdata->hg.numMyVertices);
#endif

    // rearrange myParts and malloc partPtr
    zdata->hg.partPtr.resize(parts_sizes.size()+1);
    size_t count=0;
    zdata->hg.partPtr[count] = 0;
    zdata->hg.myParts.resize(parts_sizes.size());
    for (const auto &[key, value] : parts_sizes) {
        zdata->hg.myParts[count] = key;
        zdata->hg.partPtr[count+1] = value.first;
        count++;
    }
#ifdef DEBUG
    if (partPtr_new.size() != zdata->hg.partPtr.size() || myParts_new.size() != zdata->hg.myParts.size()) {
        printf("-- partPtr_new.size() {%ld} != zdata->hg.partPtr.size() {%ld} || myParts_new.size() {%ld} != zdata->hg.myParts.size() {%ld}\n", partPtr_new.size(), zdata->hg.partPtr.size(), myParts_new.size(), zdata->hg.myParts.size());
        *ierr = ZOLTAN_FATAL;
    } else {
        for (size_t i=0; i<myParts_new.size(); i++) {
            if (myParts_new[i] != zdata->hg.myParts[i]) {
                printf("myParts_new[%ld] {%d} != zdata->hg.myParts[%ld] {%d}\n", i, myParts_new[i], i, zdata->hg.myParts[i]);
                *ierr = ZOLTAN_FATAL;
            }
        }
        for (size_t i=0; i<partPtr_new.size(); i++) {
            if (partPtr_new[i] != zdata->hg.partPtr[i]) {
                printf("partPtr_new[%ld] {%d} != zdata->hg.partPtr[%ld] {%d}\n", i, partPtr_new[i], i, zdata->hg.partPtr[i]);
                *ierr = ZOLTAN_FATAL;
            }
        }
    }
    partPtr_new.clear();
    myParts_new.clear();
#endif
}

void print_hgraph_read_gid (ZDATA & zdata)
{
    size_t numMyParts = zdata.hg.getNumMyParts();
    for (size_t partID = 0; partID < numMyParts; partID++) {
        std::string output_filename(outputDirName + "/vtxGID_r" + std::to_string(rank) + "_p" + std::to_string(zdata.hg.myParts[partID]) + ".fa");
        FILE *f = fopen(output_filename.c_str(), "w");
        fprintf(f, "%lu\n", zdata.hg.partPtr[partID+1] - zdata.hg.partPtr[partID]);

        for (int i=zdata.hg.partPtr[partID]; i<zdata.hg.partPtr[partID+1]; i++)
            fprintf(f, "%lu\n", zdata.vtxGID[i]);

        fclose(f);
    }
}

/*
 * calling zoltan to partition the vertices
 * Inputs:
 *     argc, argv: inputs for the executable to initalize Zoltan
 *     hg: structure stores the data for the hypergraph located in current
 *     mpi proc, including data for the vertices
 *     parts: the number of partitions
 */
int perform_zoltan_partition_hg (int argc, char **argv,
        HGRAPH_DATA & hg, const int parts)
{
    int rc;
    float version;
    struct Zoltan_Struct *zz;
    int changes, numGidEntries, numLidEntries;
    ZOLTAN_VERTICES_EXCHANGE in, out;

    // Logging the total number of pins (size of the hypergraph)
    size_t mynumofpins = hg.getNumAllNbors();
    size_t totalnumofpins;
    MPI_Reduce(&mynumofpins, &totalnumofpins, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        printf("Total number of pins: %lu\n", totalnumofpins);
    }

#ifdef TIME_PROFILE
    double partition_time, migrate_time;
    double begin_time, end_time;
    MPI_Barrier(MPI_COMM_WORLD);
    begin_time = MPI_Wtime ();
#endif

    /******************************************************************
     ** Initialize Zoltan
     ******************************************************************/
#ifdef DEBUG
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
        printf("init zoltan\n");
#endif
    rc = Zoltan_Initialize(argc, argv, &version);
    if (rc != ZOLTAN_OK){
        fprintf(stderr, "Zoltan partitioner: Failed to initialize");
        MPI_Finalize();
        Zoltan_Destroy(&zz);
        exit(0);
    }

    /******************************************************************
     ** Create a Zoltan library structure for this instance of load
     ** balancing.  Set the parameters and query functions that will
     ** govern the library's calculation.  See the Zoltan User's
     ** Guide for the definition of these and many other parameters.
     ******************************************************************/
#ifdef DEBUG
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
        printf("create zoltan struct\n");
#endif
    zz = Zoltan_Create(MPI_COMM_WORLD);

#ifdef DEBUG
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
        printf("set parameter\n");
#endif
    /* General parameters */
#ifdef DEBUG
    Zoltan_Set_Param(zz, "DEBUG_LEVEL", "2");
#else
    Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
#endif
    //TODO: change based on partitioner input
    Zoltan_Set_Param(zz, "LB_METHOD", "HYPERGRAPH");   /* partitioning method */
    Zoltan_Set_Param(zz, "HYPERGRAPH_PACKAGE", "PHG"); /* version of method */
    Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");/* global IDs are integers */
    Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");/* local IDs are integers */
    Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL"); /* export AND import lists */
    Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "0"); /* use Zoltan default vertex weights */
    Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM", "0");/* use Zoltan default hyperedge weights */
    Zoltan_Set_Param(zz, "IMBALANCE_TOL", "1.01"); /* the load imbalance setting, the upper bound load is 1.01*(averageload) */
    char param_str[10];
    sprintf(param_str,"%d",parts);
    Zoltan_Set_Param(zz, "NUM_GLOBAL_PARTS", param_str);/* partition the vertices to `parts` parts*/
    if (rank == 0) printf("Number of parts given to Zoltan: %s\n", param_str);
    if (partitioning_seed != 0) {
        sprintf(param_str,"%d",partitioning_seed);
        Zoltan_Set_Param(zz, "SEED", param_str); /* the seeds for partitioning */
        if (rank == 0) printf("Set seeds for Zoltan as: %s\n", param_str);
    }

    /* PHG parameters  - see the Zoltan User's Guide for many more
     *   (The "PARTITION" approach asks Zoltan to partition from scratch.)
     */
    Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");

    /* Application defined query functions */
    Zoltan_Set_Num_Obj_Fn(zz, get_number_of_vertices, &hg);
    Zoltan_Set_Obj_List_Fn(zz, get_vertex_list, &hg);
    Zoltan_Set_HG_Size_CS_Fn(zz, get_hypergraph_size, &hg);
    Zoltan_Set_HG_CS_Fn(zz, get_hypergraph, &hg);

    /******************************************************************
     ** Zoltan can now partition the vertices of hypergraph.
     ******************************************************************/

#ifdef DEBUG
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
        printf("start lb partitioning\n");
#endif

    rc = Zoltan_LB_Partition(zz, /* input (all remaining fields are output) */
            &changes,        /* 1 if partitioning was changed, 0 otherwise */
            &numGidEntries,  /* Number of integers used for a global ID */
            &numLidEntries,  /* Number of integers used for a local ID */
            &in.num,      /* Number of vertices to be sent to me */
            &in.globalGids,  /* Global IDs of vertices to be sent to me */
            &in.localGids,   /* Local IDs of vertices to be sent to me */
            &in.procs,    /* Process rank for source of each incoming vertex */
            &in.toPart,   /* New partition for each incoming vertex */
            &out.num,      /* Number of vertices I must send to other processes*/
            &out.globalGids,  /* Global IDs of the vertices I must send */
            &out.localGids,   /* Local IDs of the vertices I must send */
            &out.procs,    /* Process to which I send each of the vertices */
            &out.toPart);  /* Partition to which each vertex will belong */

    if (rc != ZOLTAN_OK){
        fprintf(stderr, "Zoltan partitioner: Failed to partition");
        MPI_Finalize();
        Zoltan_Destroy(&zz);
        exit(0);
    }

#ifdef TIME_PROFILE
    end_time = MPI_Wtime ();
    MPI_Barrier(MPI_COMM_WORLD);
    partition_time = end_time - begin_time;
    begin_time = MPI_Wtime ();
#endif

#ifdef DEBUG
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
        printf("start set parameters for migrate\n");
#endif
    /* Pack and Unpack function for anto migration */
    Zoltan_Set_Obj_Size_Fn(zz, get_datasize_of_vertices, &hg);
    Zoltan_Set_Pack_Obj_Fn(zz, get_data_of_vertices, &hg);
    ZDATA zdata(hg, rank);
    Zoltan_Set_Mid_Migrate_PP_Fn(zz, reorganize_vdataparts, &zdata);
    Zoltan_Set_Unpack_Obj_Fn(zz, append_data_imported, &zdata);
    Zoltan_Set_Post_Migrate_PP_Fn(zz, sort_byparts, &zdata);
#ifdef DEBUG
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
        printf("start migrate\n");
#endif

    // Data migration
    rc = Zoltan_Migrate(zz, /* input (all remaining fields are output) */
            in.num,      /* Number of vertices to be sent to me */
            in.globalGids,  /* Global IDs of vertices to be sent to me */
            in.localGids,   /* Local IDs of vertices to be sent to me */
            in.procs,    /* Process rank for source of each incoming vertex */
            in.toPart,   /* New partition for each incoming vertex */
            out.num,      /* Number of vertices I must send to other processes*/
            out.globalGids,  /* Global IDs of the vertices I must send */
            out.localGids,   /* Local IDs of the vertices I must send */
            out.procs,    /* Process to which I send each of the vertices */
            out.toPart);  /* Partition to which each vertex will belong */
    if (rc != ZOLTAN_OK){
        fprintf(stderr, "Zoltan partitioner: Failed to migrate data");
        MPI_Finalize();
        Zoltan_Destroy(&zz);
        exit(0);
    }
#ifdef TIME_PROFILE
    end_time = MPI_Wtime ();
    MPI_Barrier(MPI_COMM_WORLD);
    migrate_time = end_time - begin_time;
#endif

#ifdef EVAL_ZOLTAN
    // NOTE: when the number of process is not equal to the number of
    // parts, the output is incorrect, as Zoltan_LB_Eval only recognizes
    // the number of parts as the number of process.
    Zoltan_Set_Obj_List_Fn(zz, get_vertex_list_post_migrate, &zdata);
    int res = Zoltan_LB_Eval(zz, 1, NULL, NULL, NULL);
    if (res) printf("Warning: Zoltan_LB_Eval returns non-zero %d\n", res);
#endif

#ifdef DEBUG
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
        printf("start cleanup\n");

    print_hgraph_read_gid(zdata);
#endif

    Zoltan_Destroy(&zz);

#ifdef TIME_PROFILE
    double global_partition_time_sum, global_migrate_time_sum;
    double global_partition_time_max, global_migrate_time_max;

    MPI_Reduce(&partition_time, &global_partition_time_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&partition_time, &global_partition_time_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&migrate_time, &global_migrate_time_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&migrate_time, &global_migrate_time_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        printf ("Average time for partitioning across all procs (secs):   %f \n",
                (double)global_partition_time_sum/(double)size);
        printf ("Maximum time for partitioning across all procs (secs):   %f \n",
                (double)global_partition_time_max);
        printf ("Average time for data migration across all procs (secs): %f \n",
                (double)global_migrate_time_sum/(double)size);
        printf ("Maximum time for data migration across all procs (secs): %f \n",
                (double)global_migrate_time_max);
    }
#endif

    return 0;
}

// TODO: reduce redandent code, combine with Hypergraph case
/* Application defined query functions */
// For graph partitioning
static int get_number_of_vertices_g(void *data, int *ierr)
{
    GRAPH_DATA *g = (GRAPH_DATA *)data;
    *ierr = ZOLTAN_OK;
    return g->numMyVertices;
}

static void get_vertex_list_g(void *data, int sizeGID, int sizeLID,
        ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
        int wgt_dim, float *obj_wgts, int *ierr)
{
    GRAPH_DATA *g= (GRAPH_DATA *)data;
    *ierr = ZOLTAN_OK;

    /* In this example, return the IDs of our vertices, but no weights.
     * Zoltan will assume equally weighted vertices.
     */

    for (int i=0; i<g->numMyVertices; i++){
        globalID[i] = g->startingGID + i;
        localID[i] = i;
    }
}

#ifdef EVAL_ZOLTAN
static void get_vertex_list_post_migrate_g(void *data, int sizeGID, int sizeLID,
        ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
        int wgt_dim, float *obj_wgts, int *ierr)
{
    ZGRAPHDATA *zdata = (ZGRAPHDATA *)data;
    *ierr = ZOLTAN_OK;

    /* In this example, return the IDs of our vertices, but no weights.
     * Zoltan will assume equally weighted vertices.
     */

    for (int i=0; i<zdata->g.numMyVertices; i++){
        globalID[i] = zdata->vtxGID[i];
        localID[i] = i;
    }
}
#endif

static void get_edge_list_size_g(void *data, int sizeGID, int sizeLID,
                      int num_obj,
             ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
             int *numEdges, int *ierr)
{
    GRAPH_DATA *g = (GRAPH_DATA *)data;
    if ( (sizeGID != 1) || (sizeLID != 1) || (num_obj != g->numMyVertices) ){
        *ierr = ZOLTAN_FATAL;
        return;
    }

    for (int i=0;  i < num_obj ; i++){
        ZOLTAN_ID_TYPE idx = localID[i];
        numEdges[i] = g->nborIndex[idx+1] - g->nborIndex[idx];
    }

    *ierr = ZOLTAN_OK;
    return;
}

static void get_edge_list_g(void *data, int sizeGID, int sizeLID,
        int num_obj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
        int *num_edges,
        ZOLTAN_ID_PTR nborGID, int *nborProc,
        int wgt_dim, float *ewgts, int *ierr)
{
    int i, j, from, to;
    int *nextProc;
    ZOLTAN_ID_TYPE *nextNbor;

    GRAPH_DATA *g= (GRAPH_DATA *)data;

    // TODO: change to support weight
    if ( (sizeGID != 1) || (sizeLID != 1) ||
            (num_obj != g->numMyVertices )||
            (wgt_dim != 0)){
        *ierr = ZOLTAN_FATAL;
        return;
    }

    nextNbor = nborGID;
    nextProc = nborProc;

    for (i=0; i < num_obj; i++){

        /*
         * In this example, we are not setting edge weights.  Zoltan will
         * set each edge to weight 1.0.
         */

        to = g->nborIndex[localID[i]+1];
        from = g->nborIndex[localID[i]];
        if ((to - from) != num_edges[i]){
            *ierr = ZOLTAN_FATAL;
            return;
        }

        for (j=from; j < to; j++){
            *nextNbor++ = g->nborGID[j];
            *nextProc++ = g->getProcID(g->nborGID[j]);
        }
    }

    *ierr = ZOLTAN_OK;
    return;
}

// For data migrate
static int get_datasize_of_vertices_g(
  void *data,
  int num_gid_entries,
  int num_lid_entries,
  ZOLTAN_ID_PTR global_id,
  ZOLTAN_ID_PTR local_id,
  int *ierr
)
{
    GRAPH_DATA *g= (GRAPH_DATA *)data;
    *ierr = ZOLTAN_OK;
    ZOLTAN_ID_TYPE the_local_id = *local_id;
    if (the_local_id >= g->numMyVertices) {
        printf("(the_local_id >= g->numMyVertices)\n");
        *ierr = ZOLTAN_FATAL;
    }
    return g->getVtxDataSizeByte(the_local_id);
}

static void get_data_of_vertices_g (
  void *data,
  int num_gid_entries,
  int num_lid_entries,
  ZOLTAN_ID_PTR global_id,
  ZOLTAN_ID_PTR local_id,
  int dest_proc,
  int size,
  char *buf,
  int *ierr
)
{
    GRAPH_DATA *g = (GRAPH_DATA *)data;

    ZOLTAN_ID_TYPE the_local_id = *local_id;
    if (the_local_id >= g->numMyVertices) {
        printf("(the_local_id >= g->numMyVertices)\n");
        *ierr = ZOLTAN_FATAL;
    }
    if (size >= g->getVtxDataSizeByte(the_local_id))
        *ierr = ZOLTAN_OK;
    else
        *ierr = ZOLTAN_FATAL;

    size_t start_pos = g->getVtxDataStart(the_local_id);
    int datalength=g->getVtxDataSize(the_local_id);
    for (int i=0; i<datalength; i++)
        ((read_t *)buf)[i] = g->vtxData[i + start_pos];
}

// Step 1,2,3
// Move forward the vertices that stay in the current proc and disregard
// ones that are going be transferred to other proc
//
// Step 4, 5
// preallocate space for vertices that are coming to this proc, except
// vtxData array, since we don't know how many space by only looking at
// the IDs.
static void reorganize_vdataparts_g (
  void *data,
  int num_gid_entries,
  int num_lid_entries,
  int num_import,
  ZOLTAN_ID_PTR import_global_ids,
  ZOLTAN_ID_PTR import_local_ids,
  int *import_procs,
  int *import_to_part,
  int num_export,
  ZOLTAN_ID_PTR export_global_ids,
  ZOLTAN_ID_PTR export_local_ids,
  int *export_procs,
  int *export_to_part,
  int *ierr
)
{
    ZGRAPHDATA *zdata = (ZGRAPHDATA *)data;
    // only support num id entries equal to 1
    if (num_gid_entries == 1 || num_lid_entries == 1)
        *ierr = ZOLTAN_OK;
    else
        *ierr = ZOLTAN_FATAL;

    // 1. label ones to export and initial parts
    // init my parts as myRank

    bool export_sign[zdata->g.numMyVertices] = {false};
    for (size_t i=0; i<num_export; i++) {
        // update signs for exporting to other procs
#ifdef DEBUG
        if (export_local_ids[i] >= zdata->g.numMyVertices)
            *ierr = ZOLTAN_FATAL;
#endif
        if (export_procs[i] != zdata->myRank)
            export_sign[export_local_ids[i]] = true;
        else // update myparts that is within the procs
            zdata->g.myParts[export_local_ids[i]] = export_to_part[i];
    }

    // 2. remove export ones
    // TODO: optimize for different cases
    size_t count=0;
    size_t pos=0;
    for (size_t i=0; i<zdata->g.numMyVertices; i++) {
        if (!export_sign[i]) {
            zdata->vtxGID[count] = zdata->vtxGID[i];
            zdata->g.myParts[count] = zdata->g.myParts[i];
            size_t vtxidatastart=zdata->g.getVtxDataStart(i);
            size_t vtxidataend=zdata->g.getVtxDataEnd(i);
            // set the starting position of the vertex moved
            zdata->g.vtxDataPtrs[count]=pos;
            for (size_t k=vtxidatastart; k<vtxidataend; pos++, k++)
                zdata->g.vtxData[pos] = zdata->g.vtxData[k];
            count++;
        }
    }
    // set the ending position of the last vertex remains
    zdata->g.vtxDataPtrs[count]=pos;
    zdata->g.vtxData.resize(pos); // shrink vtxData to remove unused memory

    // 3. update num of vertices after removing exported ones
    zdata->g.numMyVertices = count;

    // 4. count the number of vertices to be imported into myRank only for
    // preallocation
    for (size_t i=0; i<num_import; i++) {
        if (import_procs[i] != zdata->myRank)
            count++;
    }

    // 5. Preallocate vertex data ptrs, vtxGID, and myParts for receiving
    // imported data
    // for vtxData, we have to resize for each vertex, since we don't know
    // the size of the data right now
    zdata->g.vtxDataPtrs.resize(count+1);
    zdata->vtxGID.resize(count);
    zdata->g.myParts.resize(count);
}

// Append the data for a vertex imported from other procs
static void append_data_imported_g(
  void *data,
  int num_gid_entries,
  ZOLTAN_ID_PTR global_id,
  int size,
  char *buf,
  int *ierr
)
{
    ZGRAPHDATA *zdata = (ZGRAPHDATA *)data;

    size_t localvid = zdata->g.numMyVertices;
    zdata->g.numMyVertices++;

    zdata->vtxGID[localvid] = *global_id;
    zdata->g.myParts[localvid] = -1; // as the part id is not sent (I guess the reason is to avoid unnecessary communication)
    size_t start_pos = zdata->g.vtxData.size();
    // the starting position of the current vertex should be set already
    if (start_pos != zdata->g.vtxDataPtrs[localvid]) {
        printf("start_pos != zdata->g.vtxDataPtrs[localvid]\n");
        *ierr = ZOLTAN_FATAL;
    }
    zdata->g.vtxData.resize(start_pos + size/sizeof(read_t));
    for (int i=0; i<size/sizeof(read_t); i++)
        zdata->g.vtxData[start_pos + i] = ((read_t *)buf)[i];
    // set the ending position of current vertex (and the start position
    // of the next vertex)
    zdata->g.vtxDataPtrs[localvid+1] = start_pos + size/sizeof(read_t);
}

// Step 1,2
// Assigne part id for vertices imported and appended in the vertex array,
// but with part id as -1 (as the part ids are unknown when recieving the
// vertex global id and vertex data)
// Step 3,4,5
// Sort vertices by part id, involving updating vtxDataPtrs, vtxGID, vtxData
static void sort_byparts_g(
  void *data,
  int num_gid_entries,
  int num_lid_entries,
  int num_import,
  ZOLTAN_ID_PTR import_global_ids,
  ZOLTAN_ID_PTR import_local_ids,
  int *import_procs,
  int *import_to_part,
  int num_export,
  ZOLTAN_ID_PTR export_global_ids,
  ZOLTAN_ID_PTR export_local_ids,
  int *export_procs,
  int *export_to_part,
  int *ierr
)
{
    ZGRAPHDATA *zdata = (ZGRAPHDATA *)data;
    *ierr = ZOLTAN_OK;

    // 1. build map between global vertex id to the part assign
    size_t num_v_imported = 0;
    std::map<ZOLTAN_ID_TYPE, int> m;
    for (size_t i=0; i<num_import; i++) {
        if (import_procs[i] != zdata->myRank) {
            num_v_imported++;
            m.insert(std::pair<ZOLTAN_ID_TYPE, int>(import_global_ids[i], import_to_part[i]));
        }
    }
    if (num_v_imported != m.size()) {
        printf("num_v_imported != m.size()\n");
        *ierr = ZOLTAN_FATAL;
    }

    // 2. assign part ids for vertices imported and without part id
    for(size_t pos = zdata->g.numMyVertices - num_v_imported; pos < zdata->g.numMyVertices; pos++) {
        if (zdata->g.myParts[pos] != -1) {
            printf("zdata->g.myParts[pos] != -1\n");
            *ierr = ZOLTAN_FATAL;
        }
        if (m.find(zdata->vtxGID[pos]) != m.end())
            zdata->g.myParts[pos] = m[zdata->vtxGID[pos]];
        else {
            printf("(m.find(zdata->vtxGID[pos]) != m.end())\n");
            *ierr = ZOLTAN_FATAL;
        }
    }
    m.clear();

    // 3. generate maps between `part id` and its parts_sizes: `number of vertices` and `total vertices data size`
    std::map<int, std::pair<size_t, size_t>> parts_sizes;
    for(size_t i=0; i<zdata->g.numMyVertices; i++) {
        if (parts_sizes.find(zdata->g.myParts[i]) == parts_sizes.end()) {
            parts_sizes.insert(
                    std::pair<int, std::pair<size_t, size_t>>(
                        zdata->g.myParts[i], std::pair<size_t, size_t>(1, zdata->g.getVtxDataSize(i))));
        } else {
            parts_sizes[zdata->g.myParts[i]].first += 1;
            parts_sizes[zdata->g.myParts[i]].second += zdata->g.getVtxDataSize(i);
        }
    }

    // 4. update sizes to ptr by prefix-summing them
    std::pair<size_t, size_t> p(0, 0);
    for (const auto &[key, value] : parts_sizes) {
        p.first += value.first;
        p.second += value.second;
        parts_sizes[key].first = p.first - parts_sizes[key].first;
        parts_sizes[key].second = p.second - parts_sizes[key].second;
    }

    // 5. sort vtxDataPtrs, vtxGID, vtxData by part ids
    //    - first back up the three arrays
    //    - then read from the back ups, and write to the zdata based on
    //    the parts_sizes
    std::vector<ZOLTAN_ID_TYPE> vtxGID_copy(zdata->vtxGID);
    std::vector<size_t> vtxDataPtrs_copy(zdata->g.vtxDataPtrs);
    std::vector<read_t> vtxData_copy(zdata->g.vtxData);
#ifdef DEBUG
    std::vector<int> myParts_copy(zdata->g.myParts);
#endif

    for (size_t i=0; i<zdata->g.numMyVertices; i++) {
#ifdef DEBUG
        int mypart=myParts_copy[i];
#else
        int mypart=zdata->g.myParts[i];
#endif
        size_t newLocalID=parts_sizes[mypart].first;
        size_t newDataStart=parts_sizes[mypart].second;
        parts_sizes[mypart].first += 1;
        parts_sizes[mypart].second += vtxDataPtrs_copy[i+1]-vtxDataPtrs_copy[i];

        zdata->vtxGID[newLocalID] = vtxGID_copy[i];
        // the last element of vtxDataPtrs should stay the same and does not need to be assigned
        zdata->g.vtxDataPtrs[newLocalID] = newDataStart;
#ifdef DEBUG
        zdata->g.myParts[newLocalID] = myParts_copy[i];
#endif
        if (newDataStart != vtxDataPtrs_copy[i]) {
            for (size_t t=newDataStart, f=vtxDataPtrs_copy[i]; f<vtxDataPtrs_copy[i+1]; t++, f++)
                zdata->g.vtxData[t] = vtxData_copy[f];
        }
    }
    vtxGID_copy.clear();
    vtxDataPtrs_copy.clear();
    vtxData_copy.clear();

#ifdef DEBUG
    myParts_copy.clear();
    int previous_part_id = -1;
    std::vector<size_t> partPtr_new;
    std::vector<int> myParts_new;
    for (size_t i=0; i<zdata->g.numMyVertices; i++) {
        if (zdata->g.myParts[i] != previous_part_id) {
            myParts_new.push_back(zdata->g.myParts[i]);
            partPtr_new.push_back(i);
        }
        previous_part_id = zdata->g.myParts[i];
    }
    partPtr_new.push_back(zdata->g.numMyVertices);
#endif

    // rearrange myParts and malloc partPtr
    zdata->g.partPtr.resize(parts_sizes.size()+1);
    size_t count=0;
    zdata->g.partPtr[count] = 0;
    zdata->g.myParts.resize(parts_sizes.size());
    for (const auto &[key, value] : parts_sizes) {
        zdata->g.myParts[count] = key;
        zdata->g.partPtr[count+1] = value.first;
        count++;
    }
#ifdef DEBUG
    if (partPtr_new.size() != zdata->g.partPtr.size() || myParts_new.size() != zdata->g.myParts.size()) {
        printf("-- partPtr_new.size() {%ld} != zdata->g.partPtr.size() {%ld} || myParts_new.size() {%ld} != zdata->g.myParts.size() {%ld}\n", partPtr_new.size(), zdata->g.partPtr.size(), myParts_new.size(), zdata->g.myParts.size());
        *ierr = ZOLTAN_FATAL;
    } else {
        for (size_t i=0; i<myParts_new.size(); i++) {
            if (myParts_new[i] != zdata->g.myParts[i]) {
                printf("myParts_new[%ld] {%d} != zdata->g.myParts[%ld] {%d}\n", i, myParts_new[i], i, zdata->g.myParts[i]);
                *ierr = ZOLTAN_FATAL;
            }
        }
        for (size_t i=0; i<partPtr_new.size(); i++) {
            if (partPtr_new[i] != zdata->g.partPtr[i]) {
                printf("partPtr_new[%ld] {%d} != zdata->g.partPtr[%ld] {%d}\n", i, partPtr_new[i], i, zdata->g.partPtr[i]);
                *ierr = ZOLTAN_FATAL;
            }
        }
    }
    partPtr_new.clear();
    myParts_new.clear();
#endif
}

void print_graph_read_gid (ZGRAPHDATA & zdata)
{
    size_t numMyParts = zdata.g.getNumMyParts();
    for (size_t partID = 0; partID < numMyParts; partID++) {
        std::string output_filename(outputDirName + "/vtxGID_r" + std::to_string(rank) + "_p" + std::to_string(zdata.g.myParts[partID]) + ".fa");
        FILE *f = fopen(output_filename.c_str(), "w");
        fprintf(f, "%lu\n", zdata.g.partPtr[partID+1] - zdata.g.partPtr[partID]);

        for (int i=zdata.g.partPtr[partID]; i<zdata.g.partPtr[partID+1]; i++)
            fprintf(f, "%lu\n", zdata.vtxGID[i]);

        fclose(f);
    }
}

/*
 * calling zoltan with ParMetis to partition the vertices
 * Inputs:
 *     argc, argv: inputs for the executable to initalize Zoltan
 *     g: structure stores the data for the graph located in current
 *     mpi proc, including data for the vertices
 *     its vertices and the
 *     parts: the number of partitions
 */
// TODO: reduce code redandency
int perform_zoltan_partition_g (int argc, char **argv,
        GRAPH_DATA & g, const int parts)
{
    int rc;
    float version;
    struct Zoltan_Struct *zz;
    int changes, numGidEntries, numLidEntries;
    ZOLTAN_VERTICES_EXCHANGE in, out;

    // Logging the total number of edges (size of the graph)
    size_t mynumofedges = g.getNumMyEdges();
    size_t totalnumofedges;
    MPI_Reduce(&mynumofedges, &totalnumofedges, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        printf("Total number of edges: %lu\n", totalnumofedges);
    }

#ifdef TIME_PROFILE
    double partition_time, migrate_time;
    double begin_time, end_time;
    MPI_Barrier(MPI_COMM_WORLD);
    begin_time = MPI_Wtime ();
#endif

    /******************************************************************
     ** Initialize Zoltan
     ******************************************************************/
#ifdef DEBUG
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
        printf("init zoltan\n");
#endif
    rc = Zoltan_Initialize(argc, argv, &version);
    if (rc != ZOLTAN_OK){
        fprintf(stderr, "Zoltan partitioner: Failed to initialize");
        MPI_Finalize();
        Zoltan_Destroy(&zz);
        exit(0);
    }

    /******************************************************************
     ** Create a Zoltan library structure for this instance of load
     ** balancing.  Set the parameters and query functions that will
     ** govern the library's calculation.  See the Zoltan User's
     ** Guide for the definition of these and many other parameters.
     ******************************************************************/
#ifdef DEBUG
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
        printf("create zoltan struct\n");
#endif
    zz = Zoltan_Create(MPI_COMM_WORLD);

#ifdef DEBUG
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
        printf("set parameter\n");
#endif
    /* General parameters */
#ifdef DEBUG
    Zoltan_Set_Param(zz, "DEBUG_LEVEL", "2");
#else
    Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
#endif
    //TODO: change based on partitioner input
    Zoltan_Set_Param(zz, "LB_METHOD", "GRAPH");   /* partitioning method */
    Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");
    Zoltan_Set_Param(zz, "GRAPH_PACKAGE", "ParMETIS"); /* version of method */
    Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");/* global IDs are integers */
    Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");/* local IDs are integers */
    Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL"); /* export AND import lists */
    Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "0"); /* use Zoltan default vertex weights */
    Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM", "0");/* use Zoltan default hyperedge weights */
    Zoltan_Set_Param(zz, "IMBALANCE_TOL", "1.01"); /* the load imbalance setting, the upper bound load is 1.01*(averageload) */
    char param_str[10];
    sprintf(param_str,"%d",parts);
    Zoltan_Set_Param(zz, "NUM_GLOBAL_PARTS", param_str);/* partition the vertices to `parts` parts*/
    if (rank == 0) printf("Number of parts given to Zoltan: %s\n", param_str);
    if (partitioning_seed != 0) {
        sprintf(param_str,"%d",partitioning_seed);
        Zoltan_Set_Param(zz, "SEED", param_str); /* the seeds for partitioning */
        if (rank == 0) printf("Set seeds for Zoltan as: %s\n", param_str);
    }

    /* Application defined query functions */
    Zoltan_Set_Num_Obj_Fn(zz, get_number_of_vertices_g, &g);
    Zoltan_Set_Obj_List_Fn(zz, get_vertex_list_g, &g);
    Zoltan_Set_Num_Edges_Multi_Fn(zz, get_edge_list_size_g, &g);
    Zoltan_Set_Edge_List_Multi_Fn(zz, get_edge_list_g, &g);

    /******************************************************************
     ** Zoltan can now partition the vertices of graph.
     ******************************************************************/

#ifdef DEBUG
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
        printf("start lb partitioning\n");
#endif

    rc = Zoltan_LB_Partition(zz, /* input (all remaining fields are output) */
            &changes,        /* 1 if partitioning was changed, 0 otherwise */
            &numGidEntries,  /* Number of integers used for a global ID */
            &numLidEntries,  /* Number of integers used for a local ID */
            &in.num,      /* Number of vertices to be sent to me */
            &in.globalGids,  /* Global IDs of vertices to be sent to me */
            &in.localGids,   /* Local IDs of vertices to be sent to me */
            &in.procs,    /* Process rank for source of each incoming vertex */
            &in.toPart,   /* New partition for each incoming vertex */
            &out.num,      /* Number of vertices I must send to other processes*/
            &out.globalGids,  /* Global IDs of the vertices I must send */
            &out.localGids,   /* Local IDs of the vertices I must send */
            &out.procs,    /* Process to which I send each of the vertices */
            &out.toPart);  /* Partition to which each vertex will belong */

    if (rc != ZOLTAN_OK){
        fprintf(stderr, "Zoltan partitioner: Failed to partition");
        MPI_Finalize();
        Zoltan_Destroy(&zz);
        exit(0);
    }

#ifdef TIME_PROFILE
    end_time = MPI_Wtime ();
    MPI_Barrier(MPI_COMM_WORLD);
    partition_time = end_time - begin_time;
    begin_time = MPI_Wtime ();
#endif

#ifdef DEBUG
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
        printf("start set parameters for migrate\n");
#endif
    /* Pack and Unpack function for anto migration */
    Zoltan_Set_Obj_Size_Fn(zz, get_datasize_of_vertices_g, &g);
    Zoltan_Set_Pack_Obj_Fn(zz, get_data_of_vertices_g, &g);
    ZGRAPHDATA zgraphdata(g, rank);
    Zoltan_Set_Mid_Migrate_PP_Fn(zz, reorganize_vdataparts_g, &zgraphdata);
    Zoltan_Set_Unpack_Obj_Fn(zz, append_data_imported_g, &zgraphdata);
    Zoltan_Set_Post_Migrate_PP_Fn(zz, sort_byparts_g, &zgraphdata);
#ifdef DEBUG
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
        printf("start migrate\n");
#endif

    // Data migration
    rc = Zoltan_Migrate(zz, /* input (all remaining fields are output) */
            in.num,      /* Number of vertices to be sent to me */
            in.globalGids,  /* Global IDs of vertices to be sent to me */
            in.localGids,   /* Local IDs of vertices to be sent to me */
            in.procs,    /* Process rank for source of each incoming vertex */
            in.toPart,   /* New partition for each incoming vertex */
            out.num,      /* Number of vertices I must send to other processes*/
            out.globalGids,  /* Global IDs of the vertices I must send */
            out.localGids,   /* Local IDs of the vertices I must send */
            out.procs,    /* Process to which I send each of the vertices */
            out.toPart);  /* Partition to which each vertex will belong */
    if (rc != ZOLTAN_OK){
        fprintf(stderr, "Zoltan partitioner: Failed to migrate data");
        MPI_Finalize();
        Zoltan_Destroy(&zz);
        exit(0);
    }
#ifdef TIME_PROFILE
    end_time = MPI_Wtime ();
    MPI_Barrier(MPI_COMM_WORLD);
    migrate_time = end_time - begin_time;
#endif

#ifdef DEBUG
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
        printf("start cleanup\n");

    print_graph_read_gid(zgraphdata);
#endif

    Zoltan_Destroy(&zz);

#ifdef TIME_PROFILE
    double global_partition_time_sum, global_migrate_time_sum;
    double global_partition_time_max, global_migrate_time_max;

    MPI_Reduce(&partition_time, &global_partition_time_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&partition_time, &global_partition_time_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&migrate_time, &global_migrate_time_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&migrate_time, &global_migrate_time_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        printf ("Average time for partitioning across all procs (secs):   %f \n",
                (double)global_partition_time_sum/(double)size);
        printf ("Maximum time for partitioning across all procs (secs):   %f \n",
                (double)global_partition_time_max);
        printf ("Average time for data migration across all procs (secs): %f \n",
                (double)global_migrate_time_sum/(double)size);
        printf ("Maximum time for data migration across all procs (secs): %f \n",
                (double)global_migrate_time_max);
    }
#endif

    return 0;
}
