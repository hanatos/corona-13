#pragma once

#include <stdint.h>
#include <stddef.h>
#include "vol/payload.h"

#define VOL_HEADER_MAGIC 0x9bae454d
#define VOL_HEADER_VERSION 8

typedef enum vol_channel_t
{
  s_vol_density=1,
  s_vol_temperature=2,
}
vol_channel_t;

typedef union vol_index_t
{
  uint32_t idx;
  struct
  {
    uint32_t i : 3;
    uint32_t j : 3;
    uint32_t k : 3;
    uint32_t l : 23;
  };
}
vol_index_t;


// ==============================================================
// the file format is as follows:
// 
// [vol_header_t] [vol_payload_t *] [vol_node_t *]
//
// because the tree structure (vol_node_t *) is kept in memory
// during construction and dumped last, when the construction
// is complete. payload data on the other hand is dumped as it
// appears, to facilitate almost out of core construction.
//
// during run-time (reading), this structure is just mmapped
// and the dynamic vol_tree_t struct will contain individual
// pointers and bounding boxes and transforms.
// ==============================================================

// node in the tree structure, does not contain any payload data.
// fixes the hierarchical grid with branching factor 512 (8*8*8).
// must be uint64_t aligned.
typedef struct vol_node_t
{
  // these could be stored in fewer bits, as the payload is
  // 4k-page aligned.
  uint64_t data_static0 : 1;  // mark all data in this block as static
  uint64_t data_offset0 : 63; // byte offset for off[0..255]
  uint64_t data_static1 : 1;
  uint64_t data_offset1 : 63; // for off[256..511]
  uint32_t off255_empty : 1;
  uint32_t node_offset0 : 31;
  uint32_t off511_empty : 1;
  uint32_t node_leaf    : 1;
  uint32_t node_offset1 : 30;
  uint32_t lh_offset0;
  uint32_t lh_offset1;
  uint8_t off[512];      // child node i is base_offsetX + off[i]*sizeof(*)
}
vol_node_t;

struct vol_payload_flux_t;
// 4096 bytes, pad to page boundaries
typedef struct vol_header_t
{
  uint32_t magic;       // file id magic number
  uint32_t version;     // file format version
  uint64_t nodes;       // byte offset in file to nodes
  float aabb[6];        // initial bounding box (can be relocated in vol_tree_t)
  float content_box[6]; // max extents of content inside this
  float voxel_size;     // side length of one voxel in wold space corresponding to aabb
  float rot[3];         // euler rotation angles (XYZ order)
  float loc[3];         // world space offset after scaling with voxel_size and rotation
  int32_t depth;        // depth of tree (1 means just root node)
  uint64_t light;       // offset to light hierarchy
  int32_t isstatic;     // whether root node payload is static
  int32_t shaderid;     // globally fixed shader id
  uint64_t end;         // byte offset to end of file
  uint8_t pad[3976];    // make sure payload will start at page boundaries
}
vol_header_t;

// this is a dynamic struct linking together the mmapped file for rendering use.
// allocated in heap memory.
typedef struct vol_tree_t
{
  float aabb[6];                // relocated aabb of full tree
  float content_box[6];         // max extents of content inside this
  float voxel_size;             // side length of one voxel in wold space corresponding to aabb
  int   transform;              // set if the matrix below is valid and needed:
  float w2o[3][4];              // world space to object space (rotation + offset) matrix

  int fd;                       // mmap stuff
  uint8_t *data;
  size_t data_size;

  vol_node_t *nodes;            // convenience pointer to start of tree structure (header + header->nodes)
  uint8_t *payload;             // convenience pointer to start of data (header + sizeof(header))
  struct vol_payload_flux_t *light; // convenience pointer to start of light hierarchy nodes

  vol_node_t *root_node;        // convenience pointers to root node
  vol_payload_t *root_payload;
  struct vol_payload_flux_t *root_light;

  vol_emission_shader_t shader; // pointer to shader used during building

  vol_header_t *header;         // alias to *data, header from file.
}
vol_tree_t;
