#pragma once
#include <ctime>

namespace coulomb {

// Global variables
inline constexpr double pi = 3.1415926535897932;
inline int K_SAVE_TIME = 0;                  // record the file numbers
inline bool FLAG_FILENAME_WITH_NUM = false;  // append number to the filename

inline bool FLAG_SAVEFLUX = false;
inline int NUM_MOVED = 0;
inline int NUM_RESAMPLE = 0;

inline double resample_spatial_ratio = 1.0;

inline double SYNC_TIME = 0.0;

// global variables on time
inline std::clock_t t0_all, t1_all, t0_coll, t1_coll, t0_adve, t1_adve,
    t0_resamp, t1_resamp;

}  // namespace coulomb

namespace coulomb2 {

// Global variables
inline double pi = 3.1415926535897932;
inline int FLAG_PRECOMPUTE_ALPHA_U = 1;
inline int K_SAVE_TIME = 0;                  // record the file numbers
inline bool FLAG_FILENAME_WITH_NUM = false;  // append number to the filename

inline bool FLAG_SAVEFLUX = false;
inline int NUM_MOVED = 0;
inline int NUM_RESAMPLE = 0;

inline double resample_spatial_ratio = 1.0;

inline double SYNC_TIME = 0.0;

// global variables on time
inline std::clock_t t0_all, t1_all, t0_coll, t1_coll, t0_adve, t1_adve,
    t0_resamp, t1_resamp;

}  // namespace coulomb2