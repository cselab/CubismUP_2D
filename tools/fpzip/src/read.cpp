#include <cstdio>
#include <cstdlib>
#include "pcdecoder.h"
#include "rcqsmodel.h"
#include "front.h"
#include "fpzip.h"
#include "codec.h"
#include "read.h"

#if FPZIP_FP == FPZIP_FP_FAST || FPZIP_FP == FPZIP_FP_SAFE
// decompress 3D array at specified precision using floating-point arithmetic
template <typename T, unsigned bits>
static void
decompress3d(
  RCdecoder* rd,   // entropy decoder
  T*         data, // flattened 3D array to decompress to
  unsigned   nx,   // number of x samples
  unsigned   ny,   // number of y samples
  unsigned   nz    // number of z samples
)
{
  // initialize decompressor
  typedef PCmap<T, bits> MAP;
  RCmodel* rm = new RCqsmodel(false, PCdecoder<T, MAP>::symbols);
  PCdecoder<T, MAP>* fd = new PCdecoder<T, MAP>(rd, &rm);
  FRONT<T> f(nx, ny);

  // decode difference between predicted (p) and actual (a) value
  unsigned x, y, z;
  for (z = 0, f.advance(0, 0, 1); z < nz; z++)
    for (y = 0, f.advance(0, 1, 0); y < ny; y++)
      for (x = 0, f.advance(1, 0, 0); x < nx; x++) {
        #if FPZIP_FP == FPZIP_FP_SAFE
        volatile T p = f(1, 1, 1);
        p += f(1, 0, 0);
        p -= f(0, 1, 1);
        p += f(0, 1, 0);
        p -= f(1, 0, 1);
        p += f(0, 0, 1);
        p -= f(1, 1, 0);
        #else
        T p = f(1, 0, 0) - f(0, 1, 1) +
              f(0, 1, 0) - f(1, 0, 1) +
              f(0, 0, 1) - f(1, 1, 0) +
              f(1, 1, 1);
        #endif
        T a = fd->decode(p);
        *data++ = a;
        f.push(a);
      }

  delete fd;
  delete rm;
}
#elif FPZIP_FP == FPZIP_FP_EMUL
#include "fpe.h"
// decompress 3D array at specified precision using floating-point emulation
template <typename T, unsigned bits>
static void
decompress3d(
  RCdecoder* rd,   // entropy encoder
  T*         data, // flattened 3D array to decompress to
  unsigned   nx,   // number of x samples
  unsigned   ny,   // number of y samples
  unsigned   nz    // number of z samples
)
{
  // initialize decompressor
  typedef PCmap<T, bits> MAP;
  typedef FPE<T> FLOAT;
  RCmodel* rm = new RCqsmodel(false, PCdecoder<T, MAP>::symbols);
  PCdecoder<T, MAP>* fd = new PCdecoder<T, MAP>(rd, &rm);
  FRONT<FLOAT> f(nx, ny);

  // decode difference between predicted (p) and actual (a) value
  unsigned x, y, z;
  for (z = 0, f.advance(0, 0, 1); z < nz; z++)
    for (y = 0, f.advance(0, 1, 0); y < ny; y++)
      for (x = 0, f.advance(1, 0, 0); x < nx; x++) {
        FLOAT p = f(1, 0, 0) - f(0, 1, 1) +
                  f(0, 1, 0) - f(1, 0, 1) +
                  f(0, 0, 1) - f(1, 1, 0) +
                  f(1, 1, 1);
        T a = fd->decode(T(p));
        *data++ = a;
        f.push(a);
      }
                                                                                
  delete fd;
  delete rm;
}
#else // FPZIP_FP_INT
// decompress 3D array at specified precision using integer arithmetic
template <typename T, unsigned bits>
static void
decompress3d(
  RCdecoder* rd,   // entropy decoder
  T*         data, // flattened 3D array to decompress to
  unsigned   nx,   // number of x samples
  unsigned   ny,   // number of y samples
  unsigned   nz    // number of z samples
)
{
  // initialize decompressor
  typedef PCmap<T, bits> TMAP;
  typedef typename TMAP::RANGE U;
  typedef PCmap<U, bits, U> UMAP;
  RCmodel* rm = new RCqsmodel(false, PCdecoder<U, UMAP>::symbols);
  PCdecoder<U, UMAP>* fd = new PCdecoder<U, UMAP>(rd, &rm);
  TMAP map;
  FRONT<U> f(nx, ny, map.forward(0));

  // decode difference between predicted (p) and actual (a) value
  unsigned x, y, z;
  for (z = 0, f.advance(0, 0, 1); z < nz; z++)
    for (y = 0, f.advance(0, 1, 0); y < ny; y++)
      for (x = 0, f.advance(1, 0, 0); x < nx; x++) {
        U p = f(1, 0, 0) - f(0, 1, 1) +
              f(0, 1, 0) - f(1, 0, 1) +
              f(0, 0, 1) - f(1, 1, 0) +
              f(1, 1, 1);
        U a = fd->decode(p);
        *data++ = map.inverse(a);
        f.push(a);
      }

  delete fd;
  delete rm;
}
#endif

// decompress 4D array
template <typename T>
static void
decompress4d(
  RCdecoder* rd,   // entropy decoder
  T*         data, // flattened 4D array to decompress to
  int*       prec, // per field precision
  unsigned   nx,   // number of x samples
  unsigned   ny,   // number of y samples
  unsigned   nz,   // number of z samples
  unsigned   nf    // number of fields
)
{
  // decompress one field at a time
  for (unsigned i = 0; i < nf; i++) {
    int bits = rd->template decode <unsigned>(32);
    if (prec)
      prec[i] = bits;
    switch (bits) {
      case subsize(T, 1): //  8-bit float, 16-bit double
        decompress3d<T, subsize(T, 1)>(rd, data, nx, ny, nz);
        break;
      case subsize(T, 2): // 16-bit float, 32-bit double
        decompress3d<T, subsize(T, 2)>(rd, data, nx, ny, nz);
        break;
      case subsize(T, 3): // 24-bit float, 48-bit double
        decompress3d<T, subsize(T, 3)>(rd, data, nx, ny, nz);
        break;
      case subsize(T, 4): // 32-bit float, 64-bit double
        decompress3d<T, subsize(T, 4)>(rd, data, nx, ny, nz);
        break;
      default:
        fprintf(stderr, "fpzip: invalid precision %d in file\n", bits);
        abort();
        break;
    }
    data += nx * ny * nz;
  }
}

void
read_header(
  RCdecoder* rd,
  unsigned   nx,
  unsigned   ny,
  unsigned   nz,
  unsigned   nf,
  int        dp
)
{
  // magic
  if (rd->decode<unsigned>(8) != 'f' ||
      rd->decode<unsigned>(8) != 'p' ||
      rd->decode<unsigned>(8) != 'z' ||
      rd->decode<unsigned>(8) != '\0') {
    fprintf(stderr, "fpzip: not an fpz stream\n");
    abort();
  }

  // format version
  if (rd->decode<unsigned>(16) != FPZ_MAJ_VERSION ||
      rd->decode<unsigned>(16) != FPZ_MIN_VERSION) {
    fprintf(stderr, "fpzip: format version not supported\n");
    abort();
  }

  // array dimensions
  if (rd->decode<unsigned>(32) != nf ||
      rd->decode<unsigned>(32) != nz ||
      rd->decode<unsigned>(32) != ny ||
      rd->decode<unsigned>(32) != nx) {
    fprintf(stderr, "fpzip: array dimensions do not match\n");
    fprintf(stderr, "fpzip: array dimensions do not match\n");
    abort();
  }

  // single or double precision
  if (rd->decode() != !!dp) {
    fprintf(stderr, "fpzip: floating-point type does not match\n");
    abort();
  }
}

void
get_header(
  RCdecoder* rd,
  unsigned   *nx,
  unsigned   *ny,
  unsigned   *nz,
  unsigned   *nf,
  int        *dp
)
{
  // magic
  if (rd->decode<unsigned>(8) != 'f' ||
      rd->decode<unsigned>(8) != 'p' ||
      rd->decode<unsigned>(8) != 'z' ||
      rd->decode<unsigned>(8) != '\0') {
    fprintf(stderr, "fpzip: not an fpz stream\n");
    abort();
  }

  // format version
  if (rd->decode<unsigned>(16) != FPZ_MAJ_VERSION ||
      rd->decode<unsigned>(16) != FPZ_MIN_VERSION) {
    fprintf(stderr, "fpzip: format version not supported\n");
    abort();
  }

  // array dimensions
  *nf = rd->decode<unsigned>(32);
  *nz = rd->decode<unsigned>(32);
  *ny = rd->decode<unsigned>(32);
  *nx = rd->decode<unsigned>(32);

  *dp = rd->decode();
}

static void
fpzip_stream_read(
  RCdecoder* rd,   // entropy decoder
  void*      data, // array to read
  int*       prec, // per field bits of precision
  int        dp,   // double precision array if nonzero
  unsigned   nx,   // number of x samples
  unsigned   ny,   // number of y samples
  unsigned   nz,   // number of z samples
  unsigned   nf    // number of fields
)
{
  rd->init();
  read_header(rd, nx, ny, nz, nf, dp);
  if (dp)
    decompress4d(rd, (double*)data, prec, nx, ny, nz, nf);
  else
    decompress4d(rd, (float*)data, prec, nx, ny, nz, nf);
}

static void
fpzip_stream_get_header(
  RCdecoder* rd,   // entropy decoder
  int        *dp,   // double precision array if nonzero
  unsigned   *nx,   // number of x samples
  unsigned   *ny,   // number of y samples
  unsigned   *nz,   // number of z samples
  unsigned   *nf    // number of fields
)
{
  rd->init();
  get_header(rd, nx, ny, nz, nf, dp);
}

// read and decompress a single or double precision 4D array from file
unsigned
fpzip_file_read(
  FILE*       file, // binary input stream
  void*       data, // array to read
  int*        prec, // per field bits of precision
  int         dp,   // double precision array if nonzero
  unsigned    nx,   // number of x samples
  unsigned    ny,   // number of y samples
  unsigned    nz,   // number of z samples
  unsigned    nf    // number of fields
)
{
  RCfiledecoder* rd = new RCfiledecoder(file);
  fpzip_stream_read(rd, data, prec, dp, nx, ny, nz, nf);
  unsigned bytes = rd->error ? 0 : rd->bytes();
  delete rd;
  return bytes;
}

// read and decompress a single or double precision 4D array from file
unsigned
fpzip_memory_read(
  const void* buffer, // pointer to compressed data
  void*       data,   // array to read
  int*        prec,   // per field bits of precision
  int         dp,     // double precision array if nonzero
  unsigned    nx,     // number of x samples
  unsigned    ny,     // number of y samples
  unsigned    nz,     // number of z samples
  unsigned    nf      // number of fields
)
{
  RCmemdecoder* rd = new RCmemdecoder(buffer);
  fpzip_stream_read(rd, data, prec, dp, nx, ny, nz, nf);
  unsigned bytes = rd->error ? 0 : rd->bytes();
  delete rd;
  return bytes;
}

// read and decompress a single or double precision 4D array from file
void
fpzip_memory_get_header(
  const void* buffer, // pointer to compressed data
  int         *dp,     // double precision array if nonzero
  unsigned    *nx,     // number of x samples
  unsigned    *ny,     // number of y samples
  unsigned    *nz,     // number of z samples
  unsigned    *nf      // number of fields
)
{
  RCmemdecoder* rd = new RCmemdecoder(buffer);
  fpzip_stream_get_header(rd, dp, nx, ny, nz, nf);
  delete rd;
}

// wrappers for fortran calls
void
fpzip_file_read_f(
  const char* path, // path to input file
  void*       data, // array to read
  int*        prec, // per field bits of precision
  const int*  dp,   // double precision array if nonzero
  const int*  nx,   // number of x samples
  const int*  ny,   // number of y samples
  const int*  nz,   // number of z samples
  const int*  nf    // number of fields
)
{
  FILE* file = fopen(path, "rb");
  if (!file) {
    fprintf(stderr, "fpzip: cannot open file %s\n", path);
    abort();
  }
  if (!fpzip_file_read(file, data, prec, *dp, *nx, *ny, *nz, *nf)) {
    fprintf(stderr, "fpzip: cannot read file %s\n", path);
    abort();
  }
  fclose(file);
}

void
fpzip_file_read_f_(
  const char* path, // path to input file
  void*       data, // array to read
  int*        prec, // per field bits of precision
  const int*  dp,   // double precision array if nonzero
  const int*  nx,   // number of x samples
  const int*  ny,   // number of y samples
  const int*  nz,   // number of z samples
  const int*  nf    // number of fields
)
{
  fpzip_file_read_f(path, data, prec, dp, nx, ny, nz, nf);
}
