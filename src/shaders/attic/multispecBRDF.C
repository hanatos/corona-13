// ==========================================================================
// $Id$
// ==========================================================================
// (C)opyright:
//
//   Max-Planck-Institut fuer Informatik
//   Stuhlsatzenhausweg 85, 66123 Saarbruecken
//   http://www.mpi-sb.mpg.de
// 
// Creator: hullin (Matthias Hullin)
// Email:   hullin@mpi-sb.mpg.de
// ==========================================================================
// $Log$
// ==========================================================================

#ifndef MULTISPECBRDF_C
#define MULTISPECBRDF_C

#include "multispecBRDF.hh"

// read from common data block (save 2x memory for split images)
multispecBRDF::multispecBRDF(char *buf)
{
  int *header = (int *)buf;
  _thetaBins = header[0];
  assert(header[1] == header[0]); // at some point i thought it
  // wouldn't make much sense if
  // the theta binning does not
  // agree

  _phiBins = header[2];

  _illumFrom = header[11];
  _illumTo = header[12];
  _camFrom = header[13];
  _camTo = header[14];
  _illumStep = header[15];
  _camStep = header[15];

  _illumBins = (_illumTo - _illumFrom) / _illumStep + 1;
  _camBins = (_camTo - _camFrom) / _camStep + 1;

  // assert(header[3] == _illumBins);
  // assert(header[4] == _camBins);


  // precompute directions to nearest neighbors along the five axes
  dir[4] = 1;
  dir[3] = _camBins;
  dir[2] = dir[3] * _illumBins;
  dir[1] = dir[2] * _phiBins;
  dir[0] = dir[1] * _thetaBins;

  for (int i=0; i<32; ++i) {
    size_t b0 = (i & 16) >> 4;
    size_t b1 = (i & 8) >> 3;
    size_t b2 = (i & 4) >> 2;
    size_t b3 = (i & 2) >> 1;
    size_t b4 = (i & 1);
    dircube[i] = b0 * dir[0] + b1 * dir[1] + b2 * dir[2] + b3 * dir[3] + b4 * dir[4];
  }
  data = (float*)(buf + 64);
}

multispecBRDF::multispecBRDF(const string &filename)
{
  int header[16];
  ifstream file (filename.c_str(), ios::in|ios::binary);
  assert(file);
  file.read((char*)&header, 64);
  /* Multispectral BRDF file header
     see also: http://people.csail.mit.edu/addy/research/brdf/

#  Name  	Description                                             here:
0  dim[0] 	Number of bins for the 0th dimension                    [[theta_in]] // _thetaBins
1  dim[1] 	Number of bins for the 1st dimension                    [[theta_out]] // _thetaBins
2  dim[2] 	Number of bins for the 2nd dimension                    [[phi_diff]] // _phiBins
3  dim[3] 	Number of bins for the 3rd dimension                    [[_illumWLCount]]
4  reserved 	                                                        [[_camWLCount]]
5  reserved 	 
6  paramType 	0 - reserved, 1 - Standard Parameterization             [[1]]
7  binType 	0 - Uniform binning, 1 - reserved                       [[0]]
8  reserved 	 
9  half_data 	0 - phi_diff range [0..2pi], 1 - phi_diff range [0..pi] [[1]]
10 num_channels always 3 (RGB)                                          [[_illumWLCount * _camWLCount]]
11 reserved 	                                                        [[_mIllumFrom]]
12 reserved 	                                                        [[_mIllumTo]]
13 reserved 	                                                        [[_mCamFrom]]
14 reserved 	                                                        [[_mCamTo]]
15 reserved 	                                                        [[_mIllumStep]]
*/      

  _thetaBins = header[0];
  assert(header[1] == header[0]); // at some point i thought it
  // wouldn't make much sense if
  // the theta binning does not
  // agree

  _phiBins = header[2];

  _illumFrom = header[11];
  _illumTo = header[12];
  _camFrom = header[13];
  _camTo = header[14];
  _illumStep = header[15];
  _camStep = header[15];

  _illumBins = (_illumTo - _illumFrom) / _illumStep + 1;
  _camBins = (_camTo - _camFrom) / _camStep + 1;

  assert(header[3] == _illumBins);
  assert(header[4] == _camBins);


  // precompute directions to nearest neighbors along the five axes
  dir[4] = 1;
  dir[3] = _camBins;
  dir[2] = dir[3] * _illumBins;
  dir[1] = dir[2] * _phiBins;
  dir[0] = dir[1] * _thetaBins;

  for (int i=0; i<32; ++i) {
    size_t b0 = (i & 16) >> 4;
    size_t b1 = (i & 8) >> 3;
    size_t b2 = (i & 4) >> 2;
    size_t b3 = (i & 2) >> 1;
    size_t b4 = (i & 1);
    dircube[i] = b0 * dir[0] + b1 * dir[1] + b2 * dir[2] + b3 * dir[3] + b4 * dir[4];
  }

  size_t datasize = dir[0] * _thetaBins * sizeof(float);
  data = (float*)new char[datasize];
  file.read((char*)data, datasize);
  file.close();

};

#endif /* MULTISPECBRDF_C */







