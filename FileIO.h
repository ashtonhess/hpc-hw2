//
// Created by Ashton Hess on 9/20/22.
//

#ifndef HPC22_AAHB8F_FILEIO_H
#define HPC22_AAHB8F_FILEIO_H

#include <iostream>
using namespace std;
#include <stdlib.h>
#include <vector>
#include "DataStructures.h"
#include <math.h>

class FileIO {
private:
public:
    FileIO();
//    vector<Node>*readNodes(char*filename);
    vector<Node*>read_nc(char*argv[], int*nc_lines);
    vector<vector<Node*>>read_2d(char*argv[], int*nc_lines);
    Node**read_2d_array(char*argv[]);

    };


#endif //HPC22_AAHB8F_FILEIO_H
