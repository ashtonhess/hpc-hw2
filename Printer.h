//
// Created by Ashton Hess on 10/1/22.
//

#ifndef HW2_PRINTER_H
#define HW2_PRINTER_H

#include <iostream>
#include <vector>
#include "DataStructures.h"
using namespace std;

class Printer {
private:
public:
    void print_tree_1d(vector<Node*>tree);
    void print_tree_2d(vector<vector<Node*> >tree);
};


#endif //HW2_PRINTER_H
