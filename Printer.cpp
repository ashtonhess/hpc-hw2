//
// Created by Ashton Hess on 10/1/22.
//

#include "Printer.h"

void Printer::print_tree_1d(vector<Node*>tree){
    for (long unsigned int i = 0; i < tree.size(); ++i) {
        cout<<"x: "<<tree.at(i)->x;
        cout<<"\ty: "<<tree.at(i)->y;
        cout<<"\tcol: "<<tree.at(i)->col;
        cout<<"\tgen: "<<tree.at(i)->generation;
        cout<<"\ttheta: "<<tree.at(i)->theta<<endl;
    }
}

void Printer::print_tree_2d(vector<vector<Node*> >tree){
    int cols = tree.size();
    int rows = tree.at(0).size();
    cout<<"print_tree_2d():\tCols: "<<cols<<"\tRows: "<<rows<<endl;
}