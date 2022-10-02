//
// Created by Ashton Hess on 9/20/22.
//

#include "FileIO.h"

FileIO::FileIO(){}

vector<vector<Node*>>FileIO::read_2d(char*argv[], int*nc_lines){
    double x,y;
    int N = stoi(argv[2]);
    vector<vector<Node*>>rVec;
    if (FILE *fp = fopen(argv[1], "r")) {
        for(int i = 0; i < N*7; i++){
            fscanf(fp, "%lf,%lf\n", &x, &y);
            //auto*data=(Node*)malloc(sizeof(Node);
            Node *data = new Node();
            data->x=x;
            data->y=y;
            data->col = N%i;
            data->generation=N/i;
            rVec[N%i][N/i]= data;
        }
//        while (fscanf(fp, "%lf,%lf\n", &x, &y) == 2) {
//            auto*data=(Node*)malloc(sizeof(Node));
//            data->x=x;
//            data->y=y;
//            data->N/;
//            data->generation=gen;
////            rVec.push_back(*data);
//            rVec[col][gen]=data;
//            *nc_lines+=1;
//            col++;
//            if(col>=stoi(argv[2])){
//                col=0;
//                gen++;
//            }
//        }
        fclose(fp);
        return rVec;
    }else{
        cout<<"Error opening "<<argv[1]<<endl;
        return vector<vector<Node*>>();
    }
}

vector<Node*>FileIO::read_nc(char*argv[], int*nc_lines){
    double x,y;
    vector<Node*>rVec;
    if (FILE *fp = fopen(argv[1], "r")) {
        int col=0;
        int gen=0;
        while (fscanf(fp, "%lf,%lf\n", &x, &y) == 2) {
            auto*data=(Node*)malloc(sizeof(Node));
            data->x=x;
            data->y=y;
            data->col=col;
            data->generation=gen;
            if(gen!=0){
                data->theta = atan2(rVec[(gen-1)*stoi(argv[2])+col]->y - y, rVec[(gen-1)*stoi(argv[2])+col]->x -x)*180/M_PI;
            }
            rVec.push_back(data);
            *nc_lines+=1;
            col++;
            if(col>=stoi(argv[2])){
                col=0;
                gen++;
            }
        }
        fclose(fp);
        return rVec;
    }else{
        cout<<"Error opening "<<argv[1]<<endl;
        return vector<Node*>();
    }
}

Node**FileIO::read_2d_array(char*argv[]){
    int N = stoi(argv[2]);
    int generations=stoi(argv[3]);
    Node**tree;
    tree = (Node**)malloc(sizeof(Node*)*N);
    for (int i = 0; i < N; ++i) {
        tree[i]=(Node*)malloc(sizeof(Node)*generations+7);//HARDCODING A 7 HERE FOR GIVEN GENERATIONS... MAY BE UNSAFE...
//        tree[i]= new Node();
    }

    double x,y;

    if (FILE *fp = fopen(argv[1], "r")) {
        int col=0;
        int gen=0;
        int counter=0;
        while (fscanf(fp, "%lf,%lf\n", &x, &y) == 2) {
            col = counter%N;
            if(counter%N==0 && counter!=0){
                gen++;
            }
            tree[col][gen].col=col;
            tree[col][gen].generation=gen;
            tree[col][gen].x=x;
            tree[col][gen].y=y;
            if(gen!=0) {
                tree[col][gen].theta = atan2(tree[col][gen - 1].y - y, tree[col][gen - 1].x - x) * 180 / M_PI;
            }
//            data->theta = atan2(rVec[(gen-1)*stoi(argv[2])+col]->y - y, rVec[(gen-1)*stoi(argv[2])+col]->x -x)*180/M_PI;
            counter++;
//            cout<<counter<<endl;
        }
        fclose(fp);
        return tree;
    }else{
        cout<<"Error opening "<<argv[1]<<endl;
        return NULL;
    }

}
