#include <stdlib.h>
#include <unistd.h>
#include <iostream>
using namespace std;
#include <math.h>
#include "DataStructures.h"
#include "Printer.h"
#include <vector>
#include "FileIO.h"
#include <chrono>
#include <algorithm>


#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#include <sys/shm.h>
#include <unistd.h>
#include <cstdlib>

std::vector<std::vector<Node*>> convert_2d(std::vector<Node*> nodes, int N);
void grow_tree(std::vector<std::vector<Node*>>&tree, char*argv[], int lines);
void grow_tree_array(Node**tree, char*argv[], int lines);
void grow_tree_array_parallel(Node**tree, char*argv[], int lines, int start_col, int end_col);
void grow_tree_array_parallel2(Node*shtree, char*argv[], int lines, int start_col, int end_col);
void writeTreeToCsv(std::vector<std::vector<Node*>> tree);
void writeTreeToCsv_array(Node** tree, int cols, int gens);
void collision_check(Node**tree, int cols, int gens, int check_col, double last_deltaX);
void collision_check2(Node**tree, int cols, int gens, int check_col, double last_deltaX);
void writeTreeToCsv_array_single(Node* shtree, int cols, int gens);

int main(int argc, char*argv[]) {
    cout << "Number of arguments entered: " << argc << endl;
    for (int i=0;i<argc;i++){
        cout << "Argument "<< i << ": " << argv[i] << endl;
    }
    if(argc!=5) {
        cout << "Incorrect number of arguments." << endl;
        return 0;
    }
    int N = stoi(argv[2]);
    int generations = stoi(argv[3]);

    Node**testArray;
    FileIO testfio = FileIO();

    testArray=testfio.read_2d_array(argv);
    //Just printing the stuff from the initial file... This 7 comes from how many generations are given... will change based on file.
//    for (int i = 0; i < 7; ++i) {
//        for (int j = 0; j < N; ++j) {
//            cout<<"Gen: "<<testArray[j][i].generation<<"\tCol: "<<testArray[j][i].col<<"\tX: "<<testArray[j][i].x<<"\tY: "<<testArray[j][i].y<<"\tTheta: "<<testArray[j][i].theta<<"\tJ: "<<j<<"\tI: "<<i<<endl;
//        }
//    }

    //put 2d array into shared memory.
    int shmId; 			// ID of shared memory segment
    key_t shmKey = 123460; 		// key to pass to shmget(), key_t is an IPC key type defined in sys/types
    int shmFlag = IPC_CREAT | 0666; // Flag to create with rw permissions
    if ((shmId = shmget(shmKey, sizeof(Node)*N*(generations+7), shmFlag)) < 0)
    {
        std::cerr << "Init: Failed to initialize shared memory (" << shmId << ")" << std::endl;
        exit(1);
    }
    Node*shtree = (Node*)shmat(shmId, NULL, 0);
    for (int i = 0; i < generations+7; ++i) {
        for (int j = 0; j < N; ++j) {
            shtree[i*N+j]=testArray[j][i];
        }
    }
    //ensuring that it is being placed in shared memory...
    cout<<"HERE: "<<shtree[0].x<<endl;
    cout<<"HERE: "<<testArray[0][0].x<<endl;

    int startCol1 = 0;
    int stopCol1 = N/2;
    int startCol2 = N/2;
    int stopCol2 = N;
//    grow_tree_array_parallel(testArray, argv, 350, startCol1, stopCol1);
//    grow_tree_array_parallel(testArray, argv, 350, startCol2, stopCol2);
//    grow_tree_array_parallel2(shtree, argv, 350, startCol1, stopCol1);
//    grow_tree_array_parallel2(shtree, argv, 350, startCol2, stopCol2);
//    shtree[i*N+j]=testArray[j][i];
    for (int i = 0; i < generations+7-1; ++i) {
        for (int j = 0; j < N; ++j) {
            cout<<"SHTREE: "<<"Col: "<<shtree[i*N+j].col<<"\tGen: "<<shtree[i*N+j].generation<<endl;
        }
    }
    cout<<"PARALLEL GROW DONE"<<endl<<endl<<endl;
    pid_t pid;
    pid = fork();
    if ( pid < 0 )
    {
        std::cerr << "Could not fork!!! ("<< pid <<")" << std::endl;
        exit(1);
    }
    std::cout << "I just forked without error, I see ("<< pid <<")" << std::endl;

    if ( pid == 0 ) // Child process
    {
        grow_tree_array_parallel2(shtree, argv, 350, startCol1, stopCol1);

    }else{  //Parent process
        grow_tree_array_parallel2(shtree, argv, 350, startCol2, stopCol2);

    }
    std::cout << "In the parent: " << std::endl;

    int status;	// catch the status of the child

    do  // in reality, mulptiple signals or exit status could come from the child
    {

        pid_t w = waitpid(pid, &status, WUNTRACED | WCONTINUED);
        if (w == -1)
        {
            std::cerr << "Error waiting for child process ("<< pid <<")" << std::endl;
            break;
        }

        if (WIFEXITED(status))
        {
            if (status > 0)
            {
                std::cerr << "Child process ("<< pid <<") exited with non-zero status of " << WEXITSTATUS(status) << std::endl;
                continue;
            }
            else
            {
                std::cout << "Child process ("<< pid <<") exited with status of " << WEXITSTATUS(status) << std::endl;
                continue;
            }
        }
        else if (WIFSIGNALED(status))
        {
            std::cout << "Child process ("<< pid <<") killed by signal (" << WTERMSIG(status) << ")" << std::endl;
            continue;
        }
        else if (WIFSTOPPED(status))
        {
            std::cout << "Child process ("<< pid <<") stopped by signal (" << WSTOPSIG(status) << ")" << std::endl;
            continue;
        }
        else if (WIFCONTINUED(status))
        {
            std::cout << "Child process ("<< pid <<") continued" << std::endl;
            continue;
        }
    }
    while (!WIFEXITED(status) && !WIFSIGNALED(status));
//    grow_tree_array(testArray, argv, 350);
    //printing all generations generated plus the ones that were given in initial file.
//    for (int i = 0; i < generations+7; ++i) {
//        for (int j = 0; j < N; ++j) {
//            cout<<"Gen: "<<testArray[j][i].generation<<"\tCol: "<<testArray[j][i].col<<"\tX: "<<testArray[j][i].x<<"\tY: "<<testArray[j][i].y<<"\tTheta: "<<testArray[j][i].theta<<"\tJ: "<<j<<"\tI: "<<i<<endl;
//        }
//    }
    for (int i = 0; i < generations+7-1; ++i) {
        for (int j = 0; j < N; ++j) {
            cout<<"SHTREE: "<<"Col: "<<shtree[i*N+j].col<<"\tGen: "<<shtree[i*N+j].generation<<endl;
        }
    }
    writeTreeToCsv_array_single(shtree, N, generations);

//    writeTreeToCsv_array(testArray, N, generations);




//    init_data=fio.read_nc(argv, nc_count_ptr);
//    init_data=fio.read_nc(argv, nc_count_ptr);
//    cout<<"Lines: "<<nc_count<<endl;
//    for (int i = 0; i < init_data.size(); ++i) {
//        cout<<"x: "<<init_data.at(i)->x;
//        cout<<"\ty: "<<init_data.at(i)->y;
//        cout<<"\tcol: "<<init_data.at(i)->col;
//        cout<<"\tgen: "<<init_data.at(i)->generation;
//        cout<<"\ttheta: "<<init_data.at(i)->theta<<endl;
//    }
//    vector<vector<Node*>>tree = convert_2d(init_data, N);
//    for (int i = 0; i < N; ++i) {
//        for (int j = 0; j < nc_count/N; ++j) {
//            cout<<"["<<i<<"]"<<"["<<j<<"]";
//            cout<<"\tx: "<<tree[i][j]->x;
//            cout<<"\ty: "<<tree[i][j]->y;
//            cout<<"\tcol: "<<tree[i][j]->col;
//            cout<<"\tgen: "<<tree[i][j]->generation;
//            cout<<"\ttheta: "<<tree[i][j]->theta<<endl;
//        }
//    }
//    grow_tree(tree, argv, nc_count);
//    for (int i = 0; i < N; ++i) {
//        for (int j = 0; j < nc_count/N; ++j) {
//            cout<<"AFTER GROW ["<<i<<"]"<<"["<<j<<"]";
//            cout<<"\tx: "<<tree[i][j]->x;
//            cout<<"\ty: "<<tree[i][j]->y;
//            cout<<"\tcol: "<<tree[i][j]->col;
//            cout<<"\tgen: "<<tree[i][j]->generation;
//            cout<<"\ttheta: "<<tree[i][j]->theta<<endl;
//        }
//    }
    FileIO fio_t = FileIO();
//    Printer printer = Printer();
    vector<Node*> init_data_t;
    int nc_count_t=0;
    int*nc_count_ptr_t=&nc_count_t;
    init_data_t=fio_t.read_nc(argv, nc_count_ptr_t);
//    printer.print_tree_1d(init_data_t);
    cout<<"something"<<endl;
    vector<vector<Node*>>tree_t = convert_2d(init_data_t,N);
//    printer.print_tree_2d(tree_t);
    cout<<"something2"<<endl;


//     TO RUN AND TIME (this is using vector)
    FileIO fio = FileIO();
    vector<Node*> init_data;
    int nc_count=0;
    int*nc_count_ptr=&nc_count;
    auto ustart = chrono::high_resolution_clock::now();
    init_data=fio.read_nc(argv, nc_count_ptr);
    vector<vector<Node*>>tree = convert_2d(init_data, N);
    grow_tree(tree, argv, nc_count);
    auto ustop = chrono::high_resolution_clock::now();
    auto uduration = chrono::duration_cast<chrono::seconds>(ustop-ustart);
    cout<<"Time (s): "<<uduration.count()<<endl;
//    writeTreeToCsv(tree);
}

std::vector<std::vector<Node*>> convert_2d(std::vector<Node*> nodes, int N){
    std::vector<std::vector<Node*>> cnt_tree;
    unsigned int size = nodes.size();
    if(nodes.empty()){
        return cnt_tree;
    }
    for(int i = 0; i<N; i++){
        std::vector<Node*> col;
        for(unsigned int j = 0; j<size; j+= N){
            col.push_back(nodes[i+j]);
        }
        std::reverse(col.begin(),col.end());
        cnt_tree.push_back(col);
    }
    return cnt_tree;
}

void grow_tree(std::vector<std::vector<Node*>>&tree, char*argv[], int lines){
    int N = stoi(argv[2]);
    int generations=stoi(argv[3]);
    for (int i = 0; i < generations; ++i) {
        // cout<<"Generation: "<<i<<endl;
        for (int j = 0; j < N; ++j) {
            Node*prevprev = tree[j][1];
            Node*prev=tree[j][0];
            double deltaX = prevprev->x - prev->x;
            double deltaY = prevprev->y;
//            cout<<"dX: "<<deltaX<<"\t";
//            cout<<"dY: "<<deltaY<<endl;
            Node*newNode= new Node();
            newNode->x=prev->x;
            newNode->y=0;
            newNode->theta=prev->theta;
            newNode->generation=prev->generation+1;
            newNode->col=j;
            tree[j].insert(tree[j].begin(), newNode);
            double epsilon = 5e-8;
            bool coll = false;
            //propagate changes
            //this is to propogate changes without checking for collisions
//            for (int k = 1; k < i+(lines/N); ++k) {
//                tree[j][k]->x+=deltaX;
//                tree[j][k]->y+=deltaY;
//            }
            //this is to propogate changes and check for collisions.
            for (int k = 1; k < i+(lines/N); ++k) {
                if(coll){
                    tree[j][k]->x = tree[j][k-1]->x;
                    tree[j][k]->theta=90;
                }else{
                    if(tree[j][k]->theta==90){
                        //dont change by delta x if the angle is 90
                    }else {
                        tree[j][k]->x += deltaX;//else increment x by delta x to reflect x change
                    }
                }
                if(j!=0){
                    if((tree[j][k-1]->x - tree[j][k]->x)<epsilon){
                        coll=true;
                        tree[j][k]->x-=deltaX;
                        tree[j][k]->theta=90;
                    }
                }
//                double currX=tree[j][k]->x;
//                if(currX>maxX){
//                    maxX=currX;
//                }
                tree[j][k]->y+=deltaY;
            }
        }
//        cout<<"Gen done"<<endl;
    }
//    return vector<vector<Node*>>();
}
void writeTreeToCsv(std::vector<std::vector<Node*>> tree){
    FILE *fp = fopen("nodes.csv","w");
    if(!fp){
        printf("Failed to open file\n");
        return;
    }
    for(long unsigned int i = 0; i<tree.size(); i++){
        for(long unsigned int j = 0; j<tree[i].size(); j++){
            fprintf(fp,"%e,%e\n",tree[i][j]->x,tree[i][j]->y);
        }

    }
    return;
}
void writeTreeToCsv_array(Node** tree, int cols, int gens){
    FILE *fp = fopen("nodes.csv","w");
    if(!fp){
        printf("Failed to open file\n");
        return;
    }
    for(long unsigned int i = 0; i<gens; i++){
        for(long unsigned int j = 0; j<cols; j++){
            fprintf(fp,"%e,%e\n",tree[j][i].x,tree[j][i].y);
        }
    }
    return;
}
void writeTreeToCsv_array_single(Node* shtree, int cols, int gens){
    //    shtree[i*N+j]=testArray[j][i];
    int N = 50;
    FILE *fp = fopen("nodes.csv","w");
    if(!fp){
        printf("Failed to open file\n");
        return;
    }
    for(long unsigned int i = 0; i<gens; i++){
        for(long unsigned int j = 0; j<cols; j++){
            fprintf(fp,"%e,%e\n",shtree[i*N+j].x,shtree[i*N+j].y);
        }
    }
    return;
}

void grow_tree_array(Node**tree, char*argv[], int lines){
    int N = stoi(argv[2]);
    int generations=stoi(argv[3]);
    int existingGens = lines/N;
    for (int i = 0; i < generations; ++i) {
        // cout<<"Generation: "<<i<<endl;
        double last_deltaX;
        for (int j = 0; j < N; ++j) {   //PARALLELIZE THIS.... THIS CAN ALL HAPPEN AT ONCE. AFTER, PAUSE ALL PARALLELIZATION AND CHECK FOR COLLISIONS. IF COLLISION, CORRECT. THEN CONTINUE.
            Node prevprev = tree[j][i+existingGens-2];
            Node prev=tree[j][i+existingGens-1];
            double deltaX = prevprev.x - prev.x;
            last_deltaX=deltaX;
            double deltaY = prevprev.y - prev.y;
//            cout<<"dX: "<<deltaX<<"\t";
//            cout<<"dY: "<<deltaY<<endl;
            tree[j][existingGens+i].x=prev.x;
            tree[j][existingGens+i].y=0;
            tree[j][existingGens+i].generation=prev.generation+1;
            tree[j][existingGens+i].col=j;
            tree[j][existingGens+i].theta=prev.theta;
            //propagate changes
            //this is to propogate changes without checking for collisions
            for (int k = existingGens+i-1; k > 0; --k) {
                tree[j][k].x+=deltaX;
                tree[j][k].y+=deltaY;
            }
//            for (int k = 0; k < existingGens+i-1; ++k) {
//                tree[j][k].x+=deltaX;
//                tree[j][k].y+=deltaY;
//            }
//            collision_check2(tree, N, generations, existingGens+i, last_deltaX);
        }
        cout<<"Gen "<<i<<" done."<<endl;
//        cout<<"Checking for collisions..."<<endl;
        collision_check(tree, N, generations, existingGens+i, last_deltaX);
    }
//    return vector<vector<Node*>>();
}
void grow_tree_array_parallel(Node**tree, char*argv[], int lines, int start_col, int end_col){
    int N = stoi(argv[2]);
    int generations=stoi(argv[3]);
    int existingGens = lines/N;
    for (int i = 0; i < generations; ++i) {
        // cout<<"Generation: "<<i<<endl;
        double last_deltaX;
        for (int j = start_col; j < end_col; ++j) {   //PARALLELIZE THIS.... THIS CAN ALL HAPPEN AT ONCE. AFTER, PAUSE ALL PARALLELIZATION AND CHECK FOR COLLISIONS. IF COLLISION, CORRECT. THEN CONTINUE.
            Node prevprev = tree[j][i+existingGens-2];
            Node prev=tree[j][i+existingGens-1];
            double deltaX = prevprev.x - prev.x;
            last_deltaX=deltaX;
            double deltaY = prevprev.y - prev.y;
//            cout<<"dX: "<<deltaX<<"\t";
//            cout<<"dY: "<<deltaY<<endl;
            tree[j][existingGens+i].x=prev.x;
            tree[j][existingGens+i].y=0;
            tree[j][existingGens+i].generation=prev.generation+1;
            tree[j][existingGens+i].col=j;
            tree[j][existingGens+i].theta=prev.theta;
            //propagate changes
            //this is to propogate changes without checking for collisions
            for (int k = existingGens+i-1; k > 0; --k) {
                tree[j][k].x+=deltaX;
                tree[j][k].y+=deltaY;
            }
//            for (int k = 0; k < existingGens+i-1; ++k) {
//                tree[j][k].x+=deltaX;
//                tree[j][k].y+=deltaY;
//            }
//            collision_check2(tree, N, generations, existingGens+i, last_deltaX);
        }
        cout<<"Gen "<<i<<" done."<<endl;
//        cout<<"Checking for collisions..."<<endl;
        collision_check(tree, N, generations, existingGens+i, last_deltaX);
    }
//    return vector<vector<Node*>>();
}
void grow_tree_array_parallel2(Node*shtree, char*argv[], int lines, int start_col, int end_col){
//    shtree[i*N+j]=testArray[j][i];

    int N = stoi(argv[2]);
    int generations=stoi(argv[3]);
    int existingGens = lines/N;
    for (int i = 0; i < generations; ++i) {
        // cout<<"Generation: "<<i<<endl;
//        double last_deltaX;
        for (int j = start_col; j < end_col; ++j) {   //PARALLELIZE THIS.... THIS CAN ALL HAPPEN AT ONCE. AFTER, PAUSE ALL PARALLELIZATION AND CHECK FOR COLLISIONS. IF COLLISION, CORRECT. THEN CONTINUE.
//            Node prevprev = tree[j][i+existingGens-2];
            Node prevprev = shtree[(i+existingGens-2)*N+j];
//            Node prev=tree[j][i+existingGens-1];
            Node prev=shtree[(i+existingGens-1)*N+j];

            double deltaX = prevprev.x - prev.x;
//            last_deltaX=deltaX;
            double deltaY = prevprev.y - prev.y;
//            cout<<"dX: "<<deltaX<<"\t";
//            cout<<"dY: "<<deltaY<<endl;
//            tree[j][existingGens+i].x=prev.x;
            shtree[(existingGens+i)*N+j].x=prev.x;
//            tree[j][existingGens+i].y=0;
            shtree[(existingGens+i)*N+j].y=0;
//            tree[j][existingGens+i].generation=prev.generation+1;
            shtree[(existingGens+i)*N+j].generation=prev.generation+1;
//            tree[j][existingGens+i].col=j;
            shtree[(existingGens+i)*N+j].col=j;
//            tree[j][existingGens+i].theta=prev.theta;
            shtree[(existingGens+i)*N+j].theta=prev.theta;
            //propagate changes
            //this is to propogate changes without checking for collisions
            for (int k = existingGens+i-1; k > 0; --k) {
//                tree[j][k].x+=deltaX;
                shtree[(k*N)+j].x+=deltaX;
//                tree[j][k].y+=deltaY;
                shtree[(k*N)+j].y+=deltaY;
            }
//            for (int k = 0; k < existingGens+i-1; ++k) {
//                tree[j][k].x+=deltaX;
//                tree[j][k].y+=deltaY;
//            }
//            collision_check2(tree, N, generations, existingGens+i, last_deltaX);
        }
        cout<<"Gen "<<i<<" done."<<endl;
//        cout<<"Checking for collisions..."<<endl;
//        collision_check(tree, N, generations, existingGens+i, last_deltaX);
    }
//    return vector<vector<Node*>>();
}

void collision_check(Node**tree, int cols, int gens, int current_gen, double last_deltaX){
    double epsilon = 5e-8;
    for (int i = 0; i < cols; ++i) {
//        cout<<"i: "<<i<<endl;
        Node colTop = tree[i][0];
        if((i+1)<cols){//if there is a next col... aka if we arent at the end of the cols for the generation.
            Node* nextCol = tree[i+1];
            for (int j = 0; j < current_gen; ++j) {
//                cout<<"j: "<<j<<endl;
                if((nextCol[j].x - colTop.x)<=epsilon){//if the x values are close enough to have a potential collision
//                    cout<<"POTENTIAL COLLISION"<<endl;
                    //use the distance formula to check for collision. (now take the y values into account)
                    double d = sqrt( pow((nextCol[j].x - colTop.x),2) + pow((nextCol[j].y - colTop.y),2));
//                    cout<<d<<endl;
                    if(d<=epsilon){
                        cout<<"Collision detected."<<endl;
                        cout<<"Correcting cols..."<<endl;
                        //propogate changes up... correct for collision.
                        tree[i][current_gen].theta=90;
                        for (int k = current_gen; k > 0; --k) {
                            tree[j][k].x-=last_deltaX;
                        }
                    }
                }
//                cout<<"NO COLLISION"<<endl;
            }
        }else{
            return;//return if we are at the end.
        }

    }
}
void collision_check2(Node**tree, int cols, int gens, int current_gen, double last_deltaX){
    double epsilon = 5e-8;
    for (int i = 0; i < cols; ++i) {
//        cout<<"i: "<<i<<endl;
        Node colTop = tree[i][0];
        if((i+1)<cols){//if there is a next col... aka if we arent at the end of the cols for the generation.
            Node* nextCol = tree[i+1];
            for (int j = 0; j < current_gen; ++j) {
//                cout<<"j: "<<j<<endl;
                if((nextCol[j].x - colTop.x)<=epsilon){//if the x values are close enough to have a potential collision
//                    cout<<"POTENTIAL COLLISION"<<endl;
                    //use the distance formula to check for collision. (now take the y values into account)
                    double d = sqrt( pow((nextCol[j].x - colTop.x),2) + pow((nextCol[j].y - colTop.y),2));
//                    cout<<d<<endl;
                    if(d<=epsilon){
                        cout<<"Collision detected."<<endl;
                        cout<<"Correcting cols..."<<endl;
                        //propogate changes up... correct for collision.
                        tree[i][current_gen+7-1].theta=90;
                        for (int k = 7+current_gen-1; k > 0; --k) {
                            tree[j][k].x-=last_deltaX;
                        }
                    }
                }
//                cout<<"NO COLLISION"<<endl;
            }
        }else{
            return;//return if we are at the end.
        }

    }
}
