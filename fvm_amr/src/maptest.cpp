#include <iostream>
#include <fstream>
#include <cstring>
#include <random>
#include <cmath>
#include <bitset>
#include <map>


using uint = unsigned int;
template <size_t N>
using morton = std::bitset<N>;

template <size_t N>
class compare {
    public:
        bool operator()(const morton<N> &a, const morton<N> &b) const 
        {
            return a.to_ullong() < b.to_ullong();
        }

};

template <size_t N, typename mtype>
using Mesh = std::map<morton<N>, mtype *, compare<N>>;

#define SIZE 8


class Cell {
    double value;
    bool leaf ;
    bool refine = false;

    public:
    double xmin[2], xmax[2];
    uint nx, ny;
    int level;

    Cell(double a) {value = a; leaf=true;};
    double get_value(void) {return value;};
    void set_value(double a) {value = a;};
    bool get_leaf(void) {return leaf;};
    bool get_refine(void) {return refine;};
    void off(void) {refine = false;};
    void on(void) {refine = true;};
};

template<size_t N, typename mtype>
uint key_to_level(morton<N> key, Mesh<N,mtype> mesh) {
    /* Get refinement level */
    uint lvl = N/2;
    morton<N> temp ;
    //std::cout << "\nTesting " << key << ", l=" << lvl << "\n";
    for(uint i=0; i < N-1; i+=2) {
        //std::cout << "\tat " << key[i] << key[i+1] << ", l=" << lvl;
        if ((key[i] == false) && (key[i + 1] == false)) {
            /* 00 on this level 
             * check if the 01 sibling exists
             */
            //std::cout << ", double zero";
            temp = key;
            temp.flip(i);
            //std::cout << " Checking for " << temp ;
            if (mesh.count(temp) == 0) {
                //std::cout << ", none\n";
                lvl -= 1;
            }
            else {
                //std::cout << ", found\n";
            }
        }
        else {
            //std::cout << ", Stop\n\tDone, l= " << lvl << "\n";
            return lvl;
        }

    }
    //std::cout << "Reached end, l=" << lvl << "\n";
    return lvl;

}

template<size_t N>
void ext_from_level(morton<N> key, uint level , double xmin[2], double xmax[2], uint *nx, uint *ny) {
    morton<N> xkey, ykey;
    size_t nhalf = 2*level;
    uint x,y;

    //std::cout << key;


    for(uint bit=0; bit < level; bit++) {
        xkey[level-1-bit] = key[N - 1 - 2*bit];
        ykey[level-1-bit] = key[N - 1 - 2*bit-1];
    }
    //for(uint i=0;i<N;i+=2) xkey[nhalf-i/2-1] = key[i];
    //for(uint i=1;i<N;i+=2) ykey[nhalf-i/2-1] = key[i];

    /*
    for(uint i=1;i<=nhalf;i++) {
        xkey[nhalf-i] = key[2*i];
        ykey[nhalf-i] = key[2*i+1];
    }
    */

    //std::cout << ", x: " << xkey.to_string() << "," << xkey.to_ulong();
    //std::cout << ", y: " << ykey.to_string()<< "," << ykey.to_ulong();
    x = xkey.to_ulong();
    y = ykey.to_ulong();

    *nx = x; 
    *ny = y;

    //std::cout << " => (" << level << "," << x << "," << y << ")\n";


    //double xmin[2], xmax[2];
    double dx = 1./pow(2, level);
    xmin[0] = dx*x ;
    xmax[0] = dx*(x+1);
    xmin[1] = dx*y ;
    xmax[1] = dx*(y+1);


    //std::cout << "\t dx=" << dx << " (" << xmin[0] << "," << xmin[1] << ") => (" << xmax[0] << "," << xmax[1] <<  ")\n";
}

template<size_t N>
void on_boundary(morton<N> key, uint level) {
    

}


template<size_t N, typename mtype>
void fix_neighbor(morton<N> key, uint level, morton<N> *list, Mesh<N,mtype> mesh, int dir, bool sib) {
    /* Check if neighbor is coarser or finer than key */
    /* If finer then list[0],list[1] will contain new neighbors */
    /* If coarser or same then list[1] will contain list[0] */

    int offset = N-1+2;
    morton<N> temp = list[0]; 
    if (mesh.count(temp) == 0) {
        /* Coarse */
        temp[offset - 2*level - 0] = 0;
        temp[offset - 2*level - 1] = 0;
        list[0] = temp;
        list[1] = temp;
        return;
    }
        
    if (level == N/2) {
        /* Already at max level */
        list[1] = list[0];

        return;
    }
    /* Check if sibs exist */
    if (mesh.count(temp.flip(offset - 2*(level + 1) - dir )) == 0) {
        /* No siblings so we're done */
        list[1] = list[0];
        return;
    }

    temp = list[0];
    if (sib) {
        temp[offset - 2*(level + 1) - 0] = key[offset-2*level - 0];
        temp[offset - 2*(level + 1) - 1] = key[offset-2*level - 1];
    }
    else {
        /* Repeat bits for cousins */
        temp[offset - 2*(level + 1) - 0] = temp[offset-2*level - 0];
        temp[offset - 2*(level + 1) - 1] = temp[offset-2*level - 1];

    }
    /* Flip bit in opposite direction */
    list[0] = temp;
    temp.flip(offset - 2*(level+1) - (dir + 1)%2);
    list[1] = temp;


    return;
}

template<size_t N, typename mtype>
void neighbor_list(morton<N> key, uint level, morton<N> *list, Mesh<N,mtype> mesh) {
    int i,k;
    morton<N> temp;


//    std::cout << "NEIGHBOR of " << key.to_string() << ", l=" << level << "\n";

//    std::cout << level<< "\n";

    int offset = N-1+2;
    /* First set neighbors on this level */

    /* x direction*/
    /* Sibling */
    temp = key;
    temp = temp.flip(offset - 2*level - 0 );
    list[0 + 2*0] = temp;
    fix_neighbor(key, level, &list[2*0], mesh,0, true);

    std::cout << "1\n";
    /* cousin */
    for(k=level-1; k > 1; k--) {
        if (key[offset - 2*k - 0] != key[offset - 2*level-0]) {
            break;
        }
    }
    if (k == 0) k = 1;
    std::cout << "l " << level << " k " << k << "\n";
    temp = key;
    for(i=k; i<=level; i++) {
        temp.flip(offset - 2*i - 0);
    }
    list[0 + 2*1] = temp;
    fix_neighbor(key, level, &list[2*1], mesh,0, false);



    std::cout << "1\n";


    /* y direction */
    temp = key;
    temp = temp.flip(offset - 2*level - 1 );
    list[0 + 2*2] = temp;
    fix_neighbor(key, level, &list[2*2], mesh,1,true);
    std::cout << "1\n";

    for(k=level-1; k > 1; k--) {
        if (key[offset - 2*k - 1] != key[offset - 2*level -1]) {
            break;
        }
    }
    if (k == 0) k = 1;
    temp = key;
    for(i=k; i<=level; i++) {
        temp.flip(offset - 2*i - 1);
    }
    list[0 + 2*3] = temp;
    fix_neighbor(key, level, &list[2*3], mesh,1,false);

    std::cout << "1\n";

    std::cout << key.to_string() << "\n";
    for(int i=0;i<4;i++) {
        for(k=0;k<2;k++) {
            std::cout << " " <<  list[k+2*i].to_string();
            //if (list[k + 2*i] != nullptr) {
            //    std::cout << " " << list[k + 2*i].to_string() ;
            //    if (i % 2 == 0) std::cout << " x";
            //    else std::cout << " y";
            //}
        }
        std::cout << "\n";
    }



    return;

}
template<size_t N, typename mtype>
void split(morton<N> mykey, Mesh<N,mtype> *mesh) {

    morton<N> key = mykey;

    Cell *child = mesh->at(key);
    int level = child->level;
    if (level == N/2) {
        std::cout << "Hit max level!\n";
        return;
    }

    //uint level = key_to_level(key,*mesh);
    //std::cout << "Splitting " << key.to_string() << ", level " << level << "\n";

    uint offset = N - 2*level;

    child->off();
    double value = child->get_value();

    value /= 4;

    child->set_value(value);
    child->level = level + 1;
    ext_from_level(key, child->level,child->xmin, child->xmax, &child->nx, &child->ny);

    std::cout << key << "\n";
    std::cout << "\t" << key << " "<< child->xmin[0] << " " << child->xmin[1] << " "<< child->xmax[0] << " " << child->xmax[1] << "\n";


    /* xx10 */

    key.flip(offset-1);
    child = new Cell(value);
    child->level = level + 1;
    ext_from_level(key, child->level,child->xmin, child->xmax, &child->nx, &child->ny);
    std::cout << "\t" << key << " "<< child->xmin[0] << " " << child->xmin[1] << " "<< child->xmax[0] << " " << child->xmax[1] << "\n";
    mesh->insert({key, child});

    /* xx11 */
    key.flip(offset-2);
    child = new Cell(value);
    child->level = level + 1;
    ext_from_level(key, child->level,child->xmin, child->xmax, &child->nx, &child->ny);
    std::cout << "\t" << key << " "<< child->xmin[0] << " " << child->xmin[1] << " "<< child->xmax[0] << " " << child->xmax[1] << "\n";
    mesh->insert({key, child});


    /* xx01 */
    key.flip(offset-1);
    child = new Cell(value);
    child->level = level + 1;
    ext_from_level(key, child->level,child->xmin, child->xmax, &child->nx, &child->ny);
    std::cout << "\t" << key << " "<< child->xmin[0] << " " << child->xmin[1] << " "<< child->xmax[0] << " " << child->xmax[1] << "\n";
    mesh->insert({key, child});
    return;
}

void print_mesh(Mesh<SIZE,Cell> mesh) {
    Mesh<SIZE,Cell>::iterator it;
    for (it=mesh.begin(); it!=mesh.end(); ++it) {
        std::cout << (it->first).to_string();
        std::cout << " " << it->second->level; 
        std::cout << " " << it->second->nx;
        std::cout << " " << it->second->ny;
        std::cout << " " << it->second->xmin[0];
        std::cout << " " << it->second->xmin[1];
        std::cout << " " << it->second->xmax[0];
        std::cout << " " << it->second->xmax[1];
        std::cout << "\n";
    }
}

void output_mesh(char *fname, char *fname2, Mesh<SIZE,Cell> mesh) {

    std::ofstream myfile,myfile2;
    myfile.open (fname);
    myfile2.open (fname2);

    myfile << "#i l nx ny xmin ymin xmax ymax\n";
    for (Mesh<SIZE,Cell>::iterator it=mesh.begin(); it!=mesh.end(); ++it) {
        myfile2 << (it->first).to_string() << "\n";
        myfile << (it->first).to_ulong();
        myfile << " " << it->second->level; 
        myfile << " " << it->second->nx;
        myfile << " " << it->second->ny;
        myfile << " " << it->second->xmin[0];
        myfile << " " << it->second->xmin[1];
        myfile << " " << it->second->xmax[0];
        myfile << " " << it->second->xmax[1];
        myfile << "\n";
        //std::cout << indx << "," << dat[indx] <<  ", (" << datx[indx] << ", " << datx[indx] << ")\n" ;
    }
    myfile.close();
    myfile2.close();
    return;
}
void output_neighbors(char *fname, morton<SIZE> key, morton<SIZE> *list, Mesh<SIZE,Cell> mesh) {

    Cell *cell = mesh.at(key);
    std::ofstream myfile;
    myfile.open (fname);

    myfile << "#i xc yc\n";
    myfile <<  key.to_ulong();
    myfile << " " << ( cell->xmin[0] + cell->xmax[0])/2.;
    myfile << " " << ( cell->xmin[1] + cell->xmax[1])/2.;
    myfile << "\n";
    for(uint i=0;i<4;i++) {
        for(uint j=0;j<2;j++) {
            //if (list[j + 2*i] != list[j + 2*i]) {
                cell = mesh.at(list[j + 2*i]);
                myfile <<  list[j+2*i].to_ulong();
                myfile << " " << ( cell->xmin[0] + cell->xmax[0])/2.;
                myfile << " " << ( cell->xmin[1] + cell->xmax[1])/2.;
                myfile << "\n";
            //}
        //std::cout << indx << "," << dat[indx] <<  ", (" << datx[indx] << ", " << datx[indx] << ")\n" ;
        }
    }
    myfile.close();
    return;
}

int main(int argc, char *argv[]) {
    morton<SIZE> key;

    Mesh<SIZE, Cell> mesh;

    Cell *val = new Cell(1.);
    val->level = 0;
    ext_from_level(key, val->level,val->xmin, val->xmax, &val->nx, &val->ny);
    mesh.insert({key,val});

    split(key, &mesh);

    //val = new Cell(10.);
    //key.flip(SIZE-1);
    //mesh.insert({key, val});


    //val = new Cell(4.);
    //key.flip(SIZE-1); key.flip(SIZE-2);
    //mesh.insert({key, val});


    //val = new Cell(20.);
    //key.flip(SIZE-1); 
    //mesh.insert({key, val});

    std::cout << "----------------------------\n";
    std::cout << "Before split\n";
    //print_mesh(mesh);

    Mesh<SIZE,Cell>::iterator it;
    uint ntot = 0;


   // split(key,&mesh);
   
    key.reset();
    split(key,&mesh);

    key.reset();
    key.flip(SIZE-2); 
    split(key,&mesh);
    key.flip(SIZE-1);
    split(key,&mesh);
    
    key.flip(SIZE-2);

    split(key,&mesh);

   // key.flip(SIZE-2);

   // key.reset();
   // split(key,&mesh);
    
    //key.reset();
    //key.flip(SIZE-1);
    //key.flip(SIZE-2);
    //key.flip(SIZE-4);
    //split(key,&mesh);

    /*
    key.flip(SIZE-3); key.flip(SIZE-4); 
    split(key,&mesh);


    key.flip(SIZE-5); key.flip(SIZE-6);
    split(key,&mesh);
    */

  //  ntot = 0;
  //  for (it=mesh.begin(); it!=mesh.end(); ++it) {
  //      ntot += 1;
  //      ext_from_level(it->first, key_to_level(it->first, mesh),
  //              it->second->xmin, it->second->xmax, &it->second->nx, &it->second->ny);
  //  }


    std::cout << "Output\n";
    print_mesh(mesh);
    output_mesh("mesh.txt", "names.txt",mesh);
 

    key.reset();
    key.flip(SIZE-3);
    key.flip(SIZE-4);
    split(key,&mesh);

    key.reset();
    key.flip(SIZE-2);
    key.flip(SIZE-3);
    split(key,&mesh);

    key.reset();
    key.flip(SIZE-1);
    key.flip(SIZE-2);
    split(key,&mesh);

    key.reset();
    key.flip(SIZE-1);
    key.flip(SIZE-4);
    split(key,&mesh);


    key.reset();
    key.flip(SIZE-2);
    key.flip(SIZE-3);
    key.flip(SIZE-5);
    split(key,&mesh);

    key.reset();
    key.flip(SIZE-1);
    key.flip(SIZE-2);
    split(key,&mesh);

    key.reset();
    key.flip(SIZE-1);
    key.flip(SIZE-4);
    key.flip(SIZE-6);
    split(key,&mesh);

    key.reset();
    key.flip(SIZE-3);
    key.flip(SIZE-4);
    key.flip(SIZE-5);
    key.flip(SIZE-6);
    split(key,&mesh);

 //   ntot = 0;
//    for (it=mesh.begin(); it!=mesh.end(); ++it) {
//        ntot += 1;
//        ext_from_level(it->first, key_to_level(it->first, mesh),
//                it->second->xmin, it->second->xmax, &it->second->nx, &it->second->ny);
//    }




    std::cout << "Output\n";
    output_mesh("mesh2.txt", "names2.txt",mesh);

/*
    key.reset();
    key.flip(SIZE-1);
    key.flip(SIZE-3);
    std::cout << "Finding neigh of " << key.to_string() << "\n";
    */

    morton<SIZE> neighbors[8];
    neighbor_list(key, key_to_level(key,mesh), neighbors, mesh);

    output_neighbors("neigh.txt", key, neighbors, mesh);

        



}
