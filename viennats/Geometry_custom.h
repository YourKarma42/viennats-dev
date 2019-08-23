#ifndef GEOMETRY_CUSTOM_H_
#define GEOMETRY_CUSTOM_H_

/* =========================================================================
Copyright (c)    2008-2015, Institute for Microelectronics, TU Wien.

-----------------
ViennaTS - The Vienna Topography Simulator
-----------------

Contact:         viennats@iue.tuwien.ac.at

License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */

#include <vector>
#include <string>
#include <fstream>
#include <bitset>
#include <set>
#include "message.h"
#include "vector.hpp"
#include "output.hpp"
#include "misc.hpp"
#include <stdexcept>

//json parser
#include "nlohmann/json.hpp"



namespace geometry {

    using json = nlohmann::json;



    //https://en.wikipedia.org/wiki/Midpoint_circle_algorithm
    //only 2D think of 3D algorithm
    // midpoint circle algorithm 
    //TODO: grid uebergeben 
    //doc machen
    template<int D, class LevelSetType, class index_type, class value_type>
    bool create_circle(double radius, 
                       double center[D],
                       std::list<LevelSetType>& LevelSets,
                       std::vector< std::pair< lvlset::vec<index_type,D>, value_type> > & points,
                       double g_d){

        //the position of the point is calculated at the center and then transformed to the grid position

        //TODO: check if circle radius and grid delta are high enough to create a circle

        //init recursion
        //double g_d = 0.5;
        double r_sq = (radius/g_d) * (radius/g_d);
        //TODO: umschreiben mit arrays wegen 3 dim
        double x_n_sq = r_sq;
        double y_n_sq = 0.0;

        double x_n = std::sqrt(x_n_sq);
        double y_n = std::sqrt(y_n_sq);

        double tmp_const = 1.0;

        //TODO: reserve space
        std::vector<std::pair<double, double>> circle_points;

        bool run = true;

        std::vector< std::pair< lvlset::vec<index_type,D>, value_type> > points_floor;

        std::vector< std::pair< lvlset::vec<index_type,D>, value_type> > points_ceil;



        //calculate all points of the circle by mirroring the first octant
        //circle octants counter clockwise
        //move calculated (local) coordinates to global coordinates
        while(run){


            if(!((floor(x_n) != floor(y_n)) && (ceil(x_n) != floor(y_n))))
                run = false;


            if(y_n != 0.0 ){

                circle_points.push_back(std::make_pair(-y_n + center[0], -x_n + center[1]));

                circle_points.push_back(std::make_pair(-x_n + center[0], -y_n + center[1]));

                circle_points.push_back(std::make_pair(x_n + center[0], -y_n + center[1])); 
          
                circle_points.push_back(std::make_pair(-y_n + center[0], x_n + center[1]));              

            }


            lvlset::vec<index_type,D> pos_grid1;
            lvlset::vec<index_type,D> pos_grid2;

            circle_points.push_back(std::make_pair(x_n + center[0], y_n + center[1]));

            circle_points.push_back(std::make_pair(y_n + center[0], x_n + center[1])); 

            circle_points.push_back(std::make_pair(-x_n + center[0], y_n + center[1]));   

            circle_points.push_back(std::make_pair(y_n + center[0], -x_n + center[1]));  
          

            //calculate the next point on the circle

            x_n_sq = x_n_sq - 2*y_n - tmp_const;
            y_n_sq = y_n_sq + 2*y_n + tmp_const;

            x_n = std::sqrt(x_n_sq);
            y_n = std::sqrt(y_n_sq);

        }







        // diese Schleife in die obere packen
        //iterate over all points to calculate lvl set points and values
        for(auto point: circle_points){

            lvlset::vec<index_type,D> pos_grid1;
            lvlset::vec<index_type,D> pos_grid2;

            double dist_ceil = 0.0;
            double dist_floor = 0.0;

            
            if(abs(ceil(point.first)) - abs(floor(point.first)) == 0.0){
                //the current point lies on the y line or directly on the grid
                
                if(abs(ceil(point.second)) - abs(floor(point.second)) == 0.0){
                    //the current point lies directly on a grid point

                    lvlset::vec<index_type,D> pos_grid3;

                    //std::cout << floor(point.first)- center[0] << std::endl;

                    //std::cout << floor(radius/g_d)<< std::endl;


                    if(abs(floor(point.first)- center[0]) == floor(radius/g_d)){

                        pos_grid1[0]=floor(point.first)-1;
                        pos_grid1[1]=floor(point.second);

                        pos_grid2[0]=floor(point.first)+1;
                        pos_grid2[1]=floor(point.second);

                        if(point.first > center[0]){
                            dist_floor = -0.75;
                            dist_ceil  = 0.75;
                        }else{
                            dist_floor = 0.75;
                            dist_ceil  = -0.75;
                        }

                    }else if (abs(floor(point.second)- center[1]) == floor(radius/g_d)){

                        pos_grid1[0]=floor(point.first);
                        pos_grid1[1]=floor(point.second)-1;

                        pos_grid2[0]=floor(point.first);
                        pos_grid2[1]=floor(point.second)+1;

                        if(point.second > center[1]){
                            dist_floor = -0.75;
                            dist_ceil  = 0.75;
                        }else{
                            dist_floor = 0.75;
                            dist_ceil  = -0.75;
                        }

                    }else{
                        //HACK nix gut
                        pos_grid1[0]=floor(point.first)-1;
                        pos_grid1[1]=floor(point.second);

                        pos_grid2[0]=floor(point.first)+1;
                        pos_grid2[1]=floor(point.second);

                        if(point.first > center[0]){
                            dist_floor = -0.75;
                            dist_ceil  = 0.75;
                        }else{
                            dist_floor = 0.75;
                            dist_ceil  = -0.75;
                        }


                    }



                    pos_grid3[0]=floor(point.first);
                    pos_grid3[1]=floor(point.second);

                    //points.push_back(std::make_pair(LevelSets.begin()->grid().global_indices_2_local_indices(pos_grid1), dist_floor));
                    //points.push_back(std::make_pair(LevelSets.begin()->grid().global_indices_2_local_indices(pos_grid2), dist_ceil));
                    //points.push_back(std::make_pair(LevelSets.begin()->grid().global_indices_2_local_indices(pos_grid3), 0.0));

                    points_floor.push_back(std::make_pair(pos_grid1, dist_floor));
                    points_ceil.push_back(std::make_pair(pos_grid2, dist_ceil));
                    points_ceil.push_back(std::make_pair(pos_grid3, 0.0));
                    //points.push_back(std::make_pair(pos_grid3, 0.0));

                    std::cout << "cord1:" << pos_grid1 << "    " << dist_floor << std::endl;
                    std::cout << "cord2:" << pos_grid2 << "    " << dist_ceil << std::endl;
                    
                    

                }else{


                    //the current point lies on the y line


                    pos_grid1[0]=floor(point.first);
                    pos_grid1[1]=floor(point.second);

                    dist_floor = abs(abs(floor(point.second))-abs(point.second));

                    pos_grid2[0]=ceil(point.first);
                    pos_grid2[1]=ceil(point.second);

                    dist_ceil = abs(abs(ceil(point.second))-abs(point.second));

                    if(point.second < center[1]){
                        
                        dist_ceil = dist_ceil*(-1);

                    }else{
                        
                        dist_floor = dist_floor*(-1);
                    }

                    //points.push_back(std::make_pair(LevelSets.begin()->grid().global_indices_2_local_indices(pos_grid1), dist_floor));
                    //points.push_back(std::make_pair(LevelSets.begin()->grid().global_indices_2_local_indices(pos_grid2), dist_ceil));

                    points_floor.push_back(std::make_pair(pos_grid1, dist_floor));
                    points_ceil.push_back(std::make_pair(pos_grid2, dist_ceil));




                }

            }else{
                //the current point lies on the x line

                pos_grid1[0]=floor(point.first);
                pos_grid1[1]=floor(point.second);

                dist_floor = abs(abs(floor(point.first))-abs(point.first));

                pos_grid2[0]=ceil(point.first);
                pos_grid2[1]=ceil(point.second);

                dist_ceil = abs(abs(ceil(point.first))-abs(point.first));

                if(point.first < center[0]){
                        
                    dist_ceil = dist_ceil*(-1);

                }else{
                        
                    dist_floor = dist_floor*(-1);
                }

                //points.push_back(std::make_pair(LevelSets.begin()->grid().global_indices_2_local_indices(pos_grid1), dist_floor));
                //points.push_back(std::make_pair(LevelSets.begin()->grid().global_indices_2_local_indices(pos_grid2), dist_ceil));

                points_floor.push_back(std::make_pair(pos_grid1, dist_floor));
                points_ceil.push_back(std::make_pair(pos_grid2, dist_ceil));


            }  
            
        }


        


        std::vector< std::pair< lvlset::vec<index_type,D>, value_type> > all_points = points_floor;

        std::vector< std::pair< lvlset::vec<index_type,D>, value_type> > deleted_points;

        all_points.insert(all_points.end(), points_ceil.begin(), points_ceil.end());
       
        bool first = true;

        for(auto p: all_points){
            if(p.first[0] == 9)
                std::cout << p.first << " ";
        }
        std::cout << std::endl;

        auto p = all_points.begin();

        while(p != all_points.end()){

            //std::pair< lvlset::vec<index_type,D>, value_type> c_p = *p;

            first = true;
            bool del_p = false;

            auto p_i = all_points.begin();

            while(p_i != all_points.end()){        

                if(p->first == p_i->first){

                   if(!first){
                        all_points.erase(p_i);
                        del_p = true; 
                    }else{
                        p_i++;
                        first = false; 
                    }                             
                }else{                
                    p_i++;
                }
            }
            if(del_p){
                deleted_points.push_back(*p);
                all_points.erase(p);

            }else{
                p++;
            }

        }

        for(auto p: all_points){
            if(p.first[0] == 9)
                std::cout << p.first << " ";
        }
        std::cout << std::endl;
 


       // std::cout << anz << std::endl;

        std::cout << points_floor.size() << std::endl; 

        std::cout << points_ceil.size() << std::endl;  

        string test_x = "[";

        string test_y = "[";

        for(auto p : all_points){
            points.push_back(std::make_pair(LevelSets.begin()->grid().global_indices_2_local_indices(p.first), p.second));

            test_x += std::to_string(p.first[0]) + ", ";

            test_y += std::to_string(p.first[1]) + ", ";
        }



        std::cout << "x= " << test_x.substr(0, test_x.size()-2) + "]" << std::endl;

        std::cout << "y= " << test_y.substr(0, test_y.size()-2) + "]"<< std::endl;

        test_x = "[";

        test_y = "[";

        for(auto p: deleted_points){

            double tmp_radius = std::sqrt(pow(p.first[0] - center[0], 2) + pow(p.first[1] - center[1], 2));

            double dist = tmp_radius - (radius/g_d);

            points.push_back(std::make_pair(LevelSets.begin()->grid().global_indices_2_local_indices(p.first), dist));


            /*if(tmp_radius - (radius/g_d) >= 0){
                dist = tmp_radius - (radius/g_d);
            }else{
                dist = (radius/g_d) - tmp_radius;
            } */ 

            test_x += std::to_string(p.first[0]) + ", ";

            test_y += std::to_string(p.first[1]) + ", ";

        }

        //std::cout << "x= " << test_x.substr(0, test_x.size()-2) + "]" << std::endl;

        //std::cout << "y= " << test_y.substr(0, test_y.size()-2) + "]"<< std::endl;





        test_x = "[";

        test_y = "[";

        for(auto p: circle_points){

            test_x += std::to_string(p.first) + ", ";

            test_y += std::to_string(p.second) + ", ";

        }

        //std::cout << "x= " << test_x.substr(0, test_x.size()-2) + "]" << std::endl;

        //std::cout << "y= " << test_y.substr(0, test_y.size()-2) + "]"<< std::endl;

        return true;

    }



    template<int D, class LevelSetType, class index_type, class value_type>
    bool create_circle(double radius, 
                       double center[D],
                       std::list<LevelSetType>& LevelSets,
                       std::vector< std::pair< lvlset::vec<index_type,D>, value_type> > & points){



        return true;



    }


    //TODO: at the moment hardcoded for 2 dimensions
    template<int D, class GridTraitsType, class ParameterType, class LevelSetType>
    void my_geom_test_func(GridTraitsType& GridProperties, lvlset::grid_type<GridTraitsType>& grid,
                                      ParameterType& p, std::list<LevelSetType>& LevelSets){
        
        //probably have to create the lvlset 2 times for MC

        typedef typename LevelSetType::index_type index_type;
        typedef typename LevelSetType::value_type value_type;

        //vector to store points that will be inserted into the lvlset
        typedef typename std::vector< std::pair< lvlset::vec<index_type,D>, value_type> > point_vector;
        point_vector points;

        //parse input file
        std::fstream json_input_stream;
        json_input_stream.open (p.geometry_files[0], std::fstream::in);

        if(json_input_stream.is_open()){

            json input = json::parse(json_input_stream);

            auto grid_properties = input.find("grid_size");

            if(grid_properties != input.end()){

                //create the grid

                int grid_min[D]={ };
                int grid_max[D]={ };

                for(int i=0; i<D; i++){

                    grid_min[i]=((*grid_properties)["dim_"+std::to_string(i)])[0];
                    grid_max[i]=((*grid_properties)["dim_"+std::to_string(i)])[1];
                }



                //TODO: move definition of boundary conditions somewhere else?
                //TODO: write them into the JSON file?
                //!Determine boundary conditions for level set domain
                lvlset::boundary_type bnc[D];
                for (int hh = 0; hh < D; ++hh) {
                    if (p.boundary_conditions[hh].min == bnc::PERIODIC_BOUNDARY &&
                        p.boundary_conditions[hh].max == bnc::PERIODIC_BOUNDARY) {
                        bnc[hh] = lvlset::PERIODIC_BOUNDARY;
                    } else if(p.boundary_conditions[hh].min == bnc::INFINITE_BOUNDARY &&
                            p.boundary_conditions[hh].max == bnc::INFINITE_BOUNDARY) {
                        bnc[hh] = lvlset::INFINITE_BOUNDARY;
                    } else if (p.boundary_conditions[hh].min == bnc::INFINITE_BOUNDARY) {
                        bnc[hh] = lvlset::NEG_INFINITE_BOUNDARY;
                    } else if (p.boundary_conditions[hh].max == bnc::INFINITE_BOUNDARY) {
                        bnc[hh] = lvlset::POS_INFINITE_BOUNDARY;
                    } else {
                        bnc[hh] = lvlset::SYMMETRIC_BOUNDARY;
                    }
                }

                //Set the level set GridProperties
                GridProperties = GridTraitsType(grid_min, grid_max, bnc, p.grid_delta);

                //Generate the grid with the GridProperties
                grid = lvlset::grid_type<GridTraitsType>(GridProperties);

                LevelSets.push_back(LevelSetType(grid));


            }else{

                std::cout << "error: grid_sizes has to be defined" << std::endl;

                //TODO: exit
            }

            auto geometries = input.find("geometry");

            

            if( geometries != input.end()){

                bool geomertry_found = false;

                //implement bool functions somewhere here!

                std::cout << (*geometries) << std::endl;

                for(auto geom: (*geometries)){

                    
                    if(geom["object"] == "circle"){
                        double mid[2] = {geom["center"][0], geom["center"][1]};
                        create_circle<D,LevelSetType, index_type, value_type>(geom["radius"], mid, LevelSets, points, p.grid_delta);
                    }

                    geomertry_found = true;


                }


                if(!geomertry_found){
                    std::cout << "error: geometry has to be defined" << std::endl;

                    //TODO: exit   
                }
                     

            }else{

                std::cout << "error: geometry has to be defined" << std::endl;

                //TODO: exit
            }


            

        }else{
            std::cout << "error: input geometry file not found" << std::endl;
        }










      

        std::sort(points.begin(),points.end());         //sort points lexicographically

        LevelSets.back().insert_points(points);  //initialize level set function


        LevelSets.back().prune();       //remove active grid point which have no opposite signed neighbor grid point
        LevelSets.back().segment();    //distribute points evenly across threads


        //typename LevelSetType::iterator it=LevelSets.begin();

        LevelSets.begin()->expand(2);

        write_explicit_surface_vtk(*(LevelSets.begin()),"test.vtk");

        write_explicit_levelset(*(LevelSets.begin()), "test_explicit.vtp");




    }

}









  #endif /*GEOMETRY_CUSTOM_H_*/

            /*
            lvlset::vec<index_type,D> pos_grid1;
            lvlset::vec<index_type,D> pos_grid2;

            double dist1 = 0.0;
            double dist2 = 0.0;


              if((abs(x_n+center[0])-abs(floor(x_n+center[0]))) < epsilon){
                //circle point lies directly on one grid point

                lvlset::vec<index_type,D> pos_grid3;

                pos_grid1[0]=floor(x_n+center[0])-1;
                //y_n koennte flasch sein
                pos_grid1[1]=floor(y_n+center[1]);

                dist1 = -0.75;
                
                pos_grid2[0]=floor(x_n+center[0])+1;
                //y_n koennte flasch sein
                pos_grid2[1]=floor(y_n+center[1]);

                dist2 = 0.75;

                pos_grid3[0]=floor(x_n+center[0]);
                //y_n koennte flasch sein
                pos_grid3[1]=floor(y_n+center[1]);

                points.push_back(std::make_pair(LevelSets.begin()->grid().global_indices_2_local_indices(pos_grid1), dist1));
                points.push_back(std::make_pair(LevelSets.begin()->grid().global_indices_2_local_indices(pos_grid2), dist2));
                points.push_back(std::make_pair(LevelSets.begin()->grid().global_indices_2_local_indices(pos_grid3), 0.0));



    
            }else{
                //circle point lies between 2 grid points
                pos_grid1[0]=floor(x_n+center[0]);
                //y_n koennte flasch sein
                pos_grid1[1]=floor(y_n+center[1]);

                dist1 = abs(floor(x_n+center[0])) - abs(x_n+center[0]);
                
                pos_grid2[0]=ceil(x_n+center[0]);
                //y_n koennte flasch sein
                pos_grid2[1]=floor(y_n+center[1]);

                dist2 = abs(ceil(x_n+center[0])) - abs(x_n+center[0]);

                points.push_back(std::make_pair(LevelSets.begin()->grid().global_indices_2_local_indices(pos_grid1), dist1));
                points.push_back(std::make_pair(LevelSets.begin()->grid().global_indices_2_local_indices(pos_grid2), dist2));

            }*/






  /*
       for(int i = 1;i<6;i++){

            for(int j = 1; j<6; j++){

                lvlset::vec<index_type,D> pos_grid;
                double Real_dist = 0.0;

                if(i == 1 || i == 5){

                    pos_grid[0]=i;
                    pos_grid[1]=j;

                    Real_dist = 0.5;

                    if(j==1 && i==1)
                        Real_dist = 0.75;

                    if(j==1 && i==5)
                        Real_dist = 0.75;
                    

                }


                if(i == 2 || i == 4){

                    pos_grid[0]=i;
                    pos_grid[1]=j;

                    if(j==1 || j== 5){
                        Real_dist = 0.5;
                    }else{
                        Real_dist = -0.5;
                    }
                }

                if(i == 3){

                    pos_grid[0]=i;
                    pos_grid[1]=j;

                    if(j==1 || j== 5){
                        Real_dist = 0.5;
                    }else if(j!=3){
                        Real_dist = -0.5;
                    }                 
                }

                if(Real_dist != 0)                 
                    points.push_back(std::make_pair(LevelSets.begin()->grid().global_indices_2_local_indices(pos_grid), Real_dist));
                
            }
       }*/
        
/*
            for(int i = 0;i<8;i++){
                            

                lvlset::vec<index_type,D> pos_grid;
                lvlset::vec<index_type,D> pos_grid2;


                pos_grid[0]=i;
                pos_grid[1]=2;


                //points.push_back(std::make_pair(pos_grid, 0.5));
                points.push_back(std::make_pair(LevelSets.begin()->grid().global_indices_2_local_indices(pos_grid), 0.5));
                

                pos_grid2[0]=i;
                pos_grid2[1]=3;

                //points.push_back(std::make_pair(pos_grid2, -0.5));
                points.push_back(std::make_pair(LevelSets.begin()->grid().global_indices_2_local_indices(pos_grid2), -0.5));


                lvlset::vec<index_type,D> pos_grid3;
                lvlset::vec<index_type,D> pos_grid4;


                pos_grid3[0]=i;
                pos_grid3[1]=6;


                //points.push_back(std::make_pair(pos_grid, 0.5));
                points.push_back(std::make_pair(LevelSets.begin()->grid().global_indices_2_local_indices(pos_grid3), -0.5));
                

                pos_grid4[0]=i;
                pos_grid4[1]=7;

                //points.push_back(std::make_pair(pos_grid2, -0.5));
                points.push_back(std::make_pair(LevelSets.begin()->grid().global_indices_2_local_indices(pos_grid4), 0.5));


            }
*/
       