#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cfloat> // for DBL_MAX
#include <fstream>

#include "graph.hpp"
#include "cloud.hpp"
#include "MyArea.hpp"

using std::cout;
using std::endl;
using std::ifstream;

int main(int argc, char *argv[])
{
    const int d = 2; //dimension of the data file
    const int nmax = 8001;
    int best_k = 1;
    double intrac_d, interc_d, best_ratio = __DBL_MAX__;

    if (argc < 5)
    {
        cout << "./cluster-run \"input_file\" \"max_of_clusters\" \"output_file\" \"method\"" << endl;
        cout << "Method: 1 - This method consists in removing the biggest edges until it reachs to create another cluster. It aims to minimize the relation (Maximum intracluster distance)/((Minimum intercluster distance)*N_cluster^2). The N_clusters^2 was empiricly chosen, it seems to do a good job.\n\nMethod: 2 - This method is a implementation of the Maximum Standard Deviation Reduction Clustering Algorithm (MSDR), that can be found in the article: https://ieeexplore.ieee.org/document/4031882\n";
        return -1;
    }
    int max_k = atoi(argv[2]);
    int method = atoi(argv[4]);
    const char *infile = argv[1];
    const char *outfile = argv[3];
    // construct point cloud
    cloud c(d, nmax, 1);

    // open input data file
    ifstream is(infile);
    assert(is.is_open());
    
    //open output data file
    ofstream output(outfile);
    output.clear();

    //read and store data file into the cloud and close the file
    c.load(is);
    //max_k = 30;
    edge::set_cloud(&c);
    //Task 5 and 6
    graph g = graph(c);
    if(method == 1)
        std::list<std::list<edge>> *ret = g.getBestk(&c, max_k);  
    else 
        std::list<std::list<edge>> *ret = g.MSDR(&c, max_k);
    c.save(output);
    //End task 5 and 6
   
    // launch graphical interface
    //Glib::RefPtr<Gtk::Application> app = Gtk::Application::create(argc, argv, "inf442.td3");
    int fake_argc = 1;
    Glib::RefPtr<Gtk::Application> app = Gtk::Application::create(fake_argc, argv, "inf442.td3");
    Gtk::Window win;
    win.set_title("PI");
    win.set_default_size(400, 400);

    MyArea area(&c);
    win.add(area);
    area.show();
    app->run(win);

    if(!output.is_open()){
        output.close();
    }
    if(!is.is_open()){
        is.close();
    }
    return 1;
}
