main: graph neighborhood unionFind edge mymain.cpp
	g++ -std=c++11 -o mymain mymain.cpp *.o `pkg-config gtkmm-3.0 --cflags --libs`
cluster-run: graph MyArea polyfit cluster-run.cpp
	g++ -std=c++11 -o cluster-run cluster-run.cpp *.o `pkg-config gtkmm-3.0 --cflags --libs`
treatdata: graph treatdata.cpp
	g++ -std=c++11 -o treatdata treatdata.cpp *.o `pkg-config gtkmm-3.0 --cflags --libs`
graph: edge neighborhood unionFind graph.cpp graph.hpp
	g++ -c graph.cpp
cloud: point cloud.cpp cloud.hpp
	g++ -c cloud.cpp
point: point.cpp point.hpp
	g++ -c point.cpp
edge: cloud edge.cpp edge.hpp
	g++ -c edge.cpp
unionFind: unionFind.cpp unionFind.hpp
	g++ -c unionFind.cpp
neighborhood: edge neighborhood.cpp neighborhood.hpp
	g++ -c neighborhood.cpp
MyArea: cloud MyArea.cpp MyArea.hpp
	g++ -std=c++11 -c MyArea.cpp `pkg-config gtkmm-3.0 --cflags --libs`
polyregression: polyregression.cpp
	g++ -std=c++17 -o polyregression polyregression.cpp
polyfit: polyfit.c polyfit.h
	g++ -c polyfit.c

clean:
	rm -vf mymain *.o *~
