

ifneq (,$(findstring dhvi,$(shell uname -n)))#links and libs for duke cluster                                                                                                                                     
	CFLAGS_LINK=-ggdb -O3 -std=c++11 -I/datacommons/dhvi/scripts/lib/boost/boost_1_70_0 -L/datacommons/dhvi/scripts/lib/boost/boost_1_70_0/stage/lib -lboost_regex -pthread
	CFLAGS=-ggdb -w -O3 -std=c++11 -I/datacommons/dhvi/scripts/lib/boost/boost_1_70_0
else ifneq (,$(findstring Darwin,$(shell uname -s)))#links and libs for osx 
	CFLAGS_LINK=-ggdb -O3 -std=c++11 -I/usr/local/lib -L/usr/local/lib -lboost_regex -pthread
	CFLAGS=-ggdb -w -O3 -std=c++11 -I/usr/local/lib
else #links for linux   
	CFLAGS_LINK=-ggdb -O3 -std=c++11 -I/usr/lib/x86_64-linux-gnu/ -L/usr/lib/x86_64-linux-gnu/ -lboost_regex -pthread
	CFLAGS=-ggdb -w -O3 -std=c++11 -I/usr/lib/x86_64
endif

COMPILER=g++

translater: utilities.o main.o
	$(COMPILER) $(CFLAGS_LINK) utilities.o main.o -o translater

main.o: main.cpp utilities.hpp
	$(COMPILER) -c $(CFLAGS) main.cpp

utilities.o : utilities.cpp utilities.hpp 
	$(COMPILER) -c $(CFLAGS) utilities.cpp

clean: 
	rm -f translater *.o *~
