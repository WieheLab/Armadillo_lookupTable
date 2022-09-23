

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

positional_aa_matrix_maker: positional_aa_matrix_maker.cpp
	$(COMPILER) $(CFLAGS) positional_aa_matrix_maker.cpp -o positional_aa_matrix_maker

clean: 
	rm -f positional_aa_matrix_maker *~
