CC:=g++
exe:=main.cpp struct_vars.h init_param.h LRstate.h LRstate.cpp \
	C_to_P.h C_to_P.cpp cfl.h cfl.cpp fluxes.h fluxes.cpp HLLD.h HLLD.cpp \
	update_U.h update_U.cpp BC.h BC.cpp write_log.h write_log.cpp
obj:=main.o

all:
	$(CC) $(exe) -o $(obj)
	
clean:
	rm -rf $(obj)
