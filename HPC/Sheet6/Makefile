INC=-I./inc
.PHONY: build
build: ./bin/relaxation

.PHONY: clean
clean:
	rm -f ./bin/*
	
.PHONY: rebuild
rebuild: clean build

./bin/relaxation: ./src/main.c
	mpicxx "-Wall" -O3 -o $@ $^ -lm $(INC)
