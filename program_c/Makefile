.PHONY: all clean run30

CC = gcc
CFLAGS = -Wall -Wextra -O2

SRC_DIR = ./src
BIN_DIR = ./bin
SRC = $(wildcard $(SRC_DIR)/*.c)
OBJ = $(SRC:.c=.o)

TARGET = $(BIN_DIR)/pca_program

all: $(TARGET)

$(TARGET): $(OBJ)
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJ) -lm

$(SRC_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	@rm -f $(TARGET) $(SRC_DIR)/*.o *.csv

run30: 
	@./bin/pca_program ./data/ref/01_ref30.txt ./data/01_ind30

run300: 
	@./bin/pca_program ./data/ref/02_ref300.txt ./data/02_ind300

run3k: 
	@./bin/pca_program ./data/ref/03_ref3k.txt ./data/03_ind3k

run30k: 
	@./bin/pca_program ./data/ref/04_ref30k.txt ./data/04_ind30k

run300k: 
	@./bin/pca_program ./data/ref/05_ref300k.txt ./data/05_ind300k

run3M: 
	@./bin/pca_program ./data/ref/06_ref3M.txt ./data/06_ind3M

run30M: 
	@./bin/pca_program ./data/ref/07_ref30M.txt ./data/07_ind30M

run300M:
	@./bin/pca_program ./data/ref/08_ref300M.txt ./data/08_ind300M