# Compiler Mode Selection
# Usage:
#   make OPENACC=1    # for OpenACC GPU mode
#   make DEBUG=1      # to include debug symbols
#   make              # for default CPU mode
# ----- Mode Selection -----

ifeq ($(OPENACC),1)
    CC = nvc
    MODE = openacc
    CFLAGS = -acc -Minfo=accel -O3 -Wall -Wpointer-arith
    LDFLAGS = -lm
    GPU_MSG = "Compiling with NVC (OpenACC GPU mode)"

else ifeq ($(ACC),1)
    CC = gcc
    MODE = acc
    CFLAGS = -fopenacc -fopenmp -O3 -Wall -Wpointer-arith -march=native
    LDFLAGS = -lm
    GPU_MSG = "Compiling with GCC OpenACC"

else
    CC = gcc
    MODE = cpu
    CFLAGS = -O3 -march=native -fopenmp -Wall -Wpointer-arith
    LDFLAGS = -lm
    GPU_MSG = "Compiling with GCC (CPU mode)"
endif

ifeq ($(DEBUG),1)
    CFLAGS += -g
endif

# Files

SRC_DIR = src/c_header_files
SRC = $(wildcard $(SRC_DIR)/*.c)

INIT = init.c
INIT_SRC := $(wildcard $(INIT))   # Will be empty if file doesn't exist

MAIN = mg_NS_solver.c
TARGET = a.out
MODE_FILE = .build_mode


# Build Rules

all: check_mode $(TARGET)

$(TARGET): $(SRC) $(MAIN) $(INIT_SRC)
	@echo $(GPU_MSG)
	@if [ -z "$(INIT_SRC)" ]; then \
	    echo "WARNING: $(INIT) not found. Using zero-initialisation fallback."; \
	fi
	$(CC) $(SRC) $(INIT_SRC) $(MAIN) $(CFLAGS) -o $(TARGET) $(LDFLAGS)

check_mode:
	@if [ -f $(MODE_FILE) ] && [ "`cat $(MODE_FILE)`" != "$(MODE)" ]; then \
	    echo "Build mode changed → full rebuild"; \
	    rm -f $(TARGET); \
	fi
	@echo $(MODE) > $(MODE_FILE)

clean:
	rm -f $(TARGET) $(MODE_FILE)