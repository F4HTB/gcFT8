# CC = g++ 

CFLAGS = -O3 -ggdb3 -fsanitize=address -march=native
CPPFLAGS = -I. 
LDFLAGS = -lm -fsanitize=address 
LDLIBS += -lm -lasound -lfftw3 -lfftw3_threads -pthread

TARGETS = clFT8

.PHONY: run_tests all clean

all: $(TARGETS)

clFT8: clFT8.o fft/kiss_fftr.o fft/kiss_fft.o ft8/decode.o ft8/pack.o ft8/encode.o ft8/crc.o ft8/ldpc.o ft8/unpack.o ft8/text.o ft8/constants.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

clean:
	rm -f *.o ft8/*.o common/*.o fft/*.o $(TARGETS)
install:
	$(AR) rc libft8.a ft8/constants.o ft8/encode.o ft8/pack.o ft8/text.o
	install libft8.a /usr/lib/libft8.a
