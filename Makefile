CFLAGS = -std=c99 -O3 -ggdb3 -fsanitize=address -march=native -D_GNU_SOURCE
CPPFLAGS = -Isrc
LDFLAGS = -lm -fsanitize=address 
LDLIBS += -lm -lasound -lfftw3 -lfftw3_threads -pthread

TARGETS = gcFT8
OBJDIR = build

SRCS = \
	src/main.c \
	src/dsp/monitor.c \
	src/dsp/gfsk.c \
	src/protocol/ft2/ft2_waveform.c \
	src/protocol/ft8/ft8_encode.c \
	src/protocol/ft4/ft4_encode.c \
	src/protocol/ft2/ft2_encode.c \
	src/vendor/kissfft/kiss_fftr.c \
	src/vendor/kissfft/kiss_fft.c \
	src/protocol/ftx/decode.c \
	src/protocol/ftx/message.c \
	src/protocol/ftx/encode.c \
	src/protocol/ftx/crc.c \
	src/protocol/ftx/ldpc.c \
	src/protocol/ftx/text.c \
	src/protocol/ftx/constants.c \
	src/vendor/cssl/cssl.c \
	src/util/hash_table.c

OBJS = $(SRCS:src/%.c=$(OBJDIR)/%.o)

.PHONY: run_tests all clean install

all: $(TARGETS)

gcFT8: $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(OBJDIR)/%.o: src/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<

clean:
	rm -rf $(OBJDIR) *.o ft8/*.o common/*.o fft/*.o serial/*.o hash/*.o $(TARGETS) libftx.a
install: all
	$(AR) rc libftx.a $(OBJDIR)/protocol/ftx/constants.o $(OBJDIR)/protocol/ftx/encode.o $(OBJDIR)/protocol/ft8/ft8_encode.o $(OBJDIR)/protocol/ft4/ft4_encode.o $(OBJDIR)/protocol/ft2/ft2_encode.o $(OBJDIR)/protocol/ftx/message.o $(OBJDIR)/protocol/ftx/text.o
	install libftx.a /usr/lib/libftx.a
