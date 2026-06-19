CPPFLAGS ?= -Isrc
CFLAGS ?= -std=c99 -O3 -march=native -D_GNU_SOURCE
DEBUG_CFLAGS ?= -std=c99 -O0 -ggdb3 -fsanitize=address -march=native -D_GNU_SOURCE -DDEBUG=1 -Wall -Wextra -Wformat=2 -Wshadow
LDFLAGS ?=
DEBUG_LDFLAGS ?= -fsanitize=address
LDLIBS += -lm -lasound -lfftw3 -lfftw3_threads -pthread

TARGET = gcFT8
DEBUG_TARGET = gcFT8-debug
RELEASE_OBJDIR = build/release
DEBUG_OBJDIR = build/debug

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

RELEASE_OBJS = $(SRCS:src/%.c=$(RELEASE_OBJDIR)/%.o)
DEBUG_OBJS = $(SRCS:src/%.c=$(DEBUG_OBJDIR)/%.o)

.PHONY: run_tests all debug clean install

all: $(TARGET)

$(TARGET): $(RELEASE_OBJS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS)

debug: $(DEBUG_TARGET)

$(DEBUG_TARGET): $(DEBUG_OBJS)
	$(CC) $(LDFLAGS) $(DEBUG_LDFLAGS) -o $@ $^ $(LDLIBS)

$(RELEASE_OBJDIR)/%.o: src/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<

$(DEBUG_OBJDIR)/%.o: src/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) $(DEBUG_CFLAGS) -c -o $@ $<

clean:
	rm -rf build *.o ft8/*.o common/*.o fft/*.o serial/*.o hash/*.o $(TARGET) $(DEBUG_TARGET) libftx.a

install: all
	$(AR) rc libftx.a $(RELEASE_OBJDIR)/protocol/ftx/constants.o $(RELEASE_OBJDIR)/protocol/ftx/encode.o $(RELEASE_OBJDIR)/protocol/ft8/ft8_encode.o $(RELEASE_OBJDIR)/protocol/ft4/ft4_encode.o $(RELEASE_OBJDIR)/protocol/ft2/ft2_encode.o $(RELEASE_OBJDIR)/protocol/ftx/message.o $(RELEASE_OBJDIR)/protocol/ftx/text.o
	install libftx.a /usr/lib/libftx.a
