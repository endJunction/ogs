import numpy as np
from numpy import uint8, uint32, uint16, uint64
from collections import Counter
import gzip

from enum import IntEnum

from PIL import Image

bits_op = 61
bits_len = 54
addr_mask = uint64((1 << bits_len) - 1)
len_mask = uint64(((1 << (bits_op - bits_len)) - 1) << bits_len)
op_mask = uint64(((1 << (64 - bits_op)) - 1) << bits_op)

opBytesToIntsDict = {
    b'I ': 0b001 << bits_op,
    b' L': 0b010 << bits_op,
    b' S': 0b100 << bits_op,
    b' M': 0b110 << bits_op}

def encode_trace(line):
    op = opBytesToIntsDict.get(line[0:2], None)
    if op is None:
        return None

    comma_pos = line.find(b',')
    addr = int(line[3:comma_pos], base=16) | op
    #assert((addr >> 54) == 0)

    addr |= int(line[comma_pos+1:]) << bits_len
    #assert(oplen <= 127)

    return np.uint64(addr)

def decode_trace(line):
    addr = int(line & addr_mask)
    oplen = int(line & len_mask) >> bits_len
    op = int(line & op_mask) >> bits_op

    return (op, addr, oplen)

class Trace:

    def __init__(self):
        self.data=None

    def __str__(self):
        return str(len(self.data)) + " traces\n"

    def __iter__(self):
        """ Iterator object over traces. """
        return (line for line in self.data)

    def read_gzip_trace(self, filename, max_traces = float('inf')):
        # temporary list of decoded lines.
        decoded_lines = []

        zip_stream = gzip.open(filename)
        for trace_nr, stream_line in enumerate(zip_stream):
            if trace_nr >= max_traces:
                break
            line = encode_trace(stream_line)
            # skip 'undecodable' lines, those _not_ starting with I, L, S, or M
            if line is None:
                continue
            decoded_lines.append(line)

        zip_stream.close()

        if self.data is None:
            self.data = np.array(decoded_lines)
        else:
            self.data = np.append(self.data, decoded_lines)

        return self.data

    def store_binary(self, filename='traces.bin'):
        self.data.tofile(filename)

    def load_binary(self, filename='traces.bin'):
        self.data = np.fromfile(filename, dtype=np.uint64)

    def used_addresses(self, operator_mask = 0b111):
        address_counter = Counter(addressesGenerator(self.data, operator_mask))
        addresses = np.array(list(address_counter.keys()))
        addresses.sort()
        return addresses, address_counter

def addressesGenerator(data, operator_mask=0b111):
    """ Generates a sequence of numbers [addr, addr+1, ..., addr+len-1]
    """
    for line in data:
        if (int(line) >> bits_op) & operator_mask == 0:
            continue
        addr = int(line & addr_mask)
        oplen = int(line & len_mask) >> bits_len
        for i in range(addr, addr + oplen):
            yield i

def find_nonzero_ranges(data, size=256, padding=8):
    """ Works on sorted data only. """

    def set_range_end(value):
        return value  | ((1 << padding) - 1)

    def set_start_end_range(value):
        s = value >> padding << padding
        e = set_range_end(value)
        return s, e

    range_start, range_end = set_start_end_range(data[0])

    ranges = []
    def append_range(s, e):
        if len(ranges) == 0:
            ranges.append((s, e))
            return

        last_range = ranges[-1]
        if last_range[1] + 1 == s:
            # extend the last range
            ranges[-1] = (ranges[-1][0], e)
        else:
            ranges.append((s, e))

    for i in data:
        if i <= range_end + size:
            range_end = set_range_end(i)
            continue

        # Close range and prepare next
        append_range(range_start, range_end)
        range_start, range_end = set_start_end_range(i)

    # Close the last range
    append_range(range_start, range_end)

    return ranges

class RangeMapping:
    def __init__(self, ranges, row_length = 128):
        # ranges are result of find_nonzero_ranges(xx)
        self.ranges = ranges
        self.start_addresses = np.array([r[0] for r in ranges])
        self.lengths = np.array([r[1] - r[0] for r in ranges])

        # each range has a start row in the image
        self.row_length = row_length
        self.row_index = np.append(
                [0], np.cumsum((self.lengths+1) // row_length))


    def __str__(self):
        result = ""
        if len(self.ranges) <= 20:
            for r in self.ranges:
                result += hex(r[0]) + " " + hex(r[1]) + "\n"
            return result

        for i in range(10):
            result += hex(self.ranges[i][0]) + " " + hex(self.ranges[i][1]) + "\n"
        result += "...\n"
        for i in range(-10, 0):
            result += hex(self.ranges[i][0]) + " " + hex(self.ranges[i][1]) + "\n"

        return result

    def coordinates(self, addr):
        """ Returns a pair of coordinates in an 2D array of width row_length.
        """

        # Find a range index including the given address.
        index = self.start_addresses.searchsorted(addr, side='right')
        if index == 0:
            raise IndexError("too small address requested")

        # Decrement for correct start.
        index -= 1

        assert(addr >= self.ranges[index][0])
        if addr > self.ranges[index][1]:
            raise IndexError("found range does not include the requested address.")

        y, x = divmod(addr - self.start_addresses[index], self.row_length)
        y += self.row_index[index]
        return (x, y)

def create_address_image(traces, range_map):
    w = range_map.row_length
    h = range_map.row_index[-1]
    img = np.empty((w,h), np.uint32)
    img.shape = h, w
    # color format is AABBGGRR
    img.fill(0xFF000000)

    for t in traces.data:
        op, a, op_len = decode_trace(t)

        try:
            x, y = range_map.coordinates(a)
        except IndexError:
            continue

        if op == 0b001:
            img[y, x:x+op_len] += 0xA4DDAB
        elif op == 0b010:
            img[y, x:x+op_len] += 0xBA832B
        elif op == 0b100:
            img[y, x:x+op_len] += 0x61AEFD
        elif op == 0b110:
            img[y, x:x+op_len] += 0x1C19D7

    return img


def ndarray_to_rgba_image(data, filename='image.png'):
    h, w = data.shape
    pilImage = Image.frombuffer('RGBA', (w,h), data, 'raw', 'RGBA', 0, 1)
    pilImage.save(filename)
    return pilImage

