import numpy as np
from numpy import uint8, uint32, uint16, uint64
from collections import Counter
import gzip

from enum import IntEnum

from PIL import Image

bits_op = 62
bits_len = 55
addr_mask = uint64((1 << bits_len) - 1)
len_mask = uint64(((1 << (bits_op - bits_len)) - 1) << (bits_len - 1))
op_mask = uint64(((1 << (64 - bits_op)) - 1) << (bits_op - 1))

opBytesToIntsDict = {
    b'I ': 0 << bits_op,
    b' L': 1 << bits_op,
    b' S': 2 << bits_op,
    b' M': 3 << bits_op}
def decode_line(line):
    op = opBytesToIntsDict.get(line[0:2], None)
    if op is None:
        return None

    comma_pos = line.find(b',')
    addr = int(line[3:comma_pos], base=16) | op
    #assert((addr >> 54) == 0)

    addr |= int(line[comma_pos+1:]) << bits_len
    #assert(oplen <= 127)

    return np.uint64(addr)

class Trace:

    def __init__(self):
        self.data=None
        self.address_counter=None

    def __str__(self):
        return str(len(self.data)) + " traces\n"

    def read_gzip_trace(self, filename, max_traces = float('inf')):
        # temporary list of decoded lines.
        decoded_lines = []

        zip_stream = gzip.open(filename)
        for trace_nr, stream_line in enumerate(zip_stream):
            if trace_nr >= max_traces:
                break
            line = decode_line(stream_line)
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

    def used_addresses(self):
        self.address_counter = Counter(addressesGenerator(self.data))
        addresses = np.array(list(self.address_counter.keys()))
        addresses.sort()
        return addresses

def trace_to_addr
def addressesGenerator(data):
    """ Generates a sequence of numbers [addr, addr+1, ..., addr+len-1]
    for each pair (addr, len) from data.
    """
    for line in data:
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
    img.fill(0xFF000000)

    for a in addressesGenerator(traces.data):
        x, y = range_map.coordinates(a)
        img[y, x] = 0xFFFFFFFF
    return img

    ndarray_to_rgba_image(img)


def ndarray_to_rgba_image(data, filename='image.png'):
    h, w = data.shape
    pilImage = Image.frombuffer('RGBA', (w,h), data, 'raw', 'RGBA', 0, 1)
    pilImage.save(filename)
    return pilImage

def prevPower2Int(x):
    return 2 << int(np.ceil(np.log(x)/np.log(2)) - 2)

def nextPower2Int(x):
    return 2 << int(np.ceil(np.log(x)/np.log(2)) - 1)

def instructions_to_img(counts, w, h=None):
    if h is None:
        h = nextPower2Int(len(counts)/w)

    img = np.empty((w,h), np.uint32)
    img.shape = h, w
    img.fill(0xFF000000)

    max_count = counts[0:w*h].max()
    scale_factor = float(0xFFFFFF)/max_count
    for addr, count in enumerate(counts):
        if addr > h*w:
            break
        if count == 0:
            continue

        x, y = divmod(addr, w)
        img[x,y] += int(count*scale_factor)
        print(addr, x, y, count, hex(img[x,y]))

    pilImage = Image.frombuffer('RGBA', (w,h), img, 'raw', 'RGBA', 0, 1)
    pilImage.save('i.png')
    return pilImage

def count_instr_addr(time_line):
    addr_range = time_line.iaddr_range

    addr_count = np.empty(addr_range[1] - addr_range[0], np.uint32)

    offset = addr_range[0]
    for line in time_line.data:
        addr = line[0][1]-offset
        ilen = line[0][2]
        addr_count[addr:addr+ilen] += 1

    return addr_count, offset

def heap_to_img(counts, w, h=None):
    if h is None:
        h = nextPower2Int(len(counts)/w)

    img = np.empty((w,h), np.uint32)
    img.shape = h, w
    img.fill(0xFF000000)

    max_count = counts[0:w*h].max()
    scale_factor = float(0xFFFFFF)/max_count
    for addr, count in enumerate(counts):
        if addr > h*w:
            break
        if count == 0:
            continue

        x, y = divmod(addr, w)
        img[x,y] += int(count*scale_factor)
        print(addr, x, y, count, hex(img[x,y]))

    pilImage = Image.frombuffer('RGBA', (w,h), img, 'raw', 'RGBA', 0, 1)
    pilImage.save('i.png')
    return pilImage

def count_heap_addr(time_line):
    addr_range = time_line.haddr_range

    addr_count = np.empty(addr_range[1] - addr_range[0], np.uint32)

    offset = addr_range[0]
    for line in time_line.data:
        for i in range(1, len(line)):
            addr = line[i][1]-offset
            ilen = line[i][2]
            addr_count[addr:addr+ilen] += 1

    return addr_count, offset

def count_stack_addr(time_line):
    addr_range = time_line.saddr_range

    addr_count = np.empty(addr_range[1] - addr_range[0], np.uint32)

    offset = addr_range[0]
    for line in time_line.data:
        for i in range(1, len(line)):
            addr = line[i][1]-offset
            ilen = line[i][2]
            addr_count[addr:addr+ilen] += 1

    return addr_count, offset

def time_view(time_line, w, h, offset=None):

    if offset is None:
        offset = time_line.saddr_range[0]

    img = np.empty((w,h), np.uint32)
    img.shape = h, w
    img.fill(0xFF000000)

    t = 0
    h_max = 0
    for time, line in enumerate(time_line.data):
        if t >= w:
            break

        t_inc = False
        for i in range(len(line)):
            addr = (line[i][1] - offset)//8
            ilen = line[i][2]
            print(t, addr, time, hex(addr), ilen)

            if addr > h:
                continue

            img[addr, t] = 0xFFFFFFFF

            h_max = max(h_max, addr)
            t_inc = True

        if t_inc:
            t += 1

    #img = np.resize(img, (h_max, t))
    pilImage = Image.frombuffer('RGBA', (w,h), img, 'raw', 'RGBA', 0, 1)
    pilImage.save('i.png')
    return pilImage

def run(filename):
    t = Trace()
    l = t.read_gzip_trace(filename)
