package org.clothocad.core.util;

import java.util.Arrays;

public class ByteArray {
    private byte[] array;
    private int size;

    public ByteArray() {
        array = new byte[16];
        size = 0;
    }

    public byte[] getArray() {
        return Arrays.copyOf(array, size);
    }

    public void clear() {
        if (16 <= size && 2 * size < array.length)
            array = new byte[size];
        size = 0;
    }

    public void add(final byte b) {
        ensureSize();
        array[size] = b;
        size++;
        return;
    }

    private void ensureSize() {
        while (size >= array.length)
            array = Arrays.copyOf(array, 1 + 2 * array.length);
    }
}
