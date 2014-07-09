package org.clothocad.core.execution.subprocess;

import static java.nio.charset.StandardCharsets.UTF_8;

import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.ObjectMapper;
import java.io.IOException;

class JSONUtil {
    private static final ObjectMapper mapper = new ObjectMapper();

    static Object decodeUTF8(final byte[] bytes) {
        try {
            String jsonstr = new String(bytes, UTF_8);
            return mapper.readValue(jsonstr, Object.class);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    static byte[] encodeUTF8(final Object value) {
        try {
            return mapper.writeValueAsBytes(value);
        } catch (JsonProcessingException e) {
            throw new RuntimeException(e);
        }
    }
}
