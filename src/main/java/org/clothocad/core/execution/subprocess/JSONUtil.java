package org.clothocad.core.execution.subprocess;

import static java.nio.charset.StandardCharsets.UTF_8;

import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.ObjectMapper;
import java.io.IOException;

class JSONUtil {
    private static final ObjectMapper mapper = new ObjectMapper();

    static Object fromBytes(final byte[] bytes) {
        try {
            return mapper.readValue(new String(bytes, UTF_8), Object.class);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    static byte[] toBytes(final Object value) {
        try {
            return mapper.writeValueAsBytes(value);
        } catch (JsonProcessingException e) {
            throw new RuntimeException(e);
        }
    }
}
