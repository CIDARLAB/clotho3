package org.clothocad.core.execution.subprocess;

import static java.nio.charset.StandardCharsets.UTF_8;

import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.ObjectMapper;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.clothocad.core.util.JSON;

class JSONUtil {
    private static final ObjectMapper mapper = new ObjectMapper();

    //clotho.get("org.andersonlab.py_fetchRegistry");
    //clotho.run("org.andersonlab.py_fetchRegistry", ["BBa_b0015"]);
    static Object decodeUTF8(final byte[] bytes) {
        try {
            String msg = new String(bytes, UTF_8);
            LinkedHashMap out = (LinkedHashMap) mapper.readValue(msg, Object.class);
            //If the out is a wrapped JSON string, then unwrap it
            try {
                Object ret = out.get("return");
                List obs = (List) ret;
                String jsonstr = (String) obs.get(0);
                

                Object json = JSON.deserializeObjectToMap(jsonstr);
                System.out.println(json);
                
                //Repackage
                obs.remove(0);
                obs.add(0, json);
                out.put("return", obs);
                
                
            } catch(Exception err) {
                
            }
            return out;
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
