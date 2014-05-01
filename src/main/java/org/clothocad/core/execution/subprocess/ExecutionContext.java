package org.clothocad.core.execution.subprocess;

import java.io.IOException;
import java.io.OutputStream;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.clothocad.core.communication.ServerSideAPI;

class ExecutionContext {
    private final ServerSideAPI api;
    private final JSONStreamReader reader;
    private final OutputStream output;
    private final String code;
    private final List<Object> args;

    ExecutionContext(final ServerSideAPI api,
                     final JSONStreamReader reader,
                     final OutputStream output,
                     final String code,
                     final List<Object> args) {
        this.api = api;
        this.reader = reader;
        this.output = output;
        this.code = code;
        this.args = args;
    }

    /** Communicate with subprocess */
    Object
    start() {
        sendFunctionDefcall();

        while (true) {
            final Map value;
            try {
                value = (Map) reader.waitValue();
            } catch (InterruptedException e) {
                throw new RuntimeException(e);
            }
            final String valueType = (String) value.get("type");
            if ("func".equals(valueType)) {
                return handleFunctionReturn(value);
            } else if ("api".equals(valueType)) {
                handleAPICall(value);
            } else {
                /* unrecognized message */
                throw new RuntimeException();
            }
        }
    }

    private Object
    handleFunctionReturn(final Map value) {
        if (!value.containsKey("return"))
            throw new RuntimeException();
        return value.get("return");
    }

    private void
    handleAPICall(final Map in_value) {
        final Map<String, Object> out_value = new HashMap<>();
        out_value.put("type", "api");
        out_value.put(
            "return",
            APIRelayer.relay(
                api,
                (String) in_value.get("name"),
                (List) in_value.get("args")
            )
        );
        sendValue(out_value);
    }

    private void
    sendFunctionDefcall() {
        final Map<String, Object> val = new HashMap<>();
        val.put("type", "func");
        val.put("code", code);
        val.put("args", args);
        sendValue(val);
    }

    private void
    sendValue(final Object value) {
        final byte[] bytes = JSONUtil.toBytes(value);
        try {
            output.write(bytes);
            output.write(0);
            output.flush();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}
