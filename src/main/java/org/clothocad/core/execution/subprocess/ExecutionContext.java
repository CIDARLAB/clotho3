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
        final Object return_value;
        final Map<String, Object> out_value = new HashMap<>();
        final String in_value_name = (String) in_value.get("name");
        final List in_value_args = (List) in_value.get("args");
        try {
            return_value = APIRelayer.relay(api, in_value_name, in_value_args);
        } catch (final Exception e) {
            out_value.put("type", "api_error");
            out_value.put("message", e.getMessage());
            sendValue(out_value);
            return;
        }
        out_value.put("type", "api");
        out_value.put("return", return_value);
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
        final byte[] bytes = JSONUtil.encodeUTF8(value);
        try {
            output.write(bytes);
            output.write(0);
            output.flush();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}
