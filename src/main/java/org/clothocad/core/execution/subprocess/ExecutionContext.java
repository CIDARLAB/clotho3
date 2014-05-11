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

    private void
    sendFunctionDefcall() {
        final Map<String, Object> value = new HashMap<>();
        value.put("type", "func");
        value.put("code", code);
        value.put("args", args);
        sendValue(value);
    }

    private Object
    handleFunctionReturn(final Map value) {
        if (!value.containsKey("return"))
            throw new RuntimeException();
        return value.get("return");
    }

    private void
    handleAPICall(final Map value) {
        final Map<String, Object> reply = new HashMap<>();
        final APIRelayer.Callback cb = new APIRelayer.Callback() {
            @Override public void onSuccess(final Object ret) {
                reply.put("type", "api");
                reply.put("return", ret);
            }
            @Override public void onFail(final String message) {
                reply.put("type", "api_error");
                reply.put("message", message);
            }
        };
        APIRelayer.relay(api,
                         (String) value.get("name"),
                         (List) value.get("args"),
                         cb);
        sendValue(reply);
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
