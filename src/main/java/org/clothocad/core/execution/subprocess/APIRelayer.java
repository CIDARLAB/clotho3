package org.clothocad.core.execution.subprocess;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.clothocad.core.communication.ServerSideAPI;

class APIRelayer {
    private static final Map<String, Handler> HANDLERS = new HashMap<>();
    static {
        HANDLERS.put("get", new GetHandler());
        HANDLERS.put("set", new SetHandler());
        HANDLERS.put("run", new RunHandler());
    }

    static Object relay(final ServerSideAPI api,
                        final String name,
                        final List<Object> args) {
        final Handler h = HANDLERS.get(name);
        if (h == null)
            throw new RuntimeException(
                String.format("no handler for %s", name)
            );
        return h.handle(api, args);
    }

    private static interface Handler {
        Object handle(ServerSideAPI api, List<Object> args);
    }

    private static class GetHandler implements Handler {
        @Override public Object
        handle(final ServerSideAPI api, final List<Object> args) {
            if (args.size() != 1)
                throw new RuntimeException();
            return api.get(args.get(0));
        }
    }

    private static class SetHandler implements Handler {
        @Override public Object
        handle(final ServerSideAPI api, final List<Object> args) {
            if (args.size() != 1)
                throw new RuntimeException();
            return api.set((Map) args.get(0));
        }
    }

    private static class RunHandler implements Handler {
        @Override public Object
        handle(final ServerSideAPI api, final List<Object> args) {
            if (args.size() != 2)
                throw new RuntimeException();
            return api.run2((String) args.get(0), (List) args.get(1));
        }
    }
}
