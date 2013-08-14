/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.layers.execution;

import java.io.IOException;
import java.io.Reader;
import java.io.StringReader;
import java.util.HashMap;
import javax.script.AbstractScriptEngine;
import javax.script.Bindings;
import javax.script.ScriptContext;
import javax.script.ScriptEngineFactory;
import javax.script.ScriptException;
import javax.script.SimpleBindings;
import org.mozilla.javascript.*;

/**
 *
 * @author spaige
 */
class JavaScriptEngine extends AbstractScriptEngine implements HackEngine {

    private final HashMap<Object, Object> indexedProps;
    private final ScriptableObject stdObjects;

    public JavaScriptEngine() {
        indexedProps = new HashMap<>();
        Context cx = Context.enter();
        stdObjects = cx.initStandardObjects();
        Context.exit();
    }

    @Override
    public Object eval(Reader reader, ScriptContext ctxt) throws ScriptException {
        Object ret;

        Context cx = Context.enter();
        try {
            Scriptable scope = getRuntimeScope(ctxt);
            String filename = "<Unknown source>";

            ret = cx.evaluateReader(scope, reader, filename, 1, null);
        } catch (RhinoException re) {
            re.printStackTrace();
            int line = (line = re.lineNumber()) == 0 ? -1 : line;
            String msg;
            if (re instanceof JavaScriptException) {
                msg = String.valueOf(((JavaScriptException) re).getValue());
            } else {
                msg = re.toString();
            }
            ScriptException se = new ScriptException(msg, re.sourceName(), line);
            se.initCause(re);
            throw se;
        } catch (IOException ee) {
            throw new ScriptException(ee);
        } finally {
            Context.exit();
        }

        return unwrapReturnValue(ret);
    }

    Scriptable getRuntimeScope(ScriptContext ctxt) {
        if (ctxt == null) {
            throw new NullPointerException("null script context");
        }

        // we create a scope for the given ScriptContext
        Scriptable newScope = new ExternalScriptable(ctxt, indexedProps);

        // Set the prototype of newScope to be 'topLevel' so that
        // JavaScript standard objects are visible from the scope.
        newScope.setPrototype(stdObjects);
        newScope.setParentScope(null);

        return newScope;
    }

    Object unwrapReturnValue(Object result) {
        if (result instanceof Wrapper) {
            result = ((Wrapper) result).unwrap();
        }

        return result instanceof Undefined ? null : result;
    }

    @Override
    public Object invokeFunction(String name, Object... args) throws ScriptException, NoSuchMethodException {
        return invoke(null, name, args);
    }

    @Override
    public Object invokeMethod(Object thiz, String name, Object... args) throws ScriptException, NoSuchMethodException {
        if (thiz == null) {
            throw new IllegalArgumentException("script object can not be null");
        }
        return invoke(thiz, name, args);
    }

    private Object invoke(Object thiz, String name, Object... args)
            throws ScriptException, NoSuchMethodException {
        Context cx = Context.enter();
        try {
            if (name == null) {
                throw new NullPointerException("method name is null");
            }

            if (thiz != null && !(thiz instanceof Scriptable)) {
                thiz = cx.toObject(thiz, stdObjects);
            }

            Scriptable engineScope = getRuntimeScope(context);
            Scriptable localScope = (thiz != null) ? (Scriptable) thiz
                    : engineScope;
            Object obj = ScriptableObject.getProperty(localScope, name);
            if (!(obj instanceof Function)) {
                throw new NoSuchMethodException("no such method: " + name);
            }

            Function func = (Function) obj;
            Scriptable scope = func.getParentScope();
            if (scope == null) {
                scope = engineScope;
            }
            Object result = func.call(cx, scope, localScope,
                    wrapArguments(args));
            return unwrapReturnValue(result);
        } catch (RhinoException re) {
            re.printStackTrace();
            int line = (line = re.lineNumber()) == 0 ? -1 : line;
            throw new ScriptException(re.toString(), re.sourceName(), line);
        } finally {
            cx.exit();
        }
    }

    Object[] wrapArguments(Object[] args) {
        if (args == null) {
            return Context.emptyArgs;
        }
        Object[] res = new Object[args.length];
        for (int i = 0; i < res.length; i++) {
            res[i] = Context.javaToJS(args[i], stdObjects);
        }
        return res;
    }

    @Override
    public Bindings createBindings() {
        return new SimpleBindings();
    }

    @Override
    public ScriptEngineFactory getFactory() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public Object eval(String script, ScriptContext context) throws ScriptException {
        return eval(new StringReader(script), context);
    }
}
