import importlib
import json
import sys

class ClothoError(Exception):
    pass

class ClothoMethod(object):
    def __init__(self, name, context):
        self._name = name
        self._context = context

    def __call__(self, *args):
        self._context.send_value(
            {"type": "api", "name": self._name, "args": args})
        reply_value = read_value()
        reply_type = reply_value["type"]
        if reply_type == "api":
            return reply_value["return"]
        elif reply_type == "api_error":
            raise ClothoError(reply_value["message"])
        raise RuntimeError("bad API reply")

class Clotho(object):
    def __init__(self, context):
        self._context = context

    def __getattr__(self, name):
        return ClothoMethod(name, self._context)

class Context(object):
    '''Execution context object

    Setup by initialization object
    '''
    def __init__(self, init_obj):
        if init_obj["type"] != "func":
            raise ValueError
        self.code = init_obj["code"]
        self.args = init_obj["args"]
        self._tmpfile = open(init_obj["tmpfile"], "wb")

    def __del__(self):
        self._tmpfile.flush();
        self._tmpfile.close();

    def send_value(self, value):
        '''Send one message to host through temporary file'''
        self._tmpfile.write(json.dumps(value).encode("UTF-8"))
        self._tmpfile.write(b'\0')
        self._tmpfile.flush()

def read_value():
    '''Read one message from host through standard input and return it'''
    buf = bytearray()
    while 1:
        c = sys.stdin.read(1)
        if c == "":
            raise ValueError("received incomplete JSON value")
        if c == "\0":
            break
        buf.append(c)
    return json.loads(bytes(buf).decode("UTF-8"))

def main():
    # fetch initialization object
    context = Context(read_value())

    # execute user code body
    scope_dict = {"clotho": Clotho(context), "ClothoError": ClothoError}
    exec(context.code, scope_dict)

    # execute user function, which is always called "run"
    user_func = scope_dict["run"]
    user_return_value = user_func(*context.args)

    context.send_value({"type": "func", "return": user_return_value})

main()
