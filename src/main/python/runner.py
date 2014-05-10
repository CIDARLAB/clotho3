import importlib
import json
import sys

class IncompleteReadError(Exception):
    pass

class ClothoError(Exception):
    pass

class ClothoMethod(object):
    def __init__(self, name):
        self._name = name

    def __call__(self, *args):
        send_val({"type": "api", "name": self._name, "args": args})
        reply_value = read_val()
        reply_type = reply_value["type"]
        if reply_type == "api":
            return reply_value["return"]
        elif reply_type == "api_error":
            raise ClothoError(reply_value["message"])
        raise RuntimeError("bad API reply")

class Clotho(object):
    def __getattr__(self, name):
        return ClothoMethod(name)

def read_val():
    buf = bytearray()
    while 1:
        c = sys.stdin.read(1)
        if c == "":
            raise IncompleteReadError
        if c == "\0":
            break
        buf.append(c)
    return json.loads(bytes(buf).decode("UTF-8"))

def send_val(val):
    sys.stdout.write(json.dumps(val).encode("UTF-8"))
    sys.stdout.write(b'\0')
    sys.stdout.flush()

def main():
    # fetch defcall object
    defcall_obj = read_val()
    assert defcall_obj["type"] == "func"

    # execute user code body
    scope_dict = {"clotho": Clotho(), "ClothoError": ClothoError}
    exec(defcall_obj["code"], scope_dict)

    # execute user function, which is always called "run"
    user_func = scope_dict["run"]
    user_return_value = user_func(*defcall_obj["args"])

    send_val({"type": "func", "return": user_return_value})

main()
