import importlib
import json
import sys

class IncompleteReadError(Exception):
    pass

class Clotho(object):
    @staticmethod
    def get(id):
        send_val({"type": "api", "name": "get", "args": [id]})
        api_return = read_val()
        assert api_return["type"] == "api"
        return api_return["return"]

    @staticmethod
    def set(val):
        send_val({"type": "api", "name": "set", "args": [val]})
        api_return = read_val()
        assert api_return["type"] == "api"
        return api_return["return"]

    @staticmethod
    def run(name, args):
        send_val({"type": "api", "name": "run", "args": [name, args]})
        api_return = read_val()
        assert api_return["type"] == "api"
        return api_return["return"]

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
    scope_dict = {"clotho": Clotho()}
    exec(defcall_obj["code"], scope_dict)

    # execute user function, which is always called "run"
    user_func = scope_dict["run"]
    user_return_value = user_func(*defcall_obj["args"])

    send_val({"type": "func", "return": user_return_value})

main()
