/* Low-level client<-->server communication layer
 *
 * This layer should not be accessed directly. Instead,
 * use libsend and librecv.
 */

var libtransport = new Object();
libtransport.transports = new Object();

libtransport.transports.WebSocket = function () {
    this.send = function (data) {
        switch (websocket.readyState) {
        case WebSocket.OPEN:
            websocket.send(data);
            break;
        default:
            alert("The socket is not open.");
        }
    };

    function init () {
        websocket = new WebSocket(WEBSOCKET_PATH);

        /* TODO: apparently this function can be racy
             if it fires before libtransport.emitSignal is overridden */
        websocket.onmessage = function (event) {
            var obj = JSON.parse(event.data);
            libtransport.handleMessage(obj);
        };
        websocket.onopen = function (event) {
            libtransport.send(libtransport.INTERNAL_CHANNEL, "", "");
        };
        websocket.onclose = function (event) {
            libtransport.setSocketID("");
        };
    };

    var WEBSOCKET_PATH = "ws://localhost:1025/servlet/websocket";
    var websocket = null;
    init();
}

libtransport.transports.XHR = function () {
    this.send = function (data) {
        var xhr_obj = new XMLHttpRequest();
        xhr_obj.onreadystatechange = function () {
            switch (this.readyState) {
            case XMLHttpRequest.DONE:
                var obj = JSON.parse(this.responseText);
                for (i in obj) {
                    libtransport.handleMessage(obj[i]);
                }
            }
        }
        xhr_obj.open(XHR_METHOD, XHR_PATH);
        xhr_obj.send(data);
    }

    function init() {
        clearInterval(pollhandle);
        var evalstr = "libtransport.send('" +
                      libtransport.INTERNAL_CHANNEL +
                      "', '', '');" ;
        pollhandle = setInterval(evalstr, POLL_PERIOD_MS);
    }

    var XHR_PATH = "/servlet/xhr";
    var XHR_METHOD = "POST";
    var POLL_PERIOD_MS = 1000;
    var pollhandle = null;
    init();
}

/* Low-level client-->server message send.
 * Called by libsend
 */
libtransport.send = function (channel, message, flags) {
    /* TODO: ensure that arguments are strings */
    var obj = {"channel": channel,
               "message": message,
               "flags" : flags,
               "socket_id" : libtransport.getSocketID()};
    
    if (obj.channel != libtransport.INTERNAL_CHANNEL && obj.socket_id == "") {
        alert("Lost socket ID");
        return;
    }
    data = JSON.stringify(obj);
    libtransport._transport.send(data);
};

/* Called by specific transport (above) on server-->client messages. */
/* TODO: window.name might not be empty string by default */
libtransport.handleMessage = function (obj) {
    if (obj.channel == libtransport.INTERNAL_CHANNEL) {
        libtransport.setSocketID(obj.socket_id);
        var out = {"ephemeral_link_page_id": window.name,
                   "page_mode": libpage.getCurrentMode()};
        libsend.call("linkPage", JSON.stringify(out));
    } else {
        libtransport.receive(obj.channel, obj.message, obj.flags);
    }
};

/* librecv provides concrete implementation.
 * Called by handleMessage() (above) on server-->client messages.
 */
libtransport.receive = function (channel, message, flags) {
    /* TODO: decide whether server-to-client flags are necessary */
    alert("libtransport.receive: Missing external implementation. why???");
};

/* libsignal provides concrete implementation.
 * Called when any purely-client-side signal needs to be emitted.
 */
libtransport.emitSignal = function (notification) {
    alert("libtransport.emitSignal: Missing external implementation");
}

libtransport.getSocketID = function () {
    return libtransport.socket_id;
}

libtransport.setSocketID = function (new_id) {
    libtransport.socket_id = new_id;
    libtransport.emitSignal("socketIDChange", "");
}

libtransport.init = function () {
    libtransport.INTERNAL_CHANNEL = "_transport";
    if (window.WebSocket) {
        libtransport._transport = new libtransport.transports.WebSocket();
    } else {
        libtransport._transport = new libtransport.transports.XHR();
    }

    libtransport.socket_id = "";
};
