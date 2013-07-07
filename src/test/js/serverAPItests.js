
module("Server API Tests");

var Message = function (channel, data, requestId) {
    this.channel = channel;
    this.data = data;
    this.requestId = requestId;
};

var send = function (message, socket, callback) {
    socket.router.add(callback, message.channel, message.requestId);
    socket.send(JSON.stringify(message));
};

var Router = function () {
    this.add = function (callback, channel, id) {
        this[channel + id] = callback;
    };

    this.onmessage = function (e) {
        var message = JSON.parse(e.data);
        var callbackKey = message.channel + message.requestId;
        if (this.hasOwnProperty(callbackKey)){
            this[callbackKey](message.data);
        }
    };
};

var getSocket = function (addr) {
    var router = new Router();
    var socket = new WebSocket(addr);
    socket.router = router;
    socket.onmessage = function (e) {router.onmessage(e)};
    //socket.send = function (message){
    //    socket.send(JSON.stringify(message));
    //}
    return socket;
};

var clothosocket = "ws://localhost:8080/websocket";

var testThroughAsync = function (name, message, callback) {
    asyncTest(name, function () {
        var socket = getSocket(clothosocket);
        socket.onopen = function () {
            send(message, socket, function (data) {
                callback(data);
                start();
            });
        };
    });
};

/*
   asyncTest("get", function () {
   var m = new Message("get", "Test Part 1", "1");
   var socket = getSocket(clothosocket);
   socket.onopen = function () {
   send(m, socket, function (data) {
   equal(data.name, "Test Part 1");
   start();
   });
   };
   });
   */

testThroughAsync("get",
        new Message("get", "Test Part 1", "2"),
        function (data) {
            equal(data.name, "Test Part 1");
        });

testThroughAsync("query Parts",
        new Message("query", {"schema":"Part"}, "3"),
        function (data) {
            equal(data.length, 4);
        });

testThroughAsync("query BasicParts",
        new Message("query", {"schema":"BasicPart"}, "4"),
        function (data) {
            equal(data.length, 3);
        });

testThroughAsync("query CompositeParts",
        new Message("query", {"schema":"CompositePart"}, "5"),
        function (data) {
            //equal(data.length, 1); -- because we unwrap the result
            equal(data.type, "COMPOSITE");
            //shitty hack to test if string
            ok(data.composition[0].substring);
        });

/*
 * TODO: some kind of promise test
 testThroughAsync("create",
 new Message("create", {"name":"Created Part", "sequence":"GGGGGG"}, "6"),
 function (data) {
 var id = data;
 var message = new Message("get", id, "7");
 send(message, this.socket, function (data) {
 equal(data.sequence, "GGGGGG");
 start();
 });
 });
 */

asyncTest("create", function () {
    var socket = getSocket(clothosocket);
    socket.onopen = function () {
        send( new Message("create", {"name":"Created Part", "sequence":"GGGGGG"}, "6"), socket, function (data) {
            var id = data;
            send(new Message("get", id, "7"), socket, function (data) {
                equal(data.sequence, "GGGGGG");
                start();
            });
        })
    };
});
