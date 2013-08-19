
module("Server API Tests");

var Message = function (channel, data, requestId) {
    this.channel = channel;
    this.data = data;
};

var getSocket = function (addr) {
    socket = new WebSocket(addr);
    socket.idx = -1;
    socket.callbacks = {};
        
    socket.oldsend = socket.send;

    socket.send = function (message, callback) {
        socket.idx += 1;
        message.requestId = String(socket.idx);
        socket.callbacks[message.channel+message.requestId] = callback;
        socket.oldsend(JSON.stringify(message));
    };
    socket.onmessage = function (e) {
        var message = JSON.parse(e.data);
        var callbackKey = message.channel+message.requestId;
        if (socket.callbacks.hasOwnProperty(callbackKey)){
            callback = socket.callbacks[callbackKey];
            delete socket.callbacks[callbackKey];
            callback(message.data)
        }
    };
    return socket;
}

var clothosocket = "ws://localhost:8080/websocket";

var testThroughAsync = function (name, message, callback) {
    asyncTest(name, function () {
        var socket = getSocket(clothosocket);
        socket.onopen = function () {
            socket.send(message, function (data) {
                callback(data);
                start();
            });
        };
    });
};

testThroughAsync("get",
        new Message("get", "Test Part 1"),
        function (data) {
            equal(data.name, "Test Part 1");
        });

testThroughAsync("query Parts",
        new Message("query", {"schema":"Part"}),
        function (data) {
            equal(data.length, 4);
        });

testThroughAsync("query BasicParts",
        new Message("query", {"schema":"BasicPart"}),
        function (data) {
            equal(data.length, 3);
        });

testThroughAsync("query CompositeParts",
        new Message("query", {"schema":"CompositePart"}),
        function (data) {
            equal(data.length, 1); 
            data = data[0]
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
        socket.send( new Message("create", {"name":"Created Part", "sequence":"GGGGGG"}), function (data) {
            var id = data;
            socket.send(new Message("get", id), function (data) {
                equal(data.sequence, "GGGGGG");
                start();
            });
        })
    };
});

asyncTest("create with schema", function () {
    var socket = getSocket(clothosocket);
    socket.onopen = function () {
        socket.send( new Message("create", {"name":"Created Part 2", "sequence":"CCCC", "schema":"BasicPart"}), function (data) {
            socket.send(new Message("get", "Created Part 2"), function (data2) {
                ok(data2.hasOwnProperty("schema"));
                ok(!(data2.hasOwnProperty("className")));
                //tear down created data
                socket.send(new Message("destroy", "Created Part 2"), function(data){
                    //TODO: ensure actually deleted
                    start();
                });

            });
        })
    };
});

asyncTest("set", function (){
    var socket = getSocket(clothosocket);
    socket.onopen = function () {
        socket.send(new Message("get", "Test Part 1"), function (data) {
            var id = data.id;
            socket.send(new Message("set", {"id":id, "name":"Set Part"}), function(data){
                socket.send(new Message("get", id), function (data){
                    equal(data.name, "Set Part");
                    socket.send(new Message("set", {"id":id, "name":"Test Part 1"}), function(data){
                        start();
                    });
                });
            });
        });
    };
});

// TODO: add tests for listener dereg
