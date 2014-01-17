
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

var clothosocket = "wss://localhost:8443/websocket";

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
            equal(data.length, 55);
        });

testThroughAsync("query BasicParts",
        new Message("query", {"schema":"BasicPart"}),
        function (data) {
            equal(data.length, 54);
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

testThroughAsync("create LabPerson", 
        new Message("create", {"schema":"52310822c2e67bf02c9c21f1"}),
        function (){
            expect(0);
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

asyncTest("login/logout", function(){
    var credentials = {"username":"testuser", "password":"password"};

    var socket = getSocket(clothosocket);
    socket.onopen = function () {
        socket.send(new Message("login" , credentials), function (data){
            equal(data, true);
            socket.send(new Message("submit", "var persistMe = 42"), function (data){
                socket.send(new Message("submit", "persistMe"), function (data){
                    equal(data, 42);
                    socket.send(new Message("logout", ""), function (data) {
                        equal(data,true);
                        socket.send(new Message("submit", "persistMe"), function (data) {
                            notEqual(data, 42, "data retained on logout?");
                            socket.send(new Message("login", credentials), function(data){
                                equal(data,true);
                                socket.send(new Message("submit", "persistMe"), function (data){
                                    equal(data, 42, "data recovered on login?");
                                    start();
                                });
                            });
                        });
                    });
                });
            });
        });
    };
});
// TODO: add tests for listener dereg

asyncTest("reload models", function(){
    var socket = getSocket(clothosocket);
    socket.onopen = function () {
        socket.send(new Message("get", "ClothoSchema"), function(data){
            id = data.id;
            socket.send(new Message("reloadModels", {}), function (response) {
                socket.send(new Message("get", "ClothoSchema"), function(data){
                    equal(data.id, id);
                    start();
                });
            });
        });
    };
});

module("Functions and Modules")

asyncTest("changing function", function(){
    var f1 = {"schema":"Function","language":"JAVASCRIPT", "code":"function(){return 1;};"};
    var f2 = {schema:"Function",language:"JAVASCRIPT", "code":"function(){return 2;};"};

    var socket = getSocket(clothosocket);
    socket.onopen = function () {
        socket.send(new Message("create", f1), function(id){
            socket.send(new Message("run", {"id":id, "args":[]}), function(result){
                equal(result, 1);
                f2.id = id;
                socket.send(new Message("set", f2), function(data){
                    socket.send(new Message("run", {id:id, args:[]}), function(result){
                        equal(result,2);
                        start();
                    });
                });
            });
        });
    };
});

testThroughAsync("functions with arguments",
        new Message("submit", "clotho.run('lowercase', ['HEY'])"),
        function(data){
            equal(data, "hey");
        });
testThroughAsync("module scoping",
        new Message("run", {id:"moduleTestFunction", args:[1]}),
        function(data){
            equal(data, 4);
        });

testThroughAsync("invoke module method",
        new Message("run", {id:"testModule", "function":"moduleMethod", args:[]}),
        function(data){
            equal(data, 2);
        });

testThroughAsync("module depends on module",
        new Message("run", {id:"testModule2","function":"moduleMethod", args:[]}),
        function(data){
            equal(data,2);
        });

testThroughAsync("function depends on function",
        new Message("run", {id:"moduleTestFunction", args:[1]}),
        function(data){
            equal(data,4);
        });

testThroughAsync("importing external library",
        new Message("submit", "load(\"yepnope.js\")"),
        function (data) {
            equal(data, null);
        });
testThroughAsync("use lodash",
        new Message("run", {id: "useLodash", args:[]}),
        function(data){
            deepEqual(data, [2,3,4]);
        });

module("SynBERC demo issues - all through submit channel")

testThroughAsync("console.log",
        new Message("submit", "clotho.run('consoleTest',[])"),
        function(data){
            equal(data, 'This worked!');
        });

testThroughAsync("running module functions - revcomp",
        new Message("submit", "clotho.run('DNA', 'revcomp', ['acgtac'])"),
        function(data){
            equal(data, 'gtacgt');
        });

testThroughAsync("running module functions - ligate",
        new Message("submit", "clotho.run('PCR', 'ligate', [['aaaaaaaaaaA^CATG_', '^CATG_Tttggttggttgg']])"),
        function(data){
            equal(data, "aaaaaaaaaaACATGTttggttggttgg");
        });

testThroughAsync("loading functions - submit", 
        new Message("submit", "var DNA = clotho.load('DNA'); DNA.revcomp('acgtacg')"),
        function (data){
            equal(data, "cgtacgt");
        });

testThroughAsync("loading functions - load in function", 
        new Message("submit", "clotho.run('clothoLoadTest', [])"),
        function (data){
            equal(data, 'cgtacgt');
        });

testThroughAsync("loading functions - global load through function",
        new Message("submit", "clotho.run('clientSetup', []); DNA.revcomp('ttttacccggg');"),
        function (data) {
            equal(data, 'cccgggtaaaa');
        });

testThroughAsync("functions on objects",
        new Message("submit", "var myObj = {}; myObj.myFunc = function(str) { return 'hey ' + str; }; myObj.myFunc('Bob')"),
        function (data) {
            equal(data, 'hey Bob');
        });
