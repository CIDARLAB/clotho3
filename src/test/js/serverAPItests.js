
module("Server API Tests");

var Message = function (channel, data, requestId) {
    this.channel = channel;
    this.data = data;
    if (channel == "submit" && typeof this.data == "string"){
        this.data = {query:this.data, tokens:[]};
    } 
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
        new Message("get", "org.clothocad.model.Part"),
        function (data) {
            equal(data.name, "Part");
        });

testThroughAsync("query Parts",
        new Message("query", {"schema":"org.clothocad.model.Part"}),
        function (data) {
            equal(data.length, 55);
        });

testThroughAsync("query BasicParts",
        new Message("query", {"schema":"org.clothocad.model.BasicPart"}),
        function (data) {
            equal(data.length, 54);
        });

testThroughAsync("query CompositeParts",
        new Message("query", {"schema":"org.clothocad.model.CompositePart"}),
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
        socket.send( new Message("create", {"name":"Created Part 2", "sequence":"CCCC", "schema":"org.clothocad.model.BasicPart"}), function (data) {
            socket.send(new Message("get", data), function (data2) {
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
        new Message("create", {"schema":"LabPerson"}),
        function (){
            expect(0);
        });

testThroughAsync("validation failure",
        new Message("validate", {schema:"org.clothocad.model.NucSeq", sequence:null}),
        function (errors){
            equal(errors.length, 1);
            equal(errors[0].message, "may not be null");
            equal(errors[0].propertyPath.currentLeafNode.name, "sequence");
        });

testThroughAsync("validation success", 
        new Message("validate", {schema:"org.clothocad.model.NucSeq", sequence:"ATCG"}),
        function (errors){
            equal(errors.length, 0);
        });

asyncTest("authoring schema constraints", function (){
    var socket = getSocket(clothosocket);
    socket.onopen = function (){
        socket.send(new Message("create", {id:"org.clothocad.model.SimpleFeature", language:"JSONSCHEMA", schema:"org.clothocad.core.schema.ClothoSchema", name:"SimpleFeature", description:"A simple and sloppy representation of a Feature or other DNA sequence", "fields":[{name:"sequence", type:"string", example:"ATACCGGA", access:"PUBLIC", constraints:[{constraintType:"javax.validation.constraints.Pattern", values:{flags:["CASE_INSENSITIVE"], regexp:"[ATUCGRYKMSWBDHVN]*"}}], description:"the sequence of the feature"}]}), function (response) {
            socket.send(new Message("validate", {schema:"org.clothocad.model.SimpleFeature", sequence:"this is not a valid sequence"}), function (errors) {
                equal(errors.length, 1);
                equal(errors[0].message, "must match \"[ATUCGRYKMSWBDHVN]*\"");
                equal(errors[0].propertyPath.currentLeafNode.name, "sequence");
                socket.send(new Message("destroy", "org.clothocad.model.SimpleFeature"), function(response){
                    start();
                });
            });
        });
    };
});

asyncTest("set", function (){
    var socket = getSocket(clothosocket);
    socket.onopen = function () {
        socket.send(new Message("get", "org.clothocad.model.Part"), function (data) {
            var id = data.id;
            socket.send(new Message("set", {"id":id, "name":"Set Part"}), function(data){
                socket.send(new Message("get", id), function (data){
                    equal(data.name, "Set Part");
                    socket.send(new Message("set", {"id":id, "name":"Part"}), function(data){
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

//this is probably not necessary now that schemas have deterministic ids
asyncTest("reload models", function(){
    var socket = getSocket(clothosocket);
    socket.onopen = function () {
        socket.send(new Message("get", "org.clothocad.core.schema.ClothoSchema"), function(data){
            id = data.id;
            socket.send(new Message("reloadModels", {}), function (response) {
                socket.send(new Message("get", "org.clothocad.core.schema.ClothoSchema"), function(data){
                    equal(data.id, id);
                    start();
                });
            });
        });
    };
});

module("Functions and Modules")

asyncTest("changing function", function(){
    var f1 = {schema:"org.clothocad.core.datums.Function",language:"JAVASCRIPT", code:"function(){return 1;};"};
    var f2 = {schema:"org.clothocad.core.datums.Function",language:"JAVASCRIPT", code:"function(){return 2;};"};

    var socket = getSocket(clothosocket);
    socket.onopen = function () {
        socket.send(new Message("create", f1), function(id){
            socket.send(new Message("run", {"id":id, args:[]}), function(result){
                equal(result, 1);
                f2.id = id;
                socket.send(new Message("set", f2), function(data){
                    socket.send(new Message("run", {"id":id, args:[]}), function(result){
                        equal(result,2);
                        start();
                    });
                });
            });
        });
    };
});

testThroughAsync("functions with arguments",
        new Message("submit", "clotho.run('org.clothocad.test.lowercase', ['HEY'])"),
        function(data){
            equal(data, "hey");
        });
testThroughAsync("module scoping",
        new Message("run", {id:"org.clothocad.test.moduleTestFunction", args:[1]}),
        function(data){
            equal(data, 4);
        });

testThroughAsync("invoke module method",
        new Message("run", {id:"org.clothocad.test.testModule", "function":"moduleMethod", args:[]}),
        function(data){
            equal(data, 2);
        });

testThroughAsync("module depends on module",
        new Message("run", {id:"org.clothocad.test.testModule2","function":"moduleMethod", args:[]}),
        function(data){
            equal(data,2);
        });

testThroughAsync("function depends on function",
        new Message("run", {id:"org.clothocad.test.moduleTestFunction", args:[1]}),
        function(data){
            equal(data,4);
        });

//TODO: test case for multiple dependencies

testThroughAsync("use lodash",
        new Message("run", {id: "org.clothocad.test.useLodash", args:[]}),
        function(data){
            deepEqual(data, [2,3,4]);
        });

asyncTest("package module", function (){
    var socket = getSocket(clothosocket);
    socket.onopen = function () {
        socket.send(new Message("run", {id:"packager", "function":"packJavaScript",
            //arg order is name, description, code, dependencies (as array)
            args:["packMe", "a module", "(function f() {var my = {}; var privateVariable = 2; function privateMethod() { return privateVariable; } my.moduleProperty = 1; my.moduleMethod = function () { return privateMethod(); }; return my; }());", []]}), function (results){
                //returns the created module by means of an id - use a get to grab the new module
                var id = results[0];
                socket.send(new Message("get", id), function (result){
            //this doesn't test anything interesting, but you can put a breakpoint here to inspect the returned packaged module
                    equal(result.name, "packMe");
                    start();
                });
            });
    };
});


module("SynBERC demo issues - all through submit channel")

testThroughAsync("console.log",
        new Message("submit", "clotho.run('org.clothocad.test.consoleTest',[])"),
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
        new Message("submit", "clotho.run('org.clothocad.test.clothoLoadTest', [])"),
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


