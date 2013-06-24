#!/usr/bin/env node

// added socket.io stuff to the bottom

var util = require('util'),
    http = require('http'),
    io = require('socket.io').listen(8090),
    fs = require('fs'),
    url = require('url'),
    events = require('events'),
    path = require('path')
    ;

var DEFAULT_PORT = 8000;


function main(argv) {
    new HttpServer({
        'GET': createServlet(StaticServlet),
        'HEAD': createServlet(StaticServlet)
    }).start(Number(argv[2]) || DEFAULT_PORT);
}

function escapeHtml(value) {
    return value.toString().
        replace('<', '&lt;').
        replace('>', '&gt;').
        replace('"', '&quot;');
}

function createServlet(Class) {
    var servlet = new Class();
    return servlet.handleRequest.bind(servlet);
}

/**
 * An Http server implementation that uses a map of methods to decide
 * action routing.
 *
 * @param {Object} Map of method => Handler function
 */
function HttpServer(handlers) {
    this.handlers = handlers;
    this.server = http.createServer(this.handleRequest_.bind(this));
}

HttpServer.prototype.start = function(port) {
    this.port = port;
    this.server.listen(port);
    util.puts('Http Server running at http://localhost:' + port + '/');
};

HttpServer.prototype.parseUrl_ = function(urlString) {
    var parsed = url.parse(urlString);
    parsed.pathname = url.resolve('/', parsed.pathname);
    return url.parse(url.format(parsed), true);
};

HttpServer.prototype.handleRequest_ = function(req, res) {
    var logEntry = req.method + ' ' + req.url;
    if (req.headers['user-agent']) {
        logEntry += ' ' + req.headers['user-agent'];
    }
    util.puts(logEntry);
    req.url = this.parseUrl_(req.url);
    var handler = this.handlers[req.method];
    if (!handler) {
        res.writeHead(501);
        res.end();
    } else {
        handler.call(this, req, res);
    }
};

/**
 * Handles static content.
 */
function StaticServlet() {}

StaticServlet.MimeMap = {
    'txt': 'text/plain',
    'html': 'text/html',
    'css': 'text/css',
    'xml': 'application/xml',
    'json': 'application/json',
    'js': 'application/javascript',
    'jpg': 'image/jpeg',
    'jpeg': 'image/jpeg',
    'gif': 'image/gif',
    'png': 'image/png',
    'svg': 'image/svg+xml'
};

StaticServlet.prototype.handleRequest = function(req, res) {
    var self = this;
    var path = ('./' + req.url.pathname).replace('//','/').replace(/%(..)/g, function(match, hex){
        return String.fromCharCode(parseInt(hex, 16));
    });
    var parts = path.split('/');
    if (parts[parts.length-1].charAt(0) === '.')
        return self.sendForbidden_(req, res, path);
    fs.stat(path, function(err, stat) {
        if (err)
            return self.sendMissing_(req, res, path);
        if (stat.isDirectory())
            return self.sendDirectory_(req, res, path);
        return self.sendFile_(req, res, path);
    });
};

StaticServlet.prototype.sendError_ = function(req, res, error) {
    res.writeHead(500, {
        'Content-Type': 'text/html'
    });
    res.write('<!doctype html>\n');
    res.write('<title>Internal Server Error</title>\n');
    res.write('<h1>Internal Server Error</h1>');
    res.write('<pre>' + escapeHtml(util.inspect(error)) + '</pre>');
    util.puts('500 Internal Server Error');
    util.puts(util.inspect(error));
};

StaticServlet.prototype.sendMissing_ = function(req, res, path) {
    path = path.substring(1);
    res.writeHead(404, {
        'Content-Type': 'text/html'
    });
    res.write('<!doctype html>\n');
    res.write('<title>404 Not Found</title>\n');
    res.write('<h1>Not Found</h1>');
    res.write(
        '<p>The requested URL ' +
            escapeHtml(path) +
            ' was not found on this server.</p>'
    );
    res.end();
    util.puts('404 Not Found: ' + path);
};

StaticServlet.prototype.sendForbidden_ = function(req, res, path) {
    path = path.substring(1);
    res.writeHead(403, {
        'Content-Type': 'text/html'
    });
    res.write('<!doctype html>\n');
    res.write('<title>403 Forbidden</title>\n');
    res.write('<h1>Forbidden</h1>');
    res.write(
        '<p>You do not have permission to access ' +
            escapeHtml(path) + ' on this server.</p>'
    );
    res.end();
    util.puts('403 Forbidden: ' + path);
};

StaticServlet.prototype.sendRedirect_ = function(req, res, redirectUrl) {
    res.writeHead(301, {
        'Content-Type': 'text/html',
        'Location': redirectUrl
    });
    res.write('<!doctype html>\n');
    res.write('<title>301 Moved Permanently</title>\n');
    res.write('<h1>Moved Permanently</h1>');
    res.write(
        '<p>The document has moved <a href="' +
            redirectUrl +
            '">here</a>.</p>'
    );
    res.end();
    util.puts('301 Moved Permanently: ' + redirectUrl);
};

StaticServlet.prototype.sendFile_ = function(req, res, path) {
    var self = this;
    var file = fs.createReadStream(path);
    res.writeHead(200, {
        'Content-Type': StaticServlet.
            MimeMap[path.split('.').pop()] || 'text/plain'
    });
    if (req.method === 'HEAD') {
        res.end();
    } else {
        file.on('data', res.write.bind(res));
        file.on('close', function() {
            res.end();
        });
        file.on('error', function(error) {
            self.sendError_(req, res, error);
        });
    }
};

StaticServlet.prototype.sendDirectory_ = function(req, res, path) {
    var self = this;
    if (path.match(/[^\/]$/)) {
        req.url.pathname += '/';
        var redirectUrl = url.format(url.parse(url.format(req.url)));
        return self.sendRedirect_(req, res, redirectUrl);
    }
    fs.readdir(path, function(err, files) {
        if (err)
            return self.sendError_(req, res, error);

        if (!files.length)
            return self.writeDirectoryIndex_(req, res, path, []);

        var remaining = files.length;
        files.forEach(function(fileName, index) {
            fs.stat(path + '/' + fileName, function(err, stat) {
                if (err)
                    return self.sendError_(req, res, err);
                if (stat.isDirectory()) {
                    files[index] = fileName + '/';
                }
                if (!(--remaining))
                    return self.writeDirectoryIndex_(req, res, path, files);
            });
        });
    });
};

StaticServlet.prototype.writeDirectoryIndex_ = function(req, res, path, files) {
    path = path.substring(1);
    res.writeHead(200, {
        'Content-Type': 'text/html'
    });
    if (req.method === 'HEAD') {
        res.end();
        return;
    }
    res.write('<!doctype html>\n');
    res.write('<title>' + escapeHtml(path) + '</title>\n');
    res.write('<style>\n');
    res.write('  ol { list-style-type: none; font-size: 1.2em; }\n');
    res.write('</style>\n');
    res.write('<h1>Directory: ' + escapeHtml(path) + '</h1>');
    res.write('<ol>');
    files.forEach(function(fileName) {
        if (fileName.charAt(0) !== '.') {
            res.write('<li><a href="' +
                escapeHtml(fileName) + '">' +
                escapeHtml(fileName) + '</a></li>');
        }
    });
    res.write('</ol>');
    res.end();
};

/** UTILITIES **/

function isUndefined(value){return typeof value == 'undefined';}

function isDefined(value){return typeof value != 'undefined';}

function isObject(value){return value != null && typeof value == 'object';}

function isString(value){return typeof value == 'string';}

function isNumber(value){return typeof value == 'number';}

function isDate(value){
    return toString.apply(value) == '[object Date]';
}

function isArray(value) {
    return toString.apply(value) == '[object Array]';
}

function isFunction(value){return typeof value == 'function';}

function isWindow(obj) {
    return obj && obj.document && obj.location && obj.alert && obj.setInterval;
}


function isScope(obj) {
    return obj && obj.$evalAsync && obj.$watch;
}


function isFile(obj) {
    return toString.apply(obj) === '[object File]';
}


function isBoolean(value) {
    return typeof value == 'boolean';
}
function equals(o1, o2) {
    if (o1 === o2) return true;
    if (o1 === null || o2 === null) return false;
    if (o1 !== o1 && o2 !== o2) return true; // NaN === NaN
    var t1 = typeof o1, t2 = typeof o2, length, key, keySet;
    if (t1 == t2) {
        if (t1 == 'object') {
            if (isArray(o1)) {
                if ((length = o1.length) == o2.length) {
                    for(key=0; key<length; key++) {
                        if (!equals(o1[key], o2[key])) return false;
                    }
                    return true;
                }
            } else if (isDate(o1)) {
                return isDate(o2) && o1.getTime() == o2.getTime();
            } else {
                if (isScope(o1) || isScope(o2) || isWindow(o1) || isWindow(o2)) return false;
                keySet = {};
                for(key in o1) {
                    if (key.charAt(0) === '$' || isFunction(o1[key])) continue;
                    if (!equals(o1[key], o2[key])) return false;
                    keySet[key] = true;
                }
                for(key in o2) {
                    if (!keySet[key] &&
                        key.charAt(0) !== '$' &&
                        o2[key] !== undefined &&
                        !isFunction(o2[key])) return false;
                }
                return true;
            }
        }
    }
    return false;
}


/** socket.io
*
* NOTE - CLIENT CODE REFERENCE
 * socket.io.js on client code is served dynamically by node (it's not a resource). However, this web-server does url rewriting, so you can't use it as normal:
 * <script src="/socket.io/socket.io.js"></script>
 * but instead you need to serve it like so:
 * <script src="http://nodeJS_server:port/socket.io/socket.io.js"></script>
 * where nodeJS_server is probably localhost
 * and port is the port used for socket.io (8090), not node itself (8000)
 * -- numbers given may change, but are defined (currently) in this file
*
 * NOTE - sending over data
 * emit() can send a string or JSON
 * send() can only send a string
 * WebSockets only sends strings
* */

// simple user / connection management for socket.io implementation
var userNames = (function () {
    var names = {};

    //check availability
    var claim = function (name) {
        if (!name || names[name]) {
            return false;
        } else {
            names[name] = true;
            return true;
        }
    };

    // find the lowest unused "guest" name and claim it
    var getGuestName = function () {
        var name,
            nextUserId = 1;

        do {
            name = 'Guest ' + nextUserId;
            nextUserId += 1;
        } while (!claim(name));

        return name;
    };

    // serialize claimed names as an array
    var get = function () {
        var res = [];
        for (var user in names) {
            res.push(user);
        }

        return res;
    };

    //remove from list
    var free = function (name) {
        if (names[name]) {
            delete names[name];
        }
    };

    return {
        claim: claim,
        free: free,
        get: get,
        getGuestName: getGuestName
    };
}());

//server collector
var collector = {};

io.sockets.on('connection', function (socket) {
    var user = userNames.getGuestName();
    console.log("CUSTOM LOG:\tsocket.io connection made, user: " + user);

    /********
     configuration
     ********/

    var model_folder = "models/";
    var api = {};

    /***** PACKAGE CONSTRUCTORS *****/

    api.pack = {};

    /* Basics */

    //simple JSON (UUID not relevant)
    api.pack.simple = function(channel, data) {
        return {
            "channel" : channel,
            "data" : data
        };
    };
    api.pack.nopack = function(data) {
        return data;
    };
    //JSON with uuid, all data in package 'data'
    api.pack.uuid = function(channel, uuid, data) {
        return {
            "channel" : channel,
            "data" : {
                "uuid": uuid,
                "data": data
            }
        };
    };

    /* For API */

    //pack for clientAPI
    api.pack.api_wrap = function(channel, data) {
        var packaged = api.pack.simple(channel, data);
        return JSON.stringify(packaged);
    };
    //pack for channel $clotho, for Clotho ServerAPI
    api.pack.clotho = function(uuid, data) {
        var packaged = api.pack.uuid('$clotho', uuid, data);
        return JSON.stringify(packaged);
    };

    /* Command Specifics */

    //pack for a collect command
    api.pack.collect = function(uuid, type, model_or_url, isURL) {
        isURL = !!isURL || false; //for templates
        return {
            "uuid" : uuid,
            "type" : type,
            "model" : model_or_url,
            "isURL" : isURL
        }
    };
    //pack for a say or alert command
    api.pack.message = function(msg) {
        return {
            "msg" : msg
        }
    };
    api.pack.display = function(uuid, data) {
        return {
            "uuid" : uuid,
            "data" : data
        };
    };

    /***** SIMPLE MESSAGES *****/

    api.api = {};
    api.api.alert = function(data) {
        var user = data.userID;
        //todo - logic to route to a specific user

        var msg = data.msg;

        socket.send(api.pack.api_wrap('alert',
            api.pack.nopack(msg)
        ));
    };
    api.api.broadcast = function(data) {
        var channel = data.channel;
        data = data.data;

        socket.send(api.pack.api_wrap('broadcast',
            api.pack.simple(channel, data)
        ));
    };
    api.api.get = function(uuid) {
        console.log("requesting model: " + uuid);
        var path = require('path').resolve(model_folder, uuid + ".json");

        //future - send templates also, not just JSON filetypes

        fs.readFile(path, 'utf8', function (err, data) {
            if (err) { console.log('Error: ' + err); return; }
            data = JSON.parse(data);
            //console.log("CUSTOM:\tgot the file!");

            //make assumption that get will be in client collector, so put in local collector so don't re-send on set()
            collector[uuid] = data;

            socket.send(api.pack.api_wrap('collect',
                api.pack.collect(uuid, "json", data)
            ));
        });
    };
    api.api.get_script = function(uuid) {
        console.log("requesting script for: " + uuid);
        var path = 'partials/' + uuid + '.js';

        socket.send(api.pack.api_wrap('collect',
            api.pack.collect(uuid, "js", path, true)
        ));
    };
    api.api.get_template = function(uuid) {
        console.log("requesting partial URL for: " + uuid);
        var path = 'partials/' + uuid + '.html';

        socket.send(api.pack.api_wrap('collect',
            api.pack.collect(uuid, "html", path, true)
        ));
    };
    api.api.get_url = function(uuid) {
        console.log("requesting model URL for: " + uuid);
        var path = 'models/' + uuid;

        socket.send(api.pack.api_wrap('collect',
            api.pack.collect(uuid, "json", path)
        ));
    };
    api.api.gradeQuiz = function (quiz) {
        console.log(quiz);

        //just send back everything is correct for now...
        var result = {
            "uuid" : quiz.uuid,
            "response" : quiz.answer,
            "result" : true
        };

        socket.send(api.pack.api_wrap('quizResult:' + quiz.uuid,
            api.pack.nopack(result)
        ));
    };
    api.api.log = function(data) {
        console.log("LOG\t: " + data);
    };
    api.api.notify = function(data) {
        console.log("notification!");
        console.log(data);
    };
    api.api.requestRecent = function() {

        var path = require('path').resolve(model_folder, 'recent.json');

        fs.readFile(path, 'utf8', function (err, data) {
            if (err) { console.log('Error: ' + err); return; }
            data = JSON.parse(data);

            socket.send(api.pack.api_wrap('displayRecent',
                api.pack.nopack(data)
            ));
        });

    };
    api.api.say = function(data) {
        var user = data.userID,
        //todo - separation for sending messages from "server" vs. "client"
        sender = data.sender || "server",
        msg = data.msg,
        timestamp = data.timestamp || Date.now();
        //todo - add classes: error, warning, success, muted, info
        var css = data.css || "info";

        var message = {
            "text" : msg,
            "from" : sender,
            "class" : css,
            "timestamp" : timestamp
        };

        socket.send(api.pack.api_wrap('say',
            api.pack.nopack(message)
        ));
    };
    api.api.set = function(data) {
        var uuid = data.uuid;
        data = data.data;

        //only save and send data to client if different from what the server has
        if (!equals(data, collector[uuid])) {
            console.log("different! storing " + uuid);

            // save to server collector
            collector[uuid] = data;

            socket.send(api.pack.api_wrap('collect',
                api.pack.collect(uuid, "json", data)
            ));
        }
    };
    api.api.show = function(data) {
        var uuid = data.uuid;
        data = data.data;

        //for testing, just send a dummy model and template
        data = {
            "model" : "inst_second",
            "template" : "show_template"
        };

        socket.send(api.pack.api_wrap('display',
            api.pack.display(uuid, data)
        ));
    };
    api.api.show_simple = function(data) {

        /*
        //demo
        var message = {
            "template" : "extensions/simple-template.html",
            "target" : "body"
            "controller" : "extensions/simple-controller.js",
            "dependencies" : [
                "extensions/simple-service.js"
            ],
            styles : {
                "background-color" : "#FF0000"
            }
        };
        */

        socket.send(api.pack.api_wrap('display_simple',
            api.pack.nopack(data)
        ));
    };

    /**** SEARCHBAR ***/

    api.searchbar = {};

    api.searchbar.submit = function (data) {


        var message = {
            "text" : data.query,
            "from" : "client",
            "class" : "muted",
            "timestamp" : Date.now()
        };

        socket.send(api.pack.api_wrap('say',
            api.pack.nopack(message)
        ));

        // more serverside logic would happen here
        console.log("submit received: " + data.query);

        var response = {
            "text" : 'Demo Server response',
            "from" : "server",
            "class" : "info",
            "timestamp" : Date.now()
        };

        socket.send(api.pack.api_wrap('say',
            api.pack.nopack(response)
        ));

    };

    api.searchbar.autocomplete = function(data) {
        var query = data.query;

        var demo = [
            {
                "uuid" : "1234567890",
                "text" : "This is a command",
                "value" : "clotho.run('rlaksd', ['sadfjsklsdf']);",
                "command" : {
                    "channel" : "run",
                    "functionId" : "rlaksd",
                    "args" : "['sadfjsklsdf']"
                },
                "type" : "command"
            },
            {
                "uuid" : "qwertyuiop",
                "text" : "This is a phrase",
                "type" : "phrase"
            },
            {
                "uuid" : "817924532",
                "text" : "Reverse Complement pca1502",
                "value" : "clotho.run('sdfsadg', '23tg3e2q');",
                "command" : {
                    "channel" : "run",
                    "functionId" : "sdfsadg",
                    "args" : "'23tg3e2q'"
                },
                "type" : "command"
            },
            {
                "uuid" : "xxxxxxxxxxxx",
                "text" : "Reverse complement this",
                "value" : "Reverse complement this",
                "type" : "phrase"
            }
        ];

        socket.send(api.pack.api_wrap('autocomplete',
            api.pack.nopack(demo)
        ));
    };
    api.searchbar.autocompleteDetail = function(data) {
        var uuid = data.uuid;
        console.log("requested detail for uuid: " + data.uuid);

        var path = require('path').resolve(model_folder, 'detail_' + uuid + ".json");

        fs.readFile(path, 'utf8', function (err, data) {
            if (err) { console.log('Error: ' + err); return; }
            data = JSON.parse(data);

            socket.send(api.pack.api_wrap('autocompleteDetail',
                api.pack.nopack(data)
            ));
        });
    };


    /*************
     * CUSTOM EVENTS
     *************/

    /***** CHAT *****/

    api.chat_start = function() {
        socket.send(api.pack.api_wrap('broadcast',
            api.pack.simple('chat:init', {name: user, users: userNames.get(), messages: []} )
        ));
    };
    //for receiving messages
    api.chat_send = function (msg, user) {
        console.log("received chat message, broadcasting: " + msg);

        //send to self
        socket.send(api.pack.api_wrap('broadcast',
            api.pack.uuid('chat:receive', user, msg)
        ));
        //send to others
        socket.broadcast.send(api.pack.api_wrap('broadcast',
            api.pack.uuid('chat:receive', user, msg)
        ));
    };


    /*********
    basic events / simple connection management
    - uses some socket.io specifics, but that'll be replaced later with Clotho specifics anyways
        - broadcast tag is the main one...
     ********/


    socket.broadcast.send(api.pack.api_wrap('broadcast',
        api.pack.uuid('chat:userJoin', user, '')
    ));

    socket.on('disconnect', function () {
        io.sockets.emit('user disconnected');
        userNames.free(user);
        socket.broadcast.send(api.pack.api_wrap('broadcast',
            api.pack.uuid('chat:userLeave', user, '')
        ));
    });


    /********
     custom events
     ********/

    /* wrapped like so:

    {
        channel : < 'api' | 'searchbar' >
        data : {
            channel : <channel>
            data : {object}
        }
    }

     */


    socket.on('message', function(msg) {
        msg = JSON.parse(msg);
        var wrapper = msg.channel;
        var channel = msg.data.channel;
        var data = msg.data.data;

        console.log("CUSTOM\t command in " + wrapper + " received: " + channel + "\tdata: " + data);
        
        if (api[wrapper][channel]) {
            //easy....
            api[wrapper][channel](data);
        }
        else if (api[channel]){
            //custom events
            api[channel](data);
        }
        else {
            console.log("!!!\tNo match... channel: " + channel);
        }
    });

});

// Must be last:
main(process.argv);
