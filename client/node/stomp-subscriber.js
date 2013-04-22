#!/usr/bin/env node

var sys = require('util');
var stomp = require('stomp');
var prettyjson = require('prettyjson');


// Set to true if you want a receipt
// of all messages sent.
var receipt = true;

// Set debug to true for more verbose output.
// login and passcode are optional (required by rabbitMQ)
var stomp_args = {
    port: 61613,
    host: 'localhost',
    login: 'admin',
    passcode: 'password'
}

var headers = {
    destination: '/topic/CLOTHO.UPDATE',
    ack: 'client',
    'activemq.prefetchSize': '10'
};

var client = new stomp.Stomp(stomp_args);
var messages = 0;

client.connect();

function message_callback(body, headers) {
    console.log('Message Callback Fired!');
    console.log('Headers: ' + sys.inspect(headers));
    console.log('Body: ' + body);
}

client.on('connected', function() {
    client.subscribe(headers, message_callback);
    console.log('Connected');
});

client.on('message', function(message) {
    //console.log("HEADERS: " + sys.inspect(message.headers));
    //console.log("BODY: " + message.body);
    console.log("Got message: " + message.headers['message-id']);
    client.ack(message.headers['message-id']);
    messages++;
});


client.on('error', function(error_frame) {
    console.log(error_frame.body);
    client.disconnect();
});

process.on('SIGINT', function() {
    client.disconnect();
    process.exit(0);
});
