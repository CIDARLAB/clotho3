#!/usr/bin/env node

var prettyjson = require('prettyjson');
var stomp = require('stomp');

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

var client = new stomp.Stomp(stomp_args);

var queue = '/queue/CLOTHO';

client.connect();

function sleep(milliSeconds) {
    var startTime = new Date().getTime();
    while (new Date().getTime() < startTime + milliSeconds);
}

client.on('connected', function() {
    num = 10;
    for (var i = 0; i < num; i++) {
        console.log("sending message...");
        
        var body = {
    		channel: 'EXECUTION',
    		action: 'SHOW',
            data: 
            	{
        			id: 'node-client'
           		}        
        };
        
        client.send({
            'destination': queue,
            'persistent': 'false',
        	'request': prettyjson.render(JSON.stringify(body))
        }, receipt);

        sleep(250);
    }
    console.log('Produced ' + num + ' messages');
});

client.on('receipt', function(receipt) {
    
    console.log("RECEIPT: " + receipt);
});

client.on('error', function(error_frame) {
    console.log(error_frame.body);
    client.disconnect();
});

process.on('SIGINT', function() {
    console.log('Produced ' + num + ' messages');
    client.disconnect();
    process.exit(0);
});
