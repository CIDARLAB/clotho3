'use strict';

// see http://www.html5rocks.com/en/tutorials/frameworks/angular-websockets/

Application.Chat.service('ChatSocket', ['Clotho', function(Clotho) {

    var users, messages;
    //need to return object to be reference, not copied
    users = {"data" : []};
    messages = {"data" : []};
    var usr_resolved;

    Clotho.emit("chat:start");

    Clotho.listen('chat:init', function(data) {
        console.log("CHATSOCKET\tChat Initialized");
        usr_resolved = data.name;
        users.data = data.users;
        messages.data = data.messages;
    }, 'ChatSocket');

    Clotho.listen('chat:userJoin', function(user) {
        users.data.push(user.uuid);
        messages.data.push({
            uuid: 'CHATROOM',
            data: 'User ' + user.uuid + ' has joined.'
        });
    }, 'ChatSocket');

    Clotho.listen('chat:userLeave', function(user) {
        var i, usr_tmp;
        for (i = 0; i < users.data.length; i++) {
            usr_tmp = users.data[i];
            if (usr_tmp === user.uuid) {
                users.data.splice(i, 1);
                break;
            }
        }
        messages.data.push({
            uuid: 'CHATROOM',
            data: 'User ' + user.uuid + ' has left.'
        });
    }, 'ChatSocket');

    Clotho.listen('chat:receive', function(chatData) {
        messages.data.push(chatData);
    }, 'ChatSocket');

    /*** simple ***/

    //send a message over the socket with a custom channel and JSON object
    var emit = function(channel, data) {
        Clotho.emit(channel, data);
    };

    //request a model by uuid
    var requestModel = function(uuid) {
        Clotho.get(uuid);
    };

    //request a model by uuid
    var requestModelURL = function(uuid) {
        Clotho.emit('get:url', uuid);
    };

    /**** chat ****/

    //send a message for chat
    var chat_send = function(msg) {
        //future - do escaping / error checking etc.
        var data = '{"msg" : "' + msg + '", "username" : "' + usr_resolved + '"}';
        Clotho.emit('chat:send', JSON.parse(data) );
    };

    return {
        users : users,
        messages : messages,
        emit : emit,
        requestModel : requestModel,
        requestModelURL : requestModelURL,
        chat_send : chat_send
    }
}]);