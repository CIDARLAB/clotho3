"use strict";

angular.module('clotho.core')
/**
 * @name Communicator
 *
 * @description
 * Service which handles communication to server, delegating to Socket if available, otherwise relying on REST.
 *
 * Handles Authentication for each.
 *
 * todo - handle authentication
 *
 * todo - handle communication FROM server once support SSE (currently logic in socket service) -- avoid circular deps
 *
 */
.service('Communicator', function (ClothoAuth, Socket, ClothoREST, Debug) {

  var Debugger = new Debug('Communicator', '#bbbb55'),
      strToJson = angular.fromJson, //strip $$-prefixed
      jsonToStr = angular.toJson,
      websocketMaxLength = 16348;

  //check if a message (string) can be sent over the socket
  function validSocketMessage (msgString) {
    return msgString.length < websocketMaxLength;
  }

  //determine whether message (string) should go over SOCKET or REST
  //todo update once REST is working
  var canUseSocket = function communicatorCanUseSocket (msgString) {
    //check if socket is ready
    //todo - smarter check for readyState, because queueing ok
    if (Socket.state() === 1) {
      return true;
    }
    if (!(validSocketMessage(msgString))) {
      return false;
    }
    //default
    return true;
  };

    /**
     * @name Communicator.send
     * @description
     * Given a properly formed message to sent to the server, delegate to Socket or REST and send.
     *
     * @param msgObject
     * @param forceRest
     *
     * @returns {Boolean} Whether message successfully sent
     */
  var send = function (msgObject, forceRest) {
    if (forceRest === true) {
      return ClothoREST.send(msgObject);
    }

    var msgString = angular.isObject(msgObject) ? strToJson(msgObject) : msgObject;
    if (canUseSocket(msgString)) {
      return Socket.send(msgString);
    } else {
      Debugger.warn('Message must be sent over rest', msgObject);
      return ClothoREST.send(msgObject);
    }
  };

    /**
     * @name Communicator.emit
     * @description
     * Sends a message on an channel with data passed as arguments. Can force message to go over rest if you want.
     *
     * Note that unsupported channels may crash the socket.
     *
     * @param channel {String} Clotho API channel (or custom)
     * @param data {Object} Clotho API data
     * @param options {Object} Clotho API options
     * @param callback {Function} callback executed directly after sent. Passed `wasSent` {Boolean} and `packaged` {Object} object sent to server. Listen to server-sent or PubSub events to get a real callback on data return.
     * @param forceRest {Boolean} `true` to force REST instead of socket
     *
     * @returns {Boolean} whether message sent
     */
  var emit = function (channel, data, options, callback, forceRest) {
    var packaged = {
      "channel" : channel,
      "data" : data,
      "options" : options
    };

    var wasSent = send(packaged, forceRest);

    if (wasSent && angular.isFunction(callback)) {
      callback(wasSent, packaged);
    }

    return wasSent;
  };

  return  {
    send: send,
    emit : emit
  };
});
